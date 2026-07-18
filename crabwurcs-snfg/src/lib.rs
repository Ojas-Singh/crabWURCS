use crabwurcs_core::{Monosaccharide, ResidueGraph};
use petgraph::graph::NodeIndex;
use petgraph::visit::EdgeRef;
use petgraph::Direction;
use std::collections::HashMap;
use thiserror::Error;

// ── Error types ────────────────────────────────────────────────────────────

#[derive(Debug, Error)]
pub enum SnfgError {
    #[error(
        "no SNFG symbol for monosaccharide: backbone_len={}, code={}",
        .0.backbone_length,
        .0.skeleton_code
    )]
    UnknownSymbol(Monosaccharide),

    #[error(transparent)]
    Core(#[from] crabwurcs_core::CoreError),
}

pub type SnfgResult<T> = Result<T, SnfgError>;

// ── SNFG shapes ────────────────────────────────────────────────────────────

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Shape {
    Circle,
    Square,
    NSquare,
    Triangle,
    Diamond,
    SplitDiamondTop,    // uronic acid – top half coloured
    SplitDiamondBottom, // uronic acid – bottom half coloured (IdoA)
    Star,
    Hexagon,
    Pentagon,
}

// ── SNFG colours (glycoshape.io / standard SNFG palette) ───────────────────

pub mod colour {
    pub const GLC: &str = "#0385AE"; // dark blue – Glc, GlcNAc
    pub const GAL: &str = "#FCC326"; // golden yellow – Gal, GalNAc
    pub const MAN: &str = "#058F60"; // dark green – Man, ManNAc
    pub const FUC: &str = "#C23537"; // dark red – Fuc
    pub const NEU5AC: &str = "#A15989"; // GlycoShape mauve – Neu5Ac
    pub const NEU5GC: &str = "#91D3E3"; // GlycoShape light blue – Neu5Gc
    pub const KDN: &str = "#5995B3"; // slate blue – KDN
    pub const IDOA: &str = "#9F6D55"; // tan/brown – IdoA
    pub const KDO: &str = "#FCC326"; // yellow – KDO
    pub const XYL: &str = "#F47421"; // orange – Xyl
    pub const GUL: &str = "#EF6130"; // orange – Gul (GlycoShape palette)
    pub const ALT: &str = "#F69EA1"; // pink – Alt (official SNFG RGB)
    pub const ALL: &str = "#A54399"; // purple – All (SNFG RGB palette)
    pub const TAL: &str = "#91D3E3"; // light blue – Tal (GlycoShape palette)
    pub const IDO: &str = "#A17A4D"; // brown – Ido (official SNFG RGB)
    pub const UNKNOWN: &str = "#999999";

    pub const STROKE: &str = "#000000";
    pub const BOND: &str = "#000000";
    pub const LINKAGE_TEXT: &str = "#000000";
}

#[derive(Debug, Clone, Copy)]
pub struct Symbol {
    pub shape: Shape,
    pub fill: &'static str,
    pub label: &'static str,
}

#[allow(clippy::too_many_arguments)]
fn hexose_family_symbol(
    fill: &'static str,
    neutral: &'static str,
    acid: &'static str,
    nac: &'static str,
    amine: &'static str,
    has_acid: bool,
    has_n_mod: bool,
    has_nac: bool,
    acid_bottom: bool,
) -> Symbol {
    if has_acid {
        return Symbol {
            shape: if acid_bottom {
                Shape::SplitDiamondBottom
            } else {
                Shape::SplitDiamondTop
            },
            fill,
            label: acid,
        };
    }
    if has_n_mod {
        return Symbol {
            shape: if has_nac {
                Shape::Square
            } else {
                Shape::NSquare
            },
            fill,
            label: if has_nac { nac } else { amine },
        };
    }
    Symbol {
        shape: Shape::Circle,
        fill,
        label: neutral,
    }
}

// ── Monosaccharide → SNFG symbol ───────────────────────────────────────────

pub fn symbol_for(residue: &Monosaccharide) -> SnfgResult<Symbol> {
    let skel = &residue.skeleton_code;
    let bare: String = skel.chars().take_while(|c| c.is_ascii_digit()).collect();
    let bare_str: &str = bare.as_str();

    let has_nac = residue
        .modifications
        .iter()
        .any(|m| m.descriptor.contains("NCC") || m.descriptor.contains("NC"));

    let has_n_mod = residue.modifications.iter().any(|m| {
        m.descriptor.contains("NCC") || m.descriptor.contains("NC") || m.descriptor.starts_with('N')
    });

    let has_ngc = residue
        .modifications
        .iter()
        .any(|m| m.descriptor.contains("NCCO") || m.descriptor.contains("NCO"));
    let has_deoxy = residue.skeleton_code.ends_with('m') || skel.contains('d');
    let has_acid = skel.contains('A');

    // Ketohexoses have only three stereogenic backbone digits, so they must
    // be classified before the aldopentose fallback below.  In particular,
    // fructofuranose is the green SNFG pentagon used by GlycoShape and
    // glycowork (the ring form is chemically meaningful, not decoration).
    if residue.anomeric_prefix.starts_with('h') && bare_str == "122" {
        return Ok(Symbol {
            shape: Shape::Pentagon,
            fill: colour::MAN,
            label: "Fru",
        });
    }

    // WURCS encodes ulosonic-acid oxidation in the leading `A*` carbon
    // descriptors rather than in the trailing skeleton string. KDO and Bac
    // belong to SNFG's flat-hexagon "unknown/other" family.
    if bare_str == "1122" && residue.anomeric_prefix.starts_with('A') {
        return Ok(Symbol {
            shape: Shape::Hexagon,
            fill: colour::KDO,
            label: "KDO",
        });
    }
    let n_positions = residue
        .modifications
        .iter()
        .filter(|modification| modification.descriptor.contains("NCC"))
        .map(|modification| modification.position.0)
        .collect::<std::collections::HashSet<_>>();
    if bare_str == "2122" && has_deoxy && n_positions.contains(&2) && n_positions.contains(&4) {
        return Ok(Symbol {
            shape: Shape::Hexagon,
            fill: colour::GLC,
            label: "Bac",
        });
    }

    // Each pair contains D/L-inverted SkeletonCodes. SNFG keeps the family
    // colour; non-default absolute configuration belongs in the figure legend.
    match bare_str {
        "2122" | "1211" => {
            return Ok(hexose_family_symbol(
                colour::GLC,
                "Glc",
                "GlcA",
                "GlcNAc",
                "GlcN",
                has_acid,
                has_n_mod,
                has_nac,
                false,
            ))
        }
        "2112" | "1221" if !has_deoxy => {
            return Ok(hexose_family_symbol(
                colour::GAL,
                "Gal",
                "GalA",
                "GalNAc",
                "GalN",
                has_acid,
                has_n_mod,
                has_nac,
                false,
            ))
        }
        "1122" | "2211" if !has_deoxy => {
            return Ok(hexose_family_symbol(
                colour::MAN,
                "Man",
                "ManA",
                "ManNAc",
                "ManN",
                has_acid,
                has_n_mod,
                has_nac,
                false,
            ))
        }
        "2212" | "1121" => {
            return Ok(hexose_family_symbol(
                colour::GUL,
                "Gul",
                "GulA",
                "GulNAc",
                "GulN",
                has_acid,
                has_n_mod,
                has_nac,
                false,
            ))
        }
        "1222" | "2111" => {
            return Ok(hexose_family_symbol(
                colour::ALT,
                "Alt",
                "AltA",
                "AltNAc",
                "AltN",
                has_acid,
                has_n_mod,
                has_nac,
                false,
            ))
        }
        "2222" | "1111" => {
            return Ok(hexose_family_symbol(
                colour::ALL,
                "All",
                "AllA",
                "AllNAc",
                "AllN",
                has_acid,
                has_n_mod,
                has_nac,
                false,
            ))
        }
        "1112" | "2221" => {
            return Ok(hexose_family_symbol(
                colour::TAL,
                "Tal",
                "TalA",
                "TalNAc",
                "TalN",
                has_acid,
                has_n_mod,
                has_nac,
                false,
            ))
        }
        "2121" | "1212" => {
            return Ok(hexose_family_symbol(
                colour::IDO,
                "Ido",
                "IdoA",
                "IdoNAc",
                "IdoN",
                has_acid,
                has_n_mod,
                has_nac,
                true,
            ))
        }
        _ => {}
    }

    // ── Hexoses ────────────────────────────────────────────────────────
    match bare_str {
        "2122" if has_acid => {
            return Ok(Symbol {
                shape: Shape::SplitDiamondTop,
                fill: colour::GLC,
                label: "GlcA",
            });
        }
        "2121" if has_acid => {
            return Ok(Symbol {
                shape: Shape::SplitDiamondBottom,
                fill: colour::IDOA,
                label: "IdoA",
            });
        }
        "2122" | "2121" if has_n_mod => {
            if has_nac {
                return Ok(Symbol {
                    shape: Shape::Square,
                    fill: colour::GLC,
                    label: "GlcNAc",
                });
            } else {
                return Ok(Symbol {
                    shape: Shape::NSquare,
                    fill: colour::GLC,
                    label: "GlcN",
                });
            }
        }
        "2122" | "2121" => {
            return Ok(Symbol {
                shape: Shape::Circle,
                fill: colour::GLC,
                label: "Glc",
            });
        }
        "2112" | "2111" if has_acid => {
            return Ok(Symbol {
                shape: Shape::SplitDiamondTop,
                fill: colour::GAL,
                label: "GalA",
            });
        }
        "2112" | "2111" if has_n_mod => {
            if has_nac {
                return Ok(Symbol {
                    shape: Shape::Square,
                    fill: colour::GAL,
                    label: "GalNAc",
                });
            } else {
                return Ok(Symbol {
                    shape: Shape::NSquare,
                    fill: colour::GAL,
                    label: "GalN",
                });
            }
        }
        "2112" | "2111" => {
            return Ok(Symbol {
                shape: Shape::Circle,
                fill: colour::GAL,
                label: "Gal",
            });
        }
        "1221" if has_acid => {
            return Ok(Symbol {
                shape: Shape::SplitDiamondTop,
                fill: colour::MAN,
                label: "ManA",
            });
        }
        "1221" if has_deoxy => {
            return Ok(Symbol {
                shape: Shape::Triangle,
                fill: colour::FUC,
                label: "Fuc",
            });
        }
        "1221" if has_n_mod => {
            if has_nac {
                return Ok(Symbol {
                    shape: Shape::Square,
                    fill: colour::MAN,
                    label: "ManNAc",
                });
            } else {
                return Ok(Symbol {
                    shape: Shape::NSquare,
                    fill: colour::MAN,
                    label: "ManN",
                });
            }
        }
        "1221" => {
            return Ok(Symbol {
                shape: Shape::Circle,
                fill: colour::MAN,
                label: "Man",
            });
        }
        "1122" if has_acid => {
            return Ok(Symbol {
                shape: Shape::SplitDiamondTop,
                fill: colour::MAN,
                label: "ManA",
            });
        }
        "1122" if has_n_mod => {
            if has_nac {
                return Ok(Symbol {
                    shape: Shape::Square,
                    fill: colour::MAN,
                    label: "ManNAc",
                });
            } else {
                return Ok(Symbol {
                    shape: Shape::NSquare,
                    fill: colour::MAN,
                    label: "ManN",
                });
            }
        }
        "1122" => {
            return Ok(Symbol {
                shape: Shape::Circle,
                fill: colour::MAN,
                label: "Man",
            });
        }
        "2211" if has_acid => {
            return Ok(Symbol {
                shape: Shape::SplitDiamondTop,
                fill: colour::MAN,
                label: "ManA",
            });
        }
        "2211" if has_deoxy => {
            return Ok(Symbol {
                shape: Shape::Triangle,
                fill: colour::FUC,
                label: "Rha",
            });
        }
        "2211" if has_n_mod => {
            if has_nac {
                return Ok(Symbol {
                    shape: Shape::Square,
                    fill: colour::MAN,
                    label: "ManNAc",
                });
            } else {
                return Ok(Symbol {
                    shape: Shape::NSquare,
                    fill: colour::MAN,
                    label: "ManN",
                });
            }
        }
        "2211" => {
            return Ok(Symbol {
                shape: Shape::Circle,
                fill: colour::MAN,
                label: "Man",
            });
        }
        "1121" => {
            return Ok(Symbol {
                shape: Shape::Circle,
                fill: colour::GUL,
                label: "Gul",
            });
        }
        "2222" => {
            return Ok(Symbol {
                shape: Shape::Circle,
                fill: colour::ALL,
                label: "All",
            });
        }
        "2221" => {
            return Ok(Symbol {
                shape: Shape::Circle,
                fill: colour::TAL,
                label: "Tal",
            });
        }
        _ if bare_str.contains('d') && bare_str.len() <= 5 && has_deoxy => {
            return Ok(Symbol {
                shape: Shape::Triangle,
                fill: colour::FUC,
                label: "Fuc",
            });
        }
        _ => {}
    }

    // ── Sialic acids (9‑carbon backbones) ──────────────────────────────
    if bare_str.len() >= 5 && (bare_str.contains("21122") || bare_str.contains("11212")) {
        if has_ngc {
            return Ok(Symbol {
                shape: Shape::Diamond,
                fill: colour::NEU5GC,
                label: "Neu5Gc",
            });
        }
        if has_nac {
            return Ok(Symbol {
                shape: Shape::Diamond,
                fill: colour::NEU5AC,
                label: "Neu5Ac",
            });
        }
        if residue.modifications.iter().any(|m| m.descriptor == "O") {
            return Ok(Symbol {
                shape: Shape::Diamond,
                fill: colour::KDN,
                label: "KDN",
            });
        }
        return Ok(Symbol {
            shape: Shape::Diamond,
            fill: colour::NEU5AC,
            label: "Sia",
        });
    }

    // ── KDO ────────────────────────────────────────────────────────────
    if bare_str.len() == 4 && bare_str.contains("1122") && has_acid {
        return Ok(Symbol {
            shape: Shape::Diamond,
            fill: colour::KDO,
            label: "KDO",
        });
    }

    // ── Pentoses ───────────────────────────────────────────────────────
    if bare_str.len() == 3 {
        return Ok(Symbol {
            shape: Shape::Star,
            fill: colour::XYL,
            label: "Xyl",
        });
    }

    // Composition WURCS commonly uses `xxxx` when stereochemistry is not
    // specified. The chemical class is still known from the backbone and
    // N-acetyl substituent, so use neutral Hex/HexNAc SNFG symbols rather
    // than a star (which denotes a pentose family).
    if skel.contains("xxxx") {
        return Ok(Symbol {
            shape: if has_n_mod {
                Shape::Square
            } else {
                Shape::Circle
            },
            fill: colour::UNKNOWN,
            label: if has_nac { "HexNAc" } else { "Hex" },
        });
    }

    // ── Fallback ───────────────────────────────────────────────────────
    let (shape, fill) = match bare_str.len() {
        5 => (Shape::Hexagon, colour::UNKNOWN),
        4 => (Shape::Pentagon, colour::UNKNOWN),
        _ => (Shape::Star, colour::UNKNOWN),
    };
    Ok(Symbol {
        shape,
        fill,
        label: "?",
    })
}

// ── Geometry constants ─────────────────────────────────────────────────────

pub const NODE_R: f64 = 25.0;
pub const H_SPACING: f64 = 100.0;
pub const V_SPACING: f64 = 100.0;
pub const BOND_W: f64 = 4.0;
pub const LABEL_SIZE: f64 = 20.0;

// ── Options ────────────────────────────────────────────────────────────────

#[derive(Debug, Clone)]
pub struct RenderOptions {
    pub colour: bool,
    pub show_labels: bool,   // show residue abbreviations inside shapes
    pub show_linkages: bool, // show linkage positions on bonds
    pub font_family: String,
    pub scale: f64,
}

impl Default for RenderOptions {
    fn default() -> Self {
        Self {
            colour: true,
            show_labels: false, // SNFG convention: shape + colour = identity
            show_linkages: true,
            font_family: "Arial, Helvetica, sans-serif".into(),
            scale: 1.0,
        }
    }
}

// ── Tree layout ────────────────────────────────────────────────────────────

#[derive(Debug, Clone, Default)]
struct LayoutInfo {
    x: f64,
    y: f64,
}

/// Recursive tree layout: post‑order, each node centred among children.
/// Children fan out vertically from the parent's y with `V_SPACING` separation.
fn compute_layout(graph: &ResidueGraph, root: NodeIndex) -> HashMap<usize, LayoutInfo> {
    let mut info = HashMap::new();
    let mut visited = std::collections::HashSet::new();
    let mut next_leaf = 0usize;
    layout_subtree(graph, root, 0, &mut next_leaf, &mut info, &mut visited);
    // WURCS compositions and undefined antennae can contain disconnected
    // components. Lay every component out instead of silently omitting it.
    for node in graph.inner().node_indices() {
        if !visited.contains(&node.index()) {
            if next_leaf > 0 {
                next_leaf += 1;
            }
            layout_subtree(graph, node, 0, &mut next_leaf, &mut info, &mut visited);
        }
    }

    resolve_fucose_collisions(graph, &mut info);

    // centre around y=0
    let min_y = info.values().map(|li| li.y).fold(f64::MAX, f64::min);
    let max_y = info.values().map(|li| li.y).fold(f64::MIN, f64::max);
    let shift = -(min_y + max_y) / 2.0;
    for li in info.values_mut() {
        li.y += shift;
    }
    info
}

fn resolve_fucose_collisions(graph: &ResidueGraph, info: &mut HashMap<usize, LayoutInfo>) {
    let mut branches = graph
        .inner()
        .edge_references()
        .filter(|edge| is_fucose(graph, edge.target()) && is_terminal(graph, edge.target()))
        .map(|edge| (edge.source(), edge.target()))
        .collect::<Vec<_>>();
    branches.sort_by(|(left_parent, _), (right_parent, _)| {
        info[&left_parent.index()]
            .y
            .total_cmp(&info[&right_parent.index()].y)
    });

    for (parent, fucose) in branches {
        let parent_layout = info[&parent.index()].clone();
        let desired_y = parent_layout.y + V_SPACING;
        let collision = info.iter().any(|(index, layout)| {
            *index != fucose.index()
                && (layout.x - parent_layout.x).abs() < f64::EPSILON
                && (layout.y - desired_y).abs() < f64::EPSILON
        });
        if collision {
            for (index, layout) in info.iter_mut() {
                if *index != fucose.index() && layout.y >= desired_y {
                    layout.y += V_SPACING;
                }
            }
        }
        info.insert(
            fucose.index(),
            LayoutInfo {
                x: parent_layout.x,
                y: desired_y,
            },
        );
    }
}

fn layout_subtree(
    graph: &ResidueGraph,
    node: NodeIndex,
    depth: usize,
    next_leaf: &mut usize,
    info: &mut HashMap<usize, LayoutInfo>,
    visited: &mut std::collections::HashSet<usize>,
) -> f64 {
    if !visited.insert(node.index()) {
        return info.get(&node.index()).map(|li| li.y).unwrap_or(0.0);
    }

    let mut children: Vec<(NodeIndex, u8)> = graph
        .inner()
        .edges_directed(node, Direction::Outgoing)
        .filter(|edge| edge.weight().repeat.is_none() && !edge.weight().cyclic)
        .map(|edge| (edge.target(), edge.weight().parent_position.0))
        .filter(|(child, _)| !visited.contains(&child.index()))
        .collect();
    // Standard SNFG branch order places higher acceptor positions above
    // lower ones in the conventional right-to-left layout (for example the
    // N-glycan α1-6 arm above α1-3, and β1-4 above β1-2).
    children.sort_by_key(|(child, position)| (std::cmp::Reverse(*position), child.index()));

    let y = if children.is_empty() {
        let y = *next_leaf as f64 * V_SPACING;
        *next_leaf += 1;
        y
    } else {
        let (fucose_children, ordinary_children): (Vec<_>, Vec<_>) = children
            .into_iter()
            .partition(|(child, _)| is_fucose(graph, *child) && is_terminal(graph, *child));

        let y = if ordinary_children.is_empty() {
            let y = *next_leaf as f64 * V_SPACING;
            *next_leaf += 1;
            y
        } else {
            let child_y = ordinary_children
                .into_iter()
                .map(|(child, _)| layout_subtree(graph, child, depth + 1, next_leaf, info, visited))
                .collect::<Vec<_>>();
            (child_y[0] + child_y[child_y.len() - 1]) / 2.0
        };

        // SNFG convention draws terminal fucose vertically below its parent.
        // It remains at the parent's horizontal depth but must never share the
        // parent's coordinates (the old behaviour overprinted the triangle on
        // Gal/GlcNAc in examples such as GS00742 and GS00169).
        for (index, (child, _)) in fucose_children.into_iter().enumerate() {
            visited.insert(child.index());
            info.insert(
                child.index(),
                LayoutInfo {
                    x: depth as f64,
                    y: y + (index + 1) as f64 * V_SPACING,
                },
            );
        }
        y
    };
    info.insert(node.index(), LayoutInfo { x: depth as f64, y });
    y
}

fn is_fucose(graph: &ResidueGraph, node: NodeIndex) -> bool {
    graph
        .residue(node)
        .and_then(|residue| symbol_for(residue).ok())
        .is_some_and(|symbol| symbol.shape == Shape::Triangle)
}

fn is_terminal(graph: &ResidueGraph, node: NodeIndex) -> bool {
    graph
        .inner()
        .edges_directed(node, Direction::Outgoing)
        .all(|edge| edge.weight().repeat.is_some() || edge.weight().cyclic)
}

// ── Linkage label ──────────────────────────────────────────────────────────

fn anomer_char(anom: crabwurcs_core::AnomericSymbol, prefix: &str) -> &'static str {
    let c = anom.to_char();
    if c != 'x' {
        return match c {
            'a' => "\u{03B1}",
            'b' => "\u{03B2}",
            'o' => "o",
            _ => "?",
        };
    }
    // fallback: use first char of anomeric_prefix
    match prefix.chars().next() {
        Some('a') | Some('A') => "\u{03B1}",
        Some('b') | Some('B') => "\u{03B2}",
        Some('o') | Some('O') => "o",
        _ => "?",
    }
}

fn linkage_label_for(
    inner: &petgraph::graph::Graph<Monosaccharide, crabwurcs_core::Linkage>,
    child: NodeIndex,
    linkage: &crabwurcs_core::Linkage,
) -> String {
    let Some(residue) = inner.node_weight(child) else {
        return "?".to_string();
    };
    let anomer = anomer_char(residue.anomeric_symbol, &residue.anomeric_prefix);
    let positions = linkage
        .parent_positions()
        .map(|position| {
            if position.0 == 0 {
                "?".to_string()
            } else {
                position.0.to_string()
            }
        })
        .collect::<Vec<_>>()
        .join("/");
    let bridge = linkage
        .map_code
        .as_deref()
        .and_then(map_bridge_label)
        .map(|label| format!(" · {label}"))
        .unwrap_or_default();
    format!("{anomer} {positions}{bridge}")
}

fn map_bridge_label(map_code: &str) -> Option<&'static str> {
    match map_code {
        "*O*" => Some("Anhydro"),
        "*OC^XO*/3CO/6=O/3C" | "*1OC^X*2/3CO/5=O/3C" => Some("Py"),
        "*OC^SO*/3CO/6=O/3C" | "*1OC^SO*2/3CO/6=O/3C" => Some("(S)Py"),
        "*OC^RO*/3CO/6=O/3C" | "*1OC^RO*2/3CO/6=O/3C" => Some("(R)Py"),
        "*OSO*/3=O/3=O" => Some("S"),
        "*NS*/3=O/3=O" => Some("NS"),
        "*OPO*/3O/3=O" | "*1OP^X*2/3O/3=O" => Some("P"),
        "*OPOPO*/5O/5=O/3O/3=O" => Some("PyrP"),
        "*1NCCOP^XO*2/6O/6=O" => Some("PEtn"),
        "*NCCOP^XOP^X*/8O/8=O/6O/6=O" => Some("PPEtn"),
        _ => None,
    }
}

fn map_modification_label(map_code: &str) -> Option<&'static str> {
    match map_code {
        "*OC" => Some("Me"),
        "*OCC/3=O" => Some("Ac"),
        "*OSO/3=O/3=O" => Some("S"),
        "*NSO/3=O/3=O" => Some("NS"),
        "*OPO/3O/3=O" => Some("P"),
        "*NCC/3=O" => Some("NAc"),
        _ => None,
    }
}

// ── SVG rendering ──────────────────────────────────────────────────────────

pub fn render_svg(graph: &ResidueGraph) -> SnfgResult<String> {
    render_svg_with_options(graph, &RenderOptions::default())
}

pub fn render_svg_with_options(graph: &ResidueGraph, opts: &RenderOptions) -> SnfgResult<String> {
    let inner = graph.inner();
    if inner.node_count() == 0 {
        return Ok(empty_svg());
    }
    if graph.is_composition() {
        return render_composition_svg(graph, opts);
    }

    let root = graph.root().unwrap_or_else(|| NodeIndex::from(0u32));
    let layout = compute_layout(graph, root);

    let max_depth = layout.values().map(|li| li.x).fold(0.0f64, f64::max);
    let min_y = layout.values().map(|li| li.y).fold(f64::MAX, f64::min);
    let max_y = layout.values().map(|li| li.y).fold(f64::MIN, f64::max);

    let s = opts.scale;
    let pad = if graph.undefined_modifications().is_empty() {
        70.0 * s
    } else {
        110.0 * s
    };
    let canvas_w = (max_depth * H_SPACING * s) + 2.0 * pad;
    let canvas_h = (max_y - min_y) * s + 2.0 * pad;

    // RL orientation: root at right (x = canvas_w - pad), children extend left
    let offset_y = -min_y * s + pad;
    let to_canvas = |li: &LayoutInfo| -> (f64, f64) {
        (canvas_w - pad - li.x * H_SPACING * s, li.y * s + offset_y)
    };

    // ── SVG header ────────────────────────────────────────────────────
    let mut svg = format!(
        r#"<svg xmlns="http://www.w3.org/2000/svg" role="img" aria-label="SNFG glycan diagram" viewBox="0 0 {w} {h}" width="{w}" height="{h}">
<style>
  .bond {{ stroke: #000; stroke-width: {bw}; fill: none; stroke-linecap: round; }}
  .uncertain {{ stroke: #555; stroke-width: {ubw}; fill: none; stroke-linecap: round; stroke-dasharray: 8 7; }}
  .link {{ font-family: {ff}; font-size: {ls}px; fill: #000; text-anchor: middle; dominant-baseline: central; }}
  .node {{ fill: none; }}
  .res-label {{ font-family: {ff}; font-size: 11px; fill: #000; text-anchor: middle; dominant-baseline: central; }}
</style>
"#,
        w = canvas_w,
        h = canvas_h,
        bw = BOND_W * s,
        ubw = 2.5 * s,
        ff = opts.font_family,
        ls = LABEL_SIZE * s,
    );

    // ── Edges and linkage labels ──────────────────────────────────────
    for edge in inner.edge_references() {
        let (Some(parent_layout), Some(child_layout)) = (
            layout.get(&edge.source().index()),
            layout.get(&edge.target().index()),
        ) else {
            continue;
        };
        let (px, py) = to_canvas(parent_layout);
        let (cx, cy) = to_canvas(child_layout);
        let class = if edge.weight().repeat.is_some() || edge.weight().cyclic {
            "uncertain"
        } else {
            "bond"
        };
        svg.push_str(&format!(
            r#"<line x1="{px}" y1="{py}" x2="{cx}" y2="{cy}" class="{class}"/>
"#,
        ));
        if opts.show_linkages {
            draw_linkage_text(
                &mut svg,
                px,
                py,
                cx,
                cy,
                &linkage_label_for(inner, edge.target(), edge.weight()),
                s,
            );
        }
    }

    for undefined in graph.undefined_linkages() {
        let Some(child_layout) = layout.get(&undefined.child.index()) else {
            continue;
        };
        let (cx, cy) = to_canvas(child_layout);
        for (candidate_index, parent) in undefined.parents.iter().enumerate() {
            let Some(parent_layout) = layout.get(&parent.residue.index()) else {
                continue;
            };
            let (px, py) = to_canvas(parent_layout);
            svg.push_str(&format!(
                r#"<line x1="{px}" y1="{py}" x2="{cx}" y2="{cy}" class="uncertain"/>
"#,
            ));
            if opts.show_linkages && candidate_index == 0 {
                draw_linkage_text(&mut svg, px, py, cx, cy, "?", s);
            }
        }
    }

    for modification in graph.undefined_modifications() {
        let candidates = modification
            .parents
            .iter()
            .filter_map(|parent| layout.get(&parent.residue.index()))
            .map(&to_canvas)
            .collect::<Vec<_>>();
        if candidates.is_empty() {
            continue;
        }
        let label_x = candidates.iter().map(|(x, _)| *x).fold(f64::MIN, f64::max) + 65.0 * s;
        let label_y = candidates.iter().map(|(_, y)| y).sum::<f64>() / candidates.len() as f64;
        for (parent_x, parent_y) in &candidates {
            svg.push_str(&format!(
                r#"<line x1="{parent_x}" y1="{parent_y}" x2="{label_x}" y2="{label_y}" class="uncertain"/>
"#,
            ));
        }
        let label = map_modification_label(&modification.map_code).unwrap_or("Sub");
        svg.push_str(&format!(
            r##"<rect x="{}" y="{}" width="{}" height="{}" rx="{}" fill="#fff"/>
<text x="{label_x}" y="{label_y}" class="link" data-undefined-modification="true">{{{label}?}}</text>
"##,
            label_x - 42.0 * s,
            label_y - 18.0 * s,
            84.0 * s,
            36.0 * s,
            6.0 * s,
        ));
    }

    // ── Nodes ──────────────────────────────────────────────────────────
    for node_idx in inner.node_indices() {
        if let (Some(li), Some(residue)) =
            (layout.get(&node_idx.index()), inner.node_weight(node_idx))
        {
            let (cx, cy) = to_canvas(li);
            let symbol = symbol_for(residue)?;
            draw_shape(&mut svg, &symbol, cx, cy, NODE_R * s, opts);

            // sulfation / modification label above the shape
            let mod_label = build_modification_label(residue);
            if !mod_label.is_empty() {
                svg.push_str(&format!(
                    "<text x=\"{x}\" y=\"{y}\" font-family=\"{ff}\" font-size=\"12px\" fill=\"#000\" text-anchor=\"middle\" dominant-baseline=\"central\">{lbl}</text>\n",
                    x = cx, y = cy - NODE_R * s - 14.0 * s,
                    ff = opts.font_family,
                    lbl = mod_label,
                ));
            }

            if opts.show_labels {
                svg.push_str(&format!(
                    r#"<text x="{x}" y="{y}" class="res-label">{lbl}</text>
"#,
                    x = cx,
                    y = cy,
                    lbl = symbol.label,
                ));
            }
        }
    }

    svg.push_str("</svg>\n");
    Ok(svg)
}

fn render_composition_svg(graph: &ResidueGraph, opts: &RenderOptions) -> SnfgResult<String> {
    let mut groups: Vec<(String, Symbol, String, usize)> = Vec::new();
    for residue in graph.inner().node_weights() {
        let key = format!("{residue:?}");
        if let Some((_, _, _, count)) = groups.iter_mut().find(|(value, _, _, _)| *value == key) {
            *count += 1;
        } else {
            groups.push((
                key,
                symbol_for(residue)?,
                build_modification_label(residue),
                1,
            ));
        }
    }

    let scale = opts.scale;
    let spacing = 155.0 * scale;
    let width = (groups.len().max(1) as f64 * spacing) + 50.0 * scale;
    let height = 155.0 * scale;
    let mut svg = format!(
        r#"<svg xmlns="http://www.w3.org/2000/svg" role="img" aria-label="SNFG glycan composition" viewBox="0 0 {width} {height}" width="{width}" height="{height}">
<style>
  .res-label {{ font-family: {font}; font-size: 11px; fill: #000; text-anchor: middle; dominant-baseline: central; }}
  .count {{ font-family: {font}; font-size: {count_size}px; font-weight: 600; fill: #000; text-anchor: middle; }}
</style>
"#,
        font = opts.font_family,
        count_size = 18.0 * scale,
    );
    for (index, (_, symbol, modification, count)) in groups.iter().enumerate() {
        let x = 75.0 * scale + index as f64 * spacing;
        let y = 62.0 * scale;
        draw_shape(&mut svg, symbol, x, y, NODE_R * scale, opts);
        if opts.show_labels {
            svg.push_str(&format!(
                r#"<text x="{x}" y="{y}" class="res-label">{}</text>
"#,
                symbol.label
            ));
        }
        if !modification.is_empty() {
            svg.push_str(&format!(
                r#"<text x="{x}" y="{}" class="res-label">{modification}</text>
"#,
                y - 42.0 * scale
            ));
        }
        svg.push_str(&format!(
            r#"<text x="{x}" y="{}" class="count">×{count}</text>
"#,
            y + 55.0 * scale
        ));
    }
    svg.push_str("</svg>\n");
    Ok(svg)
}

fn draw_linkage_text(
    svg: &mut String,
    px: f64,
    py: f64,
    cx: f64,
    cy: f64,
    label: &str,
    scale: f64,
) {
    let mx = (px + cx) / 2.0;
    let my = (py + cy) / 2.0;
    let dx = cx - px;
    let dy = cy - py;
    let len = (dx * dx + dy * dy).sqrt();
    let (ox, oy) = if len > 1.0 {
        let first = (-dy / len * 14.0 * scale, dx / len * 14.0 * scale);
        let second = (dy / len * 14.0 * scale, -dx / len * 14.0 * scale);
        if first.1 < 0.0 {
            first
        } else {
            second
        }
    } else {
        (0.0, -14.0 * scale)
    };
    let mut angle = dy.atan2(dx).to_degrees();
    if angle > 90.0 {
        angle -= 180.0;
    } else if angle < -90.0 {
        angle += 180.0;
    }
    svg.push_str(&format!(
        r#"<text x="0" y="0" class="link" transform="translate({x},{y}) rotate({angle})">{label}</text>
"#,
        x = mx + ox,
        y = my + oy,
    ));
}

fn empty_svg() -> String {
    concat!(
        r#"<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 120 40" width="120" height="40">
  <text x="10" y="25" font-family="sans-serif" font-size="11" fill=""#,
        "#999",
        r#"">(empty)</text>
</svg>
"#,
    )
    .to_string()
}

// ── Modification label ─────────────────────────────────────────────────────

/// Build a short sulfation / modification label.
/// e.g. GlcNS6S → "S6S", IdoA2S → "2S", GlcNS → "S", GlcA3S → "3S"
/// Convention: the NSquare shape already implies N-modification,
/// so N-sulfation alone is shown as "S", not "N".
fn build_modification_label(res: &Monosaccharide) -> String {
    let mut o_sulfo_positions: Vec<u8> = Vec::new();
    let mut has_n_sulfo = false;

    for m in &res.modifications {
        let desc = &m.descriptor;
        if desc.contains("NSO") {
            has_n_sulfo = true;
        }
        if desc.contains("OSO") {
            o_sulfo_positions.push(m.position.0);
        }
    }

    if o_sulfo_positions.is_empty() && !has_n_sulfo {
        return String::new();
    }

    o_sulfo_positions.sort();
    o_sulfo_positions.dedup();

    let mut parts = Vec::new();
    if has_n_sulfo {
        if o_sulfo_positions.is_empty() {
            // GlcNS with no O-sulfation → just "S" (NSquare implies N)
            parts.push("S".to_string());
        } else {
            // GlcNS6S → "S6S"
            parts.push("S".to_string());
        }
    }
    for pos in &o_sulfo_positions {
        parts.push(format!("{}S", pos));
    }
    parts.join("")
}

fn draw_shape(svg: &mut String, sym: &Symbol, cx: f64, cy: f64, r: f64, opts: &RenderOptions) {
    let fill = if opts.colour { sym.fill } else { "none" };
    let stroke = colour::STROKE;
    let sw = 2.0 * opts.scale;

    match sym.shape {
        Shape::Circle => {
            svg.push_str(&format!(
                r#"<circle cx="{x}" cy="{y}" r="{r}" fill="{fill}" stroke="{stroke}" stroke-width="{sw}"/>
"#,
                x = cx, y = cy, r = r, fill = fill, stroke = stroke, sw = sw,
            ));
        }
        Shape::Square => {
            let h = r;
            svg.push_str(&format!(
                r#"<rect x="{x}" y="{y}" width="{w}" height="{h}" fill="{fill}" stroke="{stroke}" stroke-width="{sw}"/>
"#,
                x = cx - h, y = cy - h, w = h * 2.0, h = h * 2.0,
                fill = fill, stroke = stroke, sw = sw,
            ));
        }
        Shape::NSquare => {
            // white square with coloured top-left triangle (N-modified, non-acetylated)
            let h = r;
            // white background square
            svg.push_str(&format!(
                r#"<rect x="{x}" y="{y}" width="{w}" height="{h}" fill="white" stroke="{stroke}" stroke-width="{sw}"/>
"#,
                x = cx - h, y = cy - h, w = h * 2.0, h = h * 2.0,
                stroke = stroke, sw = sw,
            ));
            // coloured top-left triangle
            svg.push_str(&format!(
                r#"<polygon points="{x1},{y1} {x2},{y2} {x3},{y3} {x4},{y4}" fill="{fill}" stroke="none"/>
"#,
                x1 = cx - h, y1 = cy - h,
                x2 = cx + h, y2 = cy - h,
                x3 = cx + h, y3 = cy + h,
                x4 = cx - h, y4 = cy - h,
                fill = fill,
            ));
            // inner dividing lines
            svg.push_str(&format!(
                r#"<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="{stroke}" stroke-width="1.5"/>
"#,
                x1 = cx - h, y1 = cy - h,
                x2 = cx + h, y2 = cy + h,
                stroke = stroke,
            ));
        }
        Shape::Triangle => {
            let h = r * 1.732; // equilateral
            let w = r;
            svg.push_str(&format!(
                r#"<polygon points="{x1},{y1} {x2},{y2} {x3},{y3}" fill="{fill}" stroke="{stroke}" stroke-width="{sw}"/>
"#,
                x1 = cx, y1 = cy - h * 0.667,
                x2 = cx - w, y2 = cy + h * 0.333,
                x3 = cx + w, y3 = cy + h * 0.333,
                fill = fill, stroke = stroke, sw = sw,
            ));
        }
        Shape::Diamond => {
            let d = r * 1.2;
            svg.push_str(&format!(
                r#"<polygon points="{x1},{y1} {x2},{y2} {x3},{y3} {x4},{y4}" fill="{fill}" stroke="{stroke}" stroke-width="{sw}"/>
"#,
                x1 = cx, y1 = cy - d,
                x2 = cx + d, y2 = cy,
                x3 = cx, y3 = cy + d,
                x4 = cx - d, y4 = cy,
                fill = fill, stroke = stroke, sw = sw,
            ));
        }
        Shape::SplitDiamondTop => {
            // GlcA, GalA, ManA — top half coloured
            let d = r * 1.2;
            svg.push_str(&format!(
                r#"<polygon points="{x1},{y1} {x2},{y2} {x3},{y3} {x4},{y4}" fill="white" stroke="{stroke}" stroke-width="{sw}"/>
"#,
                x1 = cx, y1 = cy - d, x2 = cx + d, y2 = cy,
                x3 = cx, y3 = cy + d, x4 = cx - d, y4 = cy,
                stroke = stroke, sw = sw,
            ));
            // coloured top triangle
            svg.push_str(&format!(
                r#"<polygon points="{x1},{y1} {x2},{y2} {x3},{y3} {x4},{y4}" fill="{fill}" stroke="none"/>
"#,
                x1 = cx - d, y1 = cy, x2 = cx, y2 = cy - d,
                x3 = cx + d, y3 = cy, x4 = cx - d, y4 = cy,
                fill = fill,
            ));
            // horizontal dividing line
            svg.push_str(&format!(
                r#"<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="{stroke}" stroke-width="1.5"/>
"#,
                x1 = cx - d, y1 = cy, x2 = cx + d, y2 = cy, stroke = stroke,
            ));
        }
        Shape::SplitDiamondBottom => {
            // IdoA — bottom half coloured (brown)
            let d = r * 1.2;
            svg.push_str(&format!(
                r#"<polygon points="{x1},{y1} {x2},{y2} {x3},{y3} {x4},{y4}" fill="white" stroke="{stroke}" stroke-width="{sw}"/>
"#,
                x1 = cx, y1 = cy - d, x2 = cx + d, y2 = cy,
                x3 = cx, y3 = cy + d, x4 = cx - d, y4 = cy,
                stroke = stroke, sw = sw,
            ));
            // coloured bottom triangle
            svg.push_str(&format!(
                r#"<polygon points="{x1},{y1} {x2},{y2} {x3},{y3} {x4},{y4}" fill="{fill}" stroke="none"/>
"#,
                x1 = cx - d, y1 = cy, x2 = cx, y2 = cy + d,
                x3 = cx + d, y3 = cy, x4 = cx - d, y4 = cy,
                fill = fill,
            ));
            // horizontal dividing line
            svg.push_str(&format!(
                r#"<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="{stroke}" stroke-width="1.5"/>
"#,
                x1 = cx - d, y1 = cy, x2 = cx + d, y2 = cy, stroke = stroke,
            ));
        }
        Shape::Star => {
            let outer = r * 1.2;
            let inner = r * 0.5;
            svg.push_str(&draw_regular_points(cx, cy, outer, inner, 5, 2));
            svg.push_str(&format!(
                " fill=\"{}\" stroke=\"{}\" stroke-width=\"{}\"/>",
                fill, stroke, sw
            ));
            svg.push('\n');
        }
        Shape::Hexagon => {
            svg.push_str(&draw_regular_points(cx, cy, r * 1.05, r * 1.05, 6, 1));
            svg.push_str(&format!(
                " fill=\"{}\" stroke=\"{}\" stroke-width=\"{}\"/>",
                fill, stroke, sw
            ));
            svg.push('\n');
        }
        Shape::Pentagon => {
            svg.push_str(&draw_regular_points(cx, cy, r * 1.05, r * 1.05, 5, 1));
            svg.push_str(&format!(
                " fill=\"{}\" stroke=\"{}\" stroke-width=\"{}\"/>",
                fill, stroke, sw
            ));
            svg.push('\n');
        }
    }
}

fn draw_regular_points(
    cx: f64,
    cy: f64,
    outer: f64,
    inner: f64,
    sides: usize,
    cycles: usize,
) -> String {
    let total = sides * cycles;
    let mut pts = String::from("<polygon points=\"");
    for i in 0..total {
        let angle =
            2.0 * std::f64::consts::PI * i as f64 / total as f64 - std::f64::consts::FRAC_PI_2;
        let r = if i % 2 == 0 { outer } else { inner };
        use std::fmt::Write;
        write!(&mut pts, "{}", cx + r * angle.cos()).unwrap();
        write!(&mut pts, ",{} ", cy + r * angle.sin()).unwrap();
    }
    pts.pop(); // trailing space
    pts.push('"');
    pts
}

// ── Tests ──────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn every_glycoshape_molecular_record_renders_with_known_symbols() {
        let mut records = 0usize;
        let mut unknown = Vec::new();
        let corpus_lines = include_str!("../../crabwurcs/data/glycoshape_notations.tsv")
            .lines()
            .chain(include_str!("../../crabwurcs/data/glycoshape_derived_notations.tsv").lines());
        for line in corpus_lines {
            let wurcs = line.split('\t').next().unwrap();
            let graph = crabwurcs_core::parse_wurcs(wurcs).unwrap();
            for residue in graph.inner().node_weights() {
                let symbol = symbol_for(residue).unwrap();
                if symbol.label == "?" {
                    unknown.push(format!("{residue:?}"));
                }
            }
            let svg = render_svg(&graph).unwrap();
            assert!(svg.contains("role=\"img\""), "{wurcs}");
            assert!(svg.contains("viewBox="), "{wurcs}");
            records += 1;
        }
        assert_eq!(records, 938);
        unknown.sort();
        unknown.dedup();
        assert!(unknown.is_empty(), "unknown SNFG symbols: {unknown:?}");
    }
    use crabwurcs_core::parse_wurcs;

    fn parse(wurcs: &str) -> ResidueGraph {
        parse_wurcs(wurcs).expect("parse WURCS")
    }

    #[test]
    fn test_symbol_glc() {
        let g = parse("WURCS=2.0/2,2,1/[u2122h][a2122h-1b_1-5]/1-2/a4-b1");
        let res = g.residue(g.root().unwrap()).unwrap();
        let sym = symbol_for(res).unwrap();
        assert_eq!(sym.shape, Shape::Circle);
        assert_eq!(sym.fill, colour::GLC);
    }

    #[test]
    fn test_symbol_glcnac() {
        let g = parse("WURCS=2.0/2,2,1/[u2122h_2*NCC/3=O][a2122h-1b_1-5]/1-1-2/a4-b1");
        let res = g.residue(g.root().unwrap()).unwrap();
        let sym = symbol_for(res).unwrap();
        assert_eq!(sym.shape, Shape::Square);
        assert_eq!(sym.fill, colour::GLC);
    }

    #[test]
    fn test_symbol_gal() {
        let g = parse("WURCS=2.0/2,2,1/[u2112h][a2112h-1b_1-5]/1-2/a3-b1");
        let res = g.residue(g.root().unwrap()).unwrap();
        let sym = symbol_for(res).unwrap();
        assert_eq!(sym.shape, Shape::Circle);
        assert_eq!(sym.fill, colour::GAL);
    }

    #[test]
    fn test_symbol_man() {
        let g = parse("WURCS=2.0/2,2,1/[u1122h][a1122h-1a_1-5]/1-2/a4-b1");
        let res = g.residue(g.root().unwrap()).unwrap();
        let sym = symbol_for(res).unwrap();
        assert_eq!(sym.shape, Shape::Circle);
        assert_eq!(sym.fill, colour::MAN);
    }

    #[test]
    fn rare_hexose_epimers_use_their_snfg_families() {
        for (code, label, fill) in [
            ("2111", "Alt", colour::ALT),
            ("1121", "Gul", colour::GUL),
            ("2222", "All", colour::ALL),
            ("2221", "Tal", colour::TAL),
            ("2121", "Ido", colour::IDO),
        ] {
            let g = parse(&format!("WURCS=2.0/1,1,0/[u{code}h]/1/"));
            let symbol = symbol_for(g.residue(g.root().unwrap()).unwrap()).unwrap();
            assert_eq!(symbol.shape, Shape::Circle, "{code}");
            assert_eq!(symbol.label, label, "{code}");
            assert_eq!(symbol.fill, fill, "{code}");
        }
    }

    #[test]
    fn test_symbol_fuc() {
        let g = parse("WURCS=2.0/2,2,1/[u1221m][a1221m-1a_1-5]/1-2/a3-b1");
        let res = g.residue(g.root().unwrap()).unwrap();
        let sym = symbol_for(res).unwrap();
        assert_eq!(sym.shape, Shape::Triangle);
        assert_eq!(sym.fill, colour::FUC);
    }

    #[test]
    fn test_symbol_neu5ac() {
        let g = parse("WURCS=2.0/2,2,1/[u2112h][Aad21122h-2a_2-6_5*NCC/3=O]/1-2/a3-b2");
        let children: Vec<_> = g
            .inner()
            .neighbors_directed(g.root().unwrap(), Direction::Outgoing)
            .collect();
        assert!(!children.is_empty());
        let neu = g.residue(children[0]).unwrap();
        let sym = symbol_for(neu).unwrap();
        assert_eq!(sym.shape, Shape::Diamond);
        assert_eq!(sym.fill, colour::NEU5AC);
    }

    #[test]
    fn glycoshape_sialic_acid_palette_and_fruf_symbol_are_exact() {
        let neu5gc = parse("WURCS=2.0/1,1,0/[AUd21122h_5*NCCO/3=O]/1/");
        let neu5gc = symbol_for(neu5gc.residue(neu5gc.root().unwrap()).unwrap()).unwrap();
        assert_eq!(neu5gc.shape, Shape::Diamond);
        assert_eq!(neu5gc.fill, "#91D3E3");

        let fructan = parse("WURCS=2.0/2,3,2/[hU122h][ha122h-2b_2-5]/1-2-2/a1-b2_b1-c2");
        for residue in fructan.inner().node_weights() {
            let symbol = symbol_for(residue).unwrap();
            assert_eq!(symbol.shape, Shape::Pentagon);
            assert_eq!(symbol.fill, colour::MAN);
            assert_eq!(symbol.label, "Fru");
        }
    }

    #[test]
    fn test_render_linear() {
        let g = parse("WURCS=2.0/2,2,1/[u2112h][a2112h-1b_1-5]/1-2/a3-b1");
        let svg = render_svg(&g).unwrap();
        assert!(svg.contains("<svg"));
        assert!(svg.contains("</svg>"));
        assert!(svg.contains("class=\"bond\""));
    }

    #[test]
    fn test_kdo_and_bac_use_reference_flat_hexagons() {
        let kdo = parse("WURCS=2.0/1,1,0/[AUd1122h]/1/");
        let kdo_symbol = symbol_for(kdo.residue(kdo.root().unwrap()).unwrap()).unwrap();
        assert_eq!(kdo_symbol.shape, Shape::Hexagon);
        assert_eq!(kdo_symbol.fill, colour::KDO);
        assert_eq!(kdo_symbol.label, "KDO");

        let bac = parse("WURCS=2.0/1,1,0/[u2122m_2*NCC/3=O_4*NCC/3=O]/1/");
        let bac_symbol = symbol_for(bac.residue(bac.root().unwrap()).unwrap()).unwrap();
        assert_eq!(bac_symbol.shape, Shape::Hexagon);
        assert_eq!(bac_symbol.fill, colour::GLC);
        assert_eq!(bac_symbol.label, "Bac");
    }

    #[test]
    fn test_render_empty() {
        let g = ResidueGraph::new();
        let svg = render_svg(&g).unwrap();
        assert!(svg.contains("empty"));
    }

    #[test]
    fn test_render_branched() {
        let g = parse(
            "WURCS=2.0/3,3,2/[u2112h_2*NCC/3=O][a2112h-1a_1-5_2*NCC/3=O][Aad21122h-2a_2-6_5*NCC/3=O]/1-2-3/a3-b1_a6-c2",
        );
        let svg = render_svg(&g).unwrap();
        assert!(svg.contains("<svg"));
        let bond_count = svg.matches("class=\"bond\"").count();
        assert_eq!(bond_count, 2);
        assert!(!svg.contains("<rect class=\"bg\""));
        assert!(svg.contains("rotate("));
    }

    #[test]
    fn test_render_map_bridge_labels_its_chemistry() {
        let g = parse("WURCS=2.0/2,2,1/[hxh][a2122h-1b_1-5]/1-2/a3n2-b1n1*1NCCOP^XO*2/6O/6=O");
        let svg = render_svg(&g).unwrap();
        assert!(svg.contains("PEtn"), "{svg}");
    }

    #[test]
    fn test_render_undefined_modification_with_candidate_bonds() {
        let g = parse("WURCS=2.0/2,2,1/[u2122h][u2112h]/1-2/a?|b?}*OCC/3=O");
        let svg = render_svg(&g).unwrap();
        assert!(svg.contains("data-undefined-modification=\"true\""));
        assert!(svg.contains("{Ac?}"));
        assert_eq!(svg.matches("class=\"uncertain\"").count(), 2);
    }

    #[test]
    fn test_render_complex_n_glycan() {
        let g = parse(
            "WURCS=2.0/6,8,7/[u2122h_2*NCC/3=O][a1221m-1a_1-5][a2122h-1b_1-5_2*NCC/3=O][a1122h-1b_1-5][a1122h-1a_1-5][a2112h-1b_1-5]/1-2-3-4-5-5-2-6/a3-b1_a4-c1_a6-g1_c4-d1_d3-e1_d6-f1_g4-h1",
        );
        let svg = render_svg(&g).unwrap();
        assert!(svg.contains("<svg"));
        assert!(svg.contains("\u{03B1}") || svg.contains("\u{03B2}"));
    }

    #[test]
    fn gs00955_uses_standard_n_glycan_branch_order() {
        let g = parse(
            "WURCS=2.0/6,12,11/[u2122h_2*NCC/3=O][a2122h-1b_1-5_2*NCC/3=O][a1122h-1b_1-5][a1122h-1a_1-5][a2112h-1b_1-5][a1221m-1a_1-5]/1-2-3-4-2-5-2-5-4-2-6-5/a4-b1_b4-c1_c3-d1_c6-i1_d2-e1_d4-g1_e4-f1_g4-h1_i2-j1_j3-k1_j4-l1",
        );
        let layout = compute_layout(&g, g.root().unwrap());
        let central_man = NodeIndex::new(2);
        let arm_y = |position| {
            let child = g
                .inner()
                .edges_directed(central_man, Direction::Outgoing)
                .find(|edge| edge.weight().parent_position.0 == position)
                .unwrap()
                .target();
            layout[&child.index()].y
        };
        assert!(arm_y(6) < arm_y(3), "α1-6 must be above α1-3");

        let alpha3_man = NodeIndex::new(3);
        let branch_y = |position| {
            let child = g
                .inner()
                .edges_directed(alpha3_man, Direction::Outgoing)
                .find(|edge| edge.weight().parent_position.0 == position)
                .unwrap()
                .target();
            layout[&child.index()].y
        };
        assert!(branch_y(4) < branch_y(2), "β1-4 must be above β1-2");

        let fucose: NodeIndex = NodeIndex::new(10);
        let fucose_parent: NodeIndex = NodeIndex::new(9);
        assert_eq!(
            layout[&fucose.index()].x,
            layout[&fucose_parent.index()].x,
            "core fucose must be vertical at its parent's depth"
        );
        assert_ne!(
            layout[&fucose.index()].y,
            layout[&fucose_parent.index()].y,
            "core fucose must not overprint its parent"
        );
    }

    #[test]
    fn terminal_fucose_reserves_a_lane_instead_of_covering_another_residue() {
        let g = parse(
            "WURCS=2.0/4,4,3/[u2112h_2*NCC/3=O][a2112h-1b_1-5][a2122h-1b_1-5_2*NCC/3=O][a1221m-1a_1-5]/1-2-3-4/a3-b1_a6-c1_c3-d1",
        );
        let layout = compute_layout(&g, g.root().unwrap());
        let coordinates = layout
            .values()
            .map(|position| (position.x as i32, position.y.round() as i32))
            .collect::<std::collections::HashSet<_>>();
        assert_eq!(coordinates.len(), g.node_count());
    }

    #[test]
    fn test_render_undefined_fragment_includes_all_components() {
        let g = parse(
            "WURCS=2.0/6,11,10/[a2122h-1x_1-5_2*NCC/3=O][a2122h-1b_1-5_2*NCC/3=O][a1122h-1b_1-5][a1122h-1a_1-5][a2112h-1b_1-5][a1221m-1a_1-5]/1-2-3-4-2-5-4-2-6-2-5/a4-b1_a6-i1_b4-c1_c3-d1_c6-g1_d2-e1_e4-f1_g2-h1_j4-k1_j1-d4|d6|g4|g6}",
        );
        assert_eq!(g.undefined_linkages().len(), 1);
        let svg = render_svg_with_options(
            &g,
            &RenderOptions {
                show_labels: true,
                ..RenderOptions::default()
            },
        )
        .unwrap();
        assert_eq!(svg.matches("class=\"res-label\"").count(), 11);
        assert_eq!(svg.matches("class=\"uncertain\"").count(), 2);
    }

    #[test]
    fn test_render_composition_does_not_drop_disconnected_residues() {
        let g = parse(
            "WURCS=2.0/4,15,0+/[AUd21122h_5*NCC/3=O][uxxxxh_2*NCC/3=O][uxxxxh][u1221m]/1-2-2-2-2-2-2-3-3-3-4-4-4-4-4/",
        );
        let svg = render_svg_with_options(
            &g,
            &RenderOptions {
                show_labels: true,
                ..RenderOptions::default()
            },
        )
        .unwrap();
        assert!(svg.contains("aria-label=\"SNFG glycan composition\""));
        assert_eq!(svg.matches("class=\"count\"").count(), 4);
        assert!(svg.contains("×6"));
        assert!(svg.contains("×5"));
    }
}
