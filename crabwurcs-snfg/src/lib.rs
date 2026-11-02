use crabwurcs_core::{
    Monosaccharide, MotifError, MotifMatch, RepeatCount, ResidueGraph, ResidueKind,
    classify_residue, classify_residue_configuration, find_motif_matches,
};
use petgraph::Direction;
use petgraph::graph::NodeIndex;
use petgraph::visit::EdgeRef;
use std::collections::{BTreeSet, HashMap};
use std::str::FromStr;
use thiserror::Error;

// ── Error types ────────────────────────────────────────────────────────────

#[derive(Debug, Error)]
pub enum SnfgError {
    #[error(
        "no SNFG symbol for monosaccharide: backbone_len={}, code={}",
        .0.backbone_length,
        .0.skeleton_code
    )]
    UnknownSymbol(Box<Monosaccharide>),

    #[error(transparent)]
    Core(#[from] crabwurcs_core::CoreError),

    #[error(transparent)]
    Motif(#[from] MotifError),

    #[error("failed to parse generated SVG for PNG rendering: {0}")]
    SvgParse(String),

    #[error("PNG dimensions overflow or are unsupported")]
    RasterDimensions,

    #[error("could not allocate the PNG raster surface")]
    RasterAllocation,

    #[error("failed to encode PNG: {0}")]
    PngEncoding(String),

    #[error("no residue {index} in a glycan of {count} residues")]
    InvalidResidueIndex { index: usize, count: usize },

    #[error("no bond {index} in a glycan of {count} bonds")]
    InvalidBondIndex { index: usize, count: usize },
}

pub type SnfgResult<T> = Result<T, SnfgError>;

// ── SNFG shapes ────────────────────────────────────────────────────────────

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Shape {
    Circle,
    Square,
    NSquare,
    Triangle,
    DividedTriangle,
    Diamond,
    FlatDiamond,
    SplitDiamondTop,    // uronic acid – top half coloured
    SplitDiamondBottom, // uronic acid – bottom half coloured (IdoA)
    FlatRectangle,
    Star,
    Hexagon,
    FlatHexagon,
    Pentagon,
}

// ── SNFG colours (glycoshape.io / standard SNFG palette) ───────────────────

pub mod colour {
    pub const WHITE: &str = "#FFFFFF";
    pub const BLUE: &str = "#0072BC";
    pub const GREEN: &str = "#00A651";
    pub const YELLOW: &str = "#FFD400";
    pub const ORANGE: &str = "#F47920";
    pub const PINK: &str = "#F69EA1";
    pub const PURPLE: &str = "#A54399";
    pub const LIGHT_BLUE: &str = "#8FCCE9";
    pub const BROWN: &str = "#A17A4D";
    pub const RED: &str = "#ED1C24";

    pub const GLC: &str = BLUE;
    pub const GAL: &str = YELLOW;
    pub const MAN: &str = GREEN;
    pub const FUC: &str = RED;
    pub const RHA: &str = GREEN;
    pub const NEU5AC: &str = PURPLE;
    pub const NEU5GC: &str = LIGHT_BLUE;
    pub const KDN: &str = GREEN;
    pub const IDOA: &str = BROWN;
    pub const KDO: &str = YELLOW;
    pub const XYL: &str = ORANGE;
    pub const GUL: &str = ORANGE;
    pub const ALT: &str = PINK;
    pub const ALL: &str = PURPLE;
    pub const TAL: &str = LIGHT_BLUE;
    pub const IDO: &str = BROWN;
    pub const UNKNOWN: &str = WHITE;

    pub const STROKE: &str = "#000000";
    pub const BOND: &str = "#000000";
    pub const LINKAGE_TEXT: &str = "#000000";
}

#[derive(Debug, Clone)]
pub struct Symbol {
    pub shape: Shape,
    pub fill: &'static str,
    pub label: String,
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
            label: acid.into(),
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
            label: if has_nac { nac } else { amine }.into(),
        };
    }
    Symbol {
        shape: Shape::Circle,
        fill,
        label: neutral.into(),
    }
}

fn registered_symbol(kind: ResidueKind, display_name: Option<&str>) -> Symbol {
    use ResidueKind::*;
    let (shape, fill) = match kind {
        Hex => (Shape::Circle, colour::WHITE),
        Glc => (Shape::Circle, colour::BLUE),
        Man => (Shape::Circle, colour::GREEN),
        Gal => (Shape::Circle, colour::YELLOW),
        Gul => (Shape::Circle, colour::ORANGE),
        Alt => (Shape::Circle, colour::PINK),
        All => (Shape::Circle, colour::PURPLE),
        Tal => (Shape::Circle, colour::LIGHT_BLUE),
        Ido => (Shape::Circle, colour::BROWN),

        HexNAc => (Shape::Square, colour::WHITE),
        GlcNAc => (Shape::Square, colour::BLUE),
        ManNAc => (Shape::Square, colour::GREEN),
        GalNAc => (Shape::Square, colour::YELLOW),
        GulNAc => (Shape::Square, colour::ORANGE),
        AltNAc => (Shape::Square, colour::PINK),
        AllNAc => (Shape::Square, colour::PURPLE),
        TalNAc => (Shape::Square, colour::LIGHT_BLUE),
        IdoNAc => (Shape::Square, colour::BROWN),

        HexN => (Shape::NSquare, colour::WHITE),
        GlcN => (Shape::NSquare, colour::BLUE),
        ManN => (Shape::NSquare, colour::GREEN),
        GalN => (Shape::NSquare, colour::YELLOW),
        GulN => (Shape::NSquare, colour::ORANGE),
        AltN => (Shape::NSquare, colour::PINK),
        AllN => (Shape::NSquare, colour::PURPLE),
        TalN => (Shape::NSquare, colour::LIGHT_BLUE),
        IdoN => (Shape::NSquare, colour::BROWN),

        HexA => (Shape::SplitDiamondTop, colour::WHITE),
        GlcA => (Shape::SplitDiamondTop, colour::BLUE),
        ManA => (Shape::SplitDiamondTop, colour::GREEN),
        GalA => (Shape::SplitDiamondTop, colour::YELLOW),
        GulA => (Shape::SplitDiamondTop, colour::ORANGE),
        AltA => (Shape::SplitDiamondTop, colour::PINK),
        AllA => (Shape::SplitDiamondTop, colour::PURPLE),
        TalA => (Shape::SplitDiamondTop, colour::LIGHT_BLUE),
        IdoA => (Shape::SplitDiamondBottom, colour::BROWN),

        DHex => (Shape::Triangle, colour::WHITE),
        Qui => (Shape::Triangle, colour::BLUE),
        Rha => (Shape::Triangle, colour::GREEN),
        SixDGul => (Shape::Triangle, colour::ORANGE),
        SixDAlt => (Shape::Triangle, colour::PINK),
        SixDTal => (Shape::Triangle, colour::LIGHT_BLUE),
        Fuc => (Shape::Triangle, colour::RED),

        DHexNAc => (Shape::DividedTriangle, colour::WHITE),
        QuiNAc => (Shape::DividedTriangle, colour::BLUE),
        RhaNAc => (Shape::DividedTriangle, colour::GREEN),
        SixDAltNAc => (Shape::DividedTriangle, colour::PINK),
        SixDTalNAc => (Shape::DividedTriangle, colour::LIGHT_BLUE),
        FucNAc => (Shape::DividedTriangle, colour::RED),

        DDHex => (Shape::FlatRectangle, colour::WHITE),
        Oli => (Shape::FlatRectangle, colour::BLUE),
        Tyv => (Shape::FlatRectangle, colour::GREEN),
        Abe => (Shape::FlatRectangle, colour::ORANGE),
        Par => (Shape::FlatRectangle, colour::PINK),
        Dig => (Shape::FlatRectangle, colour::PURPLE),
        Col => (Shape::FlatRectangle, colour::LIGHT_BLUE),

        Pen => (Shape::Star, colour::WHITE),
        Ara => (Shape::Star, colour::GREEN),
        Lyx => (Shape::Star, colour::YELLOW),
        Xyl => (Shape::Star, colour::ORANGE),
        Rib => (Shape::Star, colour::PINK),

        NulO => (Shape::Diamond, colour::WHITE),
        Kdn => (Shape::Diamond, colour::GREEN),
        Neu5Ac => (Shape::Diamond, colour::PURPLE),
        Neu5Gc => (Shape::Diamond, colour::LIGHT_BLUE),
        Neu => (Shape::Diamond, colour::BROWN),
        Sia => (Shape::Diamond, colour::RED),

        DDNulO => (Shape::FlatDiamond, colour::WHITE),
        Pse => (Shape::FlatDiamond, colour::GREEN),
        Leg => (Shape::FlatDiamond, colour::YELLOW),
        Aci => (Shape::FlatDiamond, colour::PINK),
        FourELeg => (Shape::FlatDiamond, colour::LIGHT_BLUE),

        Unknown => (Shape::FlatHexagon, colour::WHITE),
        Bac => (Shape::FlatHexagon, colour::BLUE),
        LDManHep => (Shape::FlatHexagon, colour::GREEN),
        Kdo => (Shape::FlatHexagon, colour::YELLOW),
        Dha => (Shape::FlatHexagon, colour::ORANGE),
        DDManHep => (Shape::FlatHexagon, colour::PINK),
        MurNAc => (Shape::FlatHexagon, colour::PURPLE),
        MurNGc => (Shape::FlatHexagon, colour::LIGHT_BLUE),
        Mur => (Shape::FlatHexagon, colour::BROWN),

        Assigned => (Shape::Pentagon, colour::WHITE),
        Api => (Shape::Pentagon, colour::BLUE),
        Fru => (Shape::Pentagon, colour::GREEN),
        Tag => (Shape::Pentagon, colour::YELLOW),
        Sor => (Shape::Pentagon, colour::ORANGE),
        Psi => (Shape::Pentagon, colour::PINK),
    };
    let label = if kind == Assigned {
        display_name
            .and_then(|name| name.chars().find(char::is_ascii_alphabetic))
            .map(|character| character.to_ascii_uppercase().to_string())
            .unwrap_or_else(|| "?".into())
    } else {
        kind.canonical_name().to_string()
    };
    Symbol { shape, fill, label }
}

// ── Monosaccharide → SNFG symbol ───────────────────────────────────────────

/// SNFG symbol of one residue.
///
/// The symbol identifies the monosaccharide whatever its absolute
/// configuration (SNFG 2.0: "less common configurations need to be stated in
/// a figure legend or by adding letters inside the symbol"), so D-Fuc is the
/// red triangle of Fuc and D-Rha the green triangle of Rha. The label of an
/// enantiomer carries its configuration (`D-Fuc`); [`render_svg`] writes the
/// configuration letter inside the symbol.
pub fn symbol_for(residue: &Monosaccharide) -> SnfgResult<Symbol> {
    if let Some((kind, inverted)) = classify_residue_configuration(residue) {
        let mut symbol = registered_symbol(kind, residue.display_name.as_deref());
        if inverted && let Some(letter) = configuration_mark(residue) {
            symbol.label = format!("{letter}-{}", symbol.label);
        }
        return Ok(symbol);
    }
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
            label: "Fru".into(),
        });
    }

    // L-Sorbose (ketopentose) also uses the 'h' prefix pattern
    if residue.anomeric_prefix.starts_with('h') && bare_str == "121" {
        return Ok(Symbol {
            shape: Shape::Star,
            fill: colour::MAN,
            label: "Sor".into(),
        });
    }

    // WURCS encodes ulosonic-acid oxidation in the leading `A*` carbon
    // descriptors rather than in the trailing skeleton string. KDO and Bac
    // belong to SNFG's flat-hexagon "unknown/other" family.
    if bare_str == "1122" && residue.anomeric_prefix.starts_with('A') {
        return Ok(Symbol {
            shape: Shape::Hexagon,
            fill: colour::KDO,
            label: "KDO".into(),
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
            label: "Bac".into(),
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
            ));
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
            ));
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
            ));
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
            ));
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
            ));
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
            ));
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
            ));
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
            ));
        }
        _ => {}
    }

    // ── Hexoses ────────────────────────────────────────────────────────
    match bare_str {
        "2122" if has_acid => {
            return Ok(Symbol {
                shape: Shape::SplitDiamondTop,
                fill: colour::GLC,
                label: "GlcA".into(),
            });
        }
        "2121" if has_acid => {
            return Ok(Symbol {
                shape: Shape::SplitDiamondBottom,
                fill: colour::IDOA,
                label: "IdoA".into(),
            });
        }
        "2122" | "2121" if has_n_mod => {
            if has_nac {
                return Ok(Symbol {
                    shape: Shape::Square,
                    fill: colour::GLC,
                    label: "GlcNAc".into(),
                });
            } else {
                return Ok(Symbol {
                    shape: Shape::NSquare,
                    fill: colour::GLC,
                    label: "GlcN".into(),
                });
            }
        }
        "2122" | "2121" => {
            return Ok(Symbol {
                shape: Shape::Circle,
                fill: colour::GLC,
                label: "Glc".into(),
            });
        }
        "2112" | "2111" if has_acid => {
            return Ok(Symbol {
                shape: Shape::SplitDiamondTop,
                fill: colour::GAL,
                label: "GalA".into(),
            });
        }
        "2112" | "2111" if has_n_mod => {
            if has_nac {
                return Ok(Symbol {
                    shape: Shape::Square,
                    fill: colour::GAL,
                    label: "GalNAc".into(),
                });
            } else {
                return Ok(Symbol {
                    shape: Shape::NSquare,
                    fill: colour::GAL,
                    label: "GalN".into(),
                });
            }
        }
        "2112" | "2111" => {
            return Ok(Symbol {
                shape: Shape::Circle,
                fill: colour::GAL,
                label: "Gal".into(),
            });
        }
        "1221" if has_acid => {
            return Ok(Symbol {
                shape: Shape::SplitDiamondTop,
                fill: colour::MAN,
                label: "ManA".into(),
            });
        }
        "1221" if has_deoxy => {
            return Ok(Symbol {
                shape: Shape::Triangle,
                fill: colour::FUC,
                label: "Fuc".into(),
            });
        }
        "1221" if has_n_mod => {
            if has_nac {
                return Ok(Symbol {
                    shape: Shape::Square,
                    fill: colour::MAN,
                    label: "ManNAc".into(),
                });
            } else {
                return Ok(Symbol {
                    shape: Shape::NSquare,
                    fill: colour::MAN,
                    label: "ManN".into(),
                });
            }
        }
        "1221" => {
            return Ok(Symbol {
                shape: Shape::Circle,
                fill: colour::MAN,
                label: "Man".into(),
            });
        }
        "1122" if has_acid => {
            return Ok(Symbol {
                shape: Shape::SplitDiamondTop,
                fill: colour::MAN,
                label: "ManA".into(),
            });
        }
        "1122" if has_n_mod => {
            if has_nac {
                return Ok(Symbol {
                    shape: Shape::Square,
                    fill: colour::MAN,
                    label: "ManNAc".into(),
                });
            } else {
                return Ok(Symbol {
                    shape: Shape::NSquare,
                    fill: colour::MAN,
                    label: "ManN".into(),
                });
            }
        }
        "1122" => {
            return Ok(Symbol {
                shape: Shape::Circle,
                fill: colour::MAN,
                label: "Man".into(),
            });
        }
        "2211" if has_acid => {
            return Ok(Symbol {
                shape: Shape::SplitDiamondTop,
                fill: colour::MAN,
                label: "ManA".into(),
            });
        }
        "2211" if has_deoxy => {
            return Ok(Symbol {
                shape: Shape::Triangle,
                fill: colour::RHA,
                label: "Rha".into(),
            });
        }
        "2211" if has_n_mod => {
            if has_nac {
                return Ok(Symbol {
                    shape: Shape::Square,
                    fill: colour::MAN,
                    label: "ManNAc".into(),
                });
            } else {
                return Ok(Symbol {
                    shape: Shape::NSquare,
                    fill: colour::MAN,
                    label: "ManN".into(),
                });
            }
        }
        "2211" => {
            return Ok(Symbol {
                shape: Shape::Circle,
                fill: colour::MAN,
                label: "Man".into(),
            });
        }
        "1121" => {
            return Ok(Symbol {
                shape: Shape::Circle,
                fill: colour::GUL,
                label: "Gul".into(),
            });
        }
        "2222" => {
            return Ok(Symbol {
                shape: Shape::Circle,
                fill: colour::ALL,
                label: "All".into(),
            });
        }
        "2221" => {
            return Ok(Symbol {
                shape: Shape::Circle,
                fill: colour::TAL,
                label: "Tal".into(),
            });
        }
        _ if bare_str.contains('d') && bare_str.len() <= 5 && has_deoxy => {
            return Ok(Symbol {
                shape: Shape::Triangle,
                fill: colour::FUC,
                label: "Fuc".into(),
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
                label: "Neu5Gc".into(),
            });
        }
        if has_nac {
            return Ok(Symbol {
                shape: Shape::Diamond,
                fill: colour::NEU5AC,
                label: "Neu5Ac".into(),
            });
        }
        if residue.modifications.iter().any(|m| m.descriptor == "O") {
            return Ok(Symbol {
                shape: Shape::Diamond,
                fill: colour::KDN,
                label: "KDN".into(),
            });
        }
        return Ok(Symbol {
            shape: Shape::Diamond,
            fill: colour::NEU5AC,
            label: "Sia".into(),
        });
    }

    // ── KDO ────────────────────────────────────────────────────────────
    if bare_str.len() == 4 && bare_str.contains("1122") && has_acid {
        return Ok(Symbol {
            shape: Shape::Diamond,
            fill: colour::KDO,
            label: "KDO".into(),
        });
    }

    // ── Pentoses ───────────────────────────────────────────────────────
    // Handle specific pentose types before generic fallback
    match bare_str {
        "211" | "122" => {
            // Ara (arabinose) - use green for furanose forms
            return Ok(Symbol {
                shape: Shape::Star,
                fill: colour::MAN, // green - same as other furanoses
                label: "Ara".into(),
            });
        }
        "212" => {
            // Xyl (xylose)
            return Ok(Symbol {
                shape: Shape::Star,
                fill: colour::XYL,
                label: "Xyl".into(),
            });
        }
        "112" | "221" => {
            // Lyx and its absolute-configuration inverse use the same family
            // colour in SNFG.
            return Ok(Symbol {
                shape: Shape::Star,
                fill: colour::YELLOW,
                label: "Lyx".into(),
            });
        }
        "111" | "222" => {
            // Rib and its absolute-configuration inverse.
            return Ok(Symbol {
                shape: Shape::Star,
                fill: colour::PINK,
                label: "Rib".into(),
            });
        }
        "121" => {
            // Other pentose (Sor, etc.)
            return Ok(Symbol {
                shape: Shape::Star,
                fill: colour::MAN, // green - matches reference
                label: "?".into(),
            });
        }
        _ if bare_str.len() == 3 => {
            // Generic pentose fallback
            return Ok(Symbol {
                shape: Shape::Star,
                fill: colour::MAN, // green - matches reference
                label: "Xyl".into(),
            });
        }
        _ => {}
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
            label: if has_nac { "HexNAc" } else { "Hex" }.into(),
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
        label: "?".into(),
    })
}

/// The configuration letter (`D` or `L`) SNFG asks to be shown inside the
/// symbol of a residue whose absolute configuration is not the default of its
/// name (D-Fuc, D-Rha, L-Man, D-Ara, L-Gul, ...); `None` otherwise.
fn configuration_mark(residue: &Monosaccharide) -> Option<char> {
    let (kind, inverted) = classify_residue_configuration(residue)?;
    if !inverted {
        return None;
    }
    kind.default_configuration()
        .map(|default| if default == 'L' { 'D' } else { 'L' })
}

/// The style rule of configuration letters, emitted only for a glycan that
/// has one so that every other drawing is unchanged.
fn configuration_label_css(graph: &ResidueGraph, font_family: &str, scale: f64) -> String {
    if !graph
        .inner()
        .node_weights()
        .any(|residue| configuration_mark(residue).is_some())
    {
        return String::new();
    }
    format!(
        "  .config-label {{ font-family: {font_family}; font-size: {}px; font-weight: bold; fill: #000; text-anchor: middle; dominant-baseline: central; }}\n",
        16.0 * scale
    )
}

fn draw_configuration_mark(svg: &mut String, residue: &Monosaccharide, cx: f64, cy: f64) {
    if let Some(letter) = configuration_mark(residue) {
        svg.push_str(&format!(
            "<text x=\"{cx}\" y=\"{cy}\" class=\"config-label\">{letter}</text>\n"
        ));
    }
}

// ── Geometry constants ─────────────────────────────────────────────────────

pub const NODE_R: f64 = 25.0;
pub const H_SPACING: f64 = 100.0;
pub const V_SPACING: f64 = 100.0;
pub const BOND_W: f64 = 4.0;
pub const LABEL_SIZE: f64 = 20.0;
pub const PNG_SCALE: u32 = 2;
/// Share of the usual residue spacing used by [`RenderOptions::compact`].
const COMPACT_SPACING: f64 = 0.7;

// ── Options ────────────────────────────────────────────────────────────────

/// Where the reducing end of a glycan is drawn, and with it the direction the
/// glycan grows in. The higher acceptor position of a branch stays above the
/// lower one in the horizontal orientations and to its right in the vertical
/// ones (an N-glycan's α1-6 arm is drawn above, or to the right of, α1-3).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Orientation {
    /// Reducing end on the right, the SNFG convention.
    #[default]
    RootRight,
    /// Reducing end on the left, growing to the right.
    RootLeft,
    /// Reducing end at the bottom, growing upwards.
    RootBottom,
    /// Reducing end at the top, growing downwards.
    RootTop,
}

impl Orientation {
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::RootRight => "root-right",
            Self::RootLeft => "root-left",
            Self::RootBottom => "root-bottom",
            Self::RootTop => "root-top",
        }
    }

    pub const fn is_vertical(self) -> bool {
        matches!(self, Self::RootBottom | Self::RootTop)
    }
}

impl FromStr for Orientation {
    type Err = String;

    fn from_str(value: &str) -> Result<Self, Self::Err> {
        match value.trim().to_ascii_lowercase().replace('_', "-").as_str() {
            "root-right" | "right-to-left" | "rl" | "horizontal" => Ok(Self::RootRight),
            "root-left" | "left-to-right" | "lr" => Ok(Self::RootLeft),
            "root-bottom" | "bottom-up" | "bt" | "vertical" => Ok(Self::RootBottom),
            "root-top" | "top-down" | "tb" => Ok(Self::RootTop),
            other => Err(format!(
                "unknown orientation: {other} (expected root-right, root-left, root-bottom or root-top)"
            )),
        }
    }
}

#[derive(Debug, Clone)]
pub struct SourceNotation {
    pub format: String,
    pub value: String,
}

impl SourceNotation {
    pub fn new(format: impl Into<String>, value: impl Into<String>) -> Self {
        Self {
            format: format.into(),
            value: value.into(),
        }
    }
}

#[derive(Debug, Clone)]
pub struct RenderOptions {
    pub colour: bool,
    pub show_labels: bool,   // show residue abbreviations inside shapes
    pub show_linkages: bool, // show linkage positions on bonds
    pub font_family: String,
    pub scale: f64,
    /// Exact source notation supplied by the caller. When absent, metadata
    /// falls back to source text retained by the parsed graph.
    pub source_notation: Option<SourceNotation>,
    /// Where the reducing end is drawn.
    pub orientation: Orientation,
    /// Tighter spacing and no linkage text, for heatmap axes and tables.
    pub compact: bool,
    /// Text drawn past the reducing end (`Asn`, `Ser/Thr`, `-OH`, `2-AB`).
    pub reducing_end: Option<String>,
    /// Draw repeat units in square brackets with their repeat count instead
    /// of a dashed bond closing the unit.
    pub repeat_brackets: bool,
    /// Draw undefined antennae and substituents above the structure over a
    /// bracket across their candidate parents (the SNFG convention for a
    /// fragment of unknown attachment) instead of dashed bonds to each one.
    pub floating_fragments: bool,
    /// Fragment-ion cleavages to mark on the bonds.
    pub cleavages: Vec<CleavageMark>,
}

impl Default for RenderOptions {
    fn default() -> Self {
        Self {
            colour: true,
            show_labels: false, // SNFG convention: shape + colour = identity
            show_linkages: true,
            font_family: "Arial, Helvetica, sans-serif".into(),
            scale: 1.0,
            source_notation: None,
            orientation: Orientation::default(),
            compact: false,
            reducing_end: None,
            repeat_brackets: false,
            floating_fragments: false,
            cleavages: Vec::new(),
        }
    }
}

/// The fragment a cleavage keeps, which decides the side its mark turns to.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CleavageSide {
    /// B and C ions keep the non-reducing end.
    NonReducing,
    /// Y and Z ions keep the reducing end.
    Reducing,
}

/// A Domon–Costello cleavage mark on one bond of the drawing.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CleavageMark {
    /// Edge index of the cleaved bond in the `ResidueGraph`.
    pub edge: usize,
    pub side: CleavageSide,
    /// Ion label (`B2`, `Y3α`).
    pub label: String,
}

// ── Tree layout ────────────────────────────────────────────────────────────

/// `x` counts residues from the reducing end (fragments of a floating band
/// sit between two of them); `y` is the lane, in layout pixels.
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
    layout_components(graph, root, &mut info, &mut visited);
    info
}

/// Lay out `root`'s component and every component not already `visited`,
/// and centre them around y = 0.
fn layout_components(
    graph: &ResidueGraph,
    root: NodeIndex,
    info: &mut HashMap<usize, LayoutInfo>,
    visited: &mut std::collections::HashSet<usize>,
) {
    let mut next_leaf = 0usize;
    layout_subtree(graph, root, 0, &mut next_leaf, info, visited);
    // WURCS compositions and undefined antennae can contain disconnected
    // components. Lay every component out instead of silently omitting it.
    for node in graph.inner().node_indices() {
        if !visited.contains(&node.index()) {
            if next_leaf > 0 {
                next_leaf += 1;
            }
            layout_subtree(graph, node, 0, &mut next_leaf, info, visited);
        }
    }

    resolve_triangle_collisions(graph, info);

    // centre around y=0
    let min_y = info.values().map(|li| li.y).fold(f64::MAX, f64::min);
    let max_y = info.values().map(|li| li.y).fold(f64::MIN, f64::max);
    let shift = -(min_y + max_y) / 2.0;
    for li in info.values_mut() {
        li.y += shift;
    }
}

/// One thing of unknown attachment in a floating band.
#[derive(Debug, Clone)]
enum FloatingItem {
    /// An undefined fragment: the node index of its reducing-end residue.
    Fragment(usize),
    /// An undefined substituent, by its MAP code and the lane of its label.
    Substituent(String, f64),
}

/// Fragments of unknown attachment that share their candidate parents, drawn
/// above the structure over one bracket.
#[derive(Debug, Clone)]
struct FloatingBand {
    /// Lane of the bracket.
    lane: f64,
    /// Depth of every candidate parent, in ascending order.
    depths: Vec<f64>,
    /// Depth of the stem that joins the items to the bracket.
    stem: f64,
    items: Vec<FloatingItem>,
}

/// Gap between a bracket and whatever is drawn next to it, in layout pixels.
const FLOATING_GAP: f64 = 60.0;
/// How far from the stem a floating fragment starts, in residues.
const FLOATING_REACH: f64 = 0.7;

/// The layout of a glycan whose undefined fragments float above the main
/// structure. Fragments attached to several candidates are grouped by their
/// candidates; everything else is laid out as usual.
fn compute_floating_layout(
    graph: &ResidueGraph,
    root: NodeIndex,
) -> (HashMap<usize, LayoutInfo>, Vec<FloatingBand>) {
    let mut fragments = std::collections::HashSet::new();
    let main = graph.subtree(root);
    for undefined in graph.undefined_linkages() {
        if !main.contains(&undefined.child.index()) {
            fragments.extend(graph.subtree(undefined.child));
        }
    }

    let mut info = HashMap::new();
    // A fragment is laid out on its own below; marking it visited keeps it
    // out of the main structure's lanes.
    let mut visited = fragments.clone();
    layout_components(graph, root, &mut info, &mut visited);

    // Group by candidate parents, keeping the order of first appearance.
    let mut groups: Vec<(BTreeSet<usize>, Vec<FloatingItem>)> = Vec::new();
    let mut add = |parents: BTreeSet<usize>, item: FloatingItem| match groups
        .iter_mut()
        .find(|(known, _)| *known == parents)
    {
        Some((_, items)) => items.push(item),
        None => groups.push((parents, vec![item])),
    };
    for undefined in graph.undefined_linkages() {
        let child = undefined.child.index();
        if fragments.contains(&child) {
            let parents = undefined.parents.iter().map(|p| p.residue.index());
            add(parents.collect(), FloatingItem::Fragment(child));
        }
    }
    for modification in graph.undefined_modifications() {
        let parents = modification.parents.iter().map(|p| p.residue.index());
        add(
            parents.collect(),
            FloatingItem::Substituent(modification.map_code.clone(), 0.0),
        );
    }

    // The main structure is centred on lane 0, so its top is at or above it.
    let mut cursor = info.values().map(|li| li.y).fold(0.0, f64::min);
    let mut bands = Vec::new();
    for (parents, items) in groups {
        let mut depths = parents
            .iter()
            .filter_map(|parent| info.get(parent).map(|li| li.x))
            .collect::<Vec<_>>();
        depths.sort_by(f64::total_cmp);
        depths.dedup();
        let (first, last) = match (depths.first(), depths.last()) {
            (Some(first), Some(last)) => (*first, *last),
            _ => (0.0, 0.0),
        };
        let stem = (first + last) / 2.0;
        let lane = cursor - FLOATING_GAP;
        // Lowest free lane of the band; items stack upwards from here.
        let mut base = lane - FLOATING_GAP;
        let mut band_items = Vec::with_capacity(items.len());
        for item in items {
            match item {
                FloatingItem::Fragment(child) => {
                    let mut fragment = HashMap::new();
                    let mut seen = std::collections::HashSet::new();
                    let mut next_leaf = 0usize;
                    layout_subtree(
                        graph,
                        NodeIndex::new(child),
                        0,
                        &mut next_leaf,
                        &mut fragment,
                        &mut seen,
                    );
                    resolve_triangle_collisions(graph, &mut fragment);
                    let lowest = fragment.values().map(|li| li.y).fold(f64::MIN, f64::max);
                    let highest = fragment.values().map(|li| li.y).fold(f64::MAX, f64::min);
                    let shift = base - lowest;
                    for (node, position) in fragment {
                        info.insert(
                            node,
                            LayoutInfo {
                                x: position.x + stem + FLOATING_REACH,
                                y: position.y + shift,
                            },
                        );
                    }
                    base += highest - lowest - V_SPACING;
                    band_items.push(FloatingItem::Fragment(child));
                }
                FloatingItem::Substituent(map_code, _) => {
                    band_items.push(FloatingItem::Substituent(map_code, base));
                    base -= V_SPACING * 0.6;
                }
            }
        }
        cursor = base + V_SPACING * 0.6;
        bands.push(FloatingBand {
            lane,
            depths,
            stem,
            items: band_items,
        });
    }
    (info, bands)
}

fn resolve_triangle_collisions(graph: &ResidueGraph, info: &mut HashMap<usize, LayoutInfo>) {
    let mut branches = graph
        .inner()
        .edge_references()
        .filter(|edge| {
            is_fucose(graph, edge.target())
                && is_terminal(graph, edge.target())
                && info.contains_key(&edge.source().index())
                && info.contains_key(&edge.target().index())
        })
        .map(|edge| {
            (
                edge.source(),
                edge.target(),
                edge.weight().parent_position.0,
            )
        })
        .collect::<Vec<_>>();
    branches.sort_by(|(left_parent, _, _), (right_parent, _, _)| {
        info[&left_parent.index()]
            .y
            .total_cmp(&info[&right_parent.index()].y)
    });

    for (parent, fucose, linkage_pos) in branches {
        let parent_layout = info[&parent.index()].clone();

        // Check if this parent has both α3 and α6 triangle children (Fuc) or α2 and α4 triangle children (Rha)
        let parent_fucose_children: Vec<u8> = graph
            .inner()
            .edges_directed(parent, Direction::Outgoing)
            .filter(|edge| is_fucose(graph, edge.target()) && is_terminal(graph, edge.target()))
            .map(|edge| edge.weight().parent_position.0)
            .collect();

        let has_both_alpha3_and_alpha6 =
            parent_fucose_children.contains(&3) && parent_fucose_children.contains(&6);
        let has_both_alpha2_and_alpha4 =
            parent_fucose_children.contains(&2) && parent_fucose_children.contains(&4);

        // Check if this is core fucose - attached to the root GlcNAc (reducing end)
        let is_core_fucose = linkage_pos == 6
            && graph.residue(parent).is_some_and(|res| {
                let skel = &res.skeleton_code;
                let bare: String = skel.chars().take_while(|c| c.is_ascii_digit()).collect();
                // Check if this is a GlcNAc (2122 with N-acetyl) AND is the root node
                bare == "2122"
                    && res
                        .modifications
                        .iter()
                        .any(|m| m.descriptor.contains("NCC"))
                    && parent == graph.root().unwrap()
            });

        // Use the same positioning logic as layout_subtree for consistency
        let desired_y = if has_both_alpha3_and_alpha6 {
            parent_layout.y
                + if linkage_pos == 6 {
                    -V_SPACING // α6 fucose goes UP when paired with α3
                } else {
                    V_SPACING // α3 and other positions go DOWN
                }
        } else if has_both_alpha2_and_alpha4 {
            parent_layout.y
                + if linkage_pos == 4 {
                    -V_SPACING // α4 rhamnose goes UP when paired with α2
                } else {
                    V_SPACING // α2 and other positions go DOWN
                }
        } else if is_core_fucose {
            parent_layout.y + -V_SPACING // Core α6 fucose defaults to UP
        } else {
            parent_layout.y + V_SPACING // Single triangle defaults to DOWN
        };

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

        // SNFG convention draws terminal fucose/rhamnose vertically aligned with parent.
        // To prevent overlap when multiple triangles are attached to the same parent,
        // position them in opposite vertical directions when both are present.
        let has_both_alpha3_and_alpha6 = fucose_children.iter().any(|(_, pos)| *pos == 3)
            && fucose_children.iter().any(|(_, pos)| *pos == 6);
        let has_both_alpha2_and_alpha4 = fucose_children.iter().any(|(_, pos)| *pos == 2)
            && fucose_children.iter().any(|(_, pos)| *pos == 4);

        for (child, linkage_pos) in fucose_children.into_iter() {
            visited.insert(child.index());
            // Check if this is core fucose - attached to the root GlcNAc (reducing end)
            let is_core_fucose = linkage_pos == 6
                && graph.residue(node).is_some_and(|res| {
                    let skel = &res.skeleton_code;
                    let bare: String = skel.chars().take_while(|c| c.is_ascii_digit()).collect();
                    // Check if this is a GlcNAc (2122 with N-acetyl) AND is the root node
                    bare == "2122"
                        && res
                            .modifications
                            .iter()
                            .any(|m| m.descriptor.contains("NCC"))
                        && node == graph.root().unwrap()
                });

            let vertical_offset = if has_both_alpha3_and_alpha6 {
                if linkage_pos == 6 {
                    -V_SPACING // α6 fucose goes UP when paired with α3
                } else {
                    V_SPACING // α3 and other positions go DOWN
                }
            } else if has_both_alpha2_and_alpha4 {
                if linkage_pos == 4 {
                    -V_SPACING // α4 rhamnose goes UP when paired with α2
                } else {
                    V_SPACING // α2 and other positions go DOWN
                }
            } else if is_core_fucose {
                -V_SPACING // Core α6 fucose defaults to UP
            } else {
                V_SPACING // Single triangle defaults to DOWN
            };
            info.insert(
                child.index(),
                LayoutInfo {
                    x: depth as f64,
                    y: y + vertical_offset,
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
    render_svg_internal(graph, opts, None)
}

/// Exact graph elements to keep vivid in a highlighted SNFG rendering.
/// Indices are the stable `petgraph` node/edge indices exposed by the source
/// `ResidueGraph`; elements not listed here are rendered with the muted SNFG
/// palette.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct HighlightSelection {
    pub node_indices: BTreeSet<usize>,
    pub edge_indices: BTreeSet<usize>,
}

impl HighlightSelection {
    /// The union of motif occurrences, for example of a positional
    /// [`MotifQuery`](crabwurcs_core::MotifQuery).
    pub fn from_matches<'a>(matches: impl IntoIterator<Item = &'a MotifMatch>) -> Self {
        let mut selection = Self::default();
        for found in matches {
            selection.node_indices.extend(&found.node_indices);
            selection.edge_indices.extend(&found.edge_indices);
        }
        selection
    }

    /// The given residues and the bonds between them, for residues known by
    /// their node index (PDB provenance, a fragment ion).
    pub fn from_residues(
        graph: &ResidueGraph,
        residues: impl IntoIterator<Item = usize>,
    ) -> SnfgResult<Self> {
        let count = graph.node_count();
        let node_indices = residues.into_iter().collect::<BTreeSet<_>>();
        if let Some(index) = node_indices.iter().find(|index| **index >= count) {
            return Err(SnfgError::InvalidResidueIndex {
                index: *index,
                count,
            });
        }
        let edge_indices = graph
            .inner()
            .edge_references()
            .filter(|edge| {
                node_indices.contains(&edge.source().index())
                    && node_indices.contains(&edge.target().index())
            })
            .map(|edge| edge.id().index())
            .collect();
        Ok(Self {
            node_indices,
            edge_indices,
        })
    }
}

/// Render an SNFG SVG with exactly the supplied graph nodes and edges
/// highlighted. This is intended for callers that have already resolved a
/// specific occurrence and must not highlight every matching motif.
pub fn render_svg_with_selection(
    graph: &ResidueGraph,
    selection: &HighlightSelection,
    opts: &RenderOptions,
) -> SnfgResult<String> {
    render_svg_internal(graph, opts, Some(selection))
}

/// Render an SNFG SVG while highlighting the union of every occurrence of
/// every supplied motif.
///
/// Matched nodes and motif-internal edges remain fully opaque. Everything
/// else uses GlycoDraw's muted SNFG palette while remaining fully opaque. This
/// prevents bonds drawn behind unmatched residues from showing through their
/// symbols. An empty motif list is identical to [`render_svg_with_options`].
pub fn render_svg_with_motifs(
    graph: &ResidueGraph,
    motifs: &[ResidueGraph],
    opts: &RenderOptions,
) -> SnfgResult<String> {
    if motifs.is_empty() {
        return render_svg_with_options(graph, opts);
    }
    let mut selection = HighlightSelection::default();
    for motif in motifs {
        for found in find_motif_matches(graph, motif)? {
            selection.node_indices.extend(found.node_indices);
            selection.edge_indices.extend(found.edge_indices);
        }
    }
    render_svg_with_selection(graph, &selection, opts)
}

/// Render an SNFG diagram as a transparent PNG at twice the SVG dimensions.
pub fn render_png(graph: &ResidueGraph) -> SnfgResult<Vec<u8>> {
    render_png_with_options(graph, &RenderOptions::default())
}

/// Render an SNFG diagram as a transparent PNG using the supplied SVG
/// rendering options.
pub fn render_png_with_options(graph: &ResidueGraph, opts: &RenderOptions) -> SnfgResult<Vec<u8>> {
    svg_to_png(&render_svg_with_options(graph, opts)?)
}

/// Render an exact-selection SNFG diagram as a transparent PNG at twice the
/// SVG dimensions.
pub fn render_png_with_selection(
    graph: &ResidueGraph,
    selection: &HighlightSelection,
    opts: &RenderOptions,
) -> SnfgResult<Vec<u8>> {
    svg_to_png(&render_svg_with_selection(graph, selection, opts)?)
}

/// Render a motif-highlighted SNFG diagram as a transparent PNG at twice the
/// SVG dimensions.
pub fn render_png_with_motifs(
    graph: &ResidueGraph,
    motifs: &[ResidueGraph],
    opts: &RenderOptions,
) -> SnfgResult<Vec<u8>> {
    svg_to_png(&render_svg_with_motifs(graph, motifs, opts)?)
}

/// Render a single SNFG symbol as a tightly-cropped, transparent-background
/// SVG.
///
/// Convenient for building legends, icons, or thumbnail galleries — for
/// example the per-residue symbol table in this crate's documentation. The
/// symbol is selected with the same [`symbol_for`] classification used by the
/// full graph renderer, so a standalone symbol always matches what
/// [`render_svg`] would draw for a one-residue graph. Pass `show_labels: true`
/// in [`RenderOptions`] to draw the residue abbreviation inside the shape.
pub fn render_symbol_svg(residue: &Monosaccharide, opts: &RenderOptions) -> SnfgResult<String> {
    let symbol = symbol_for(residue)?;
    let r = NODE_R * opts.scale;
    // Star and diamond outer vertices reach roughly 1.3 × r; leave room for
    // the stroke so the shapes are not clipped at the canvas edge.
    let extent = r * 1.3 + 4.0 * opts.scale;
    let side = extent * 2.0;
    let cx = side / 2.0;
    let cy = side / 2.0;
    let mut svg = format!(
        r#"<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" viewBox="0 0 {w} {h}">
"#,
        w = side,
        h = side,
    );
    draw_shape(&mut svg, &symbol, cx, cy, r, opts);
    if opts.show_labels {
        let x = cx;
        let y = cy;
        let ff = opts.font_family.as_str();
        let label = escape_xml_text(&symbol.label);
        svg.push_str(&format!(
            r##"<text x="{x}" y="{y}" font-family="{ff}" font-size="11px" fill="#000" text-anchor="middle" dominant-baseline="central">{label}</text>
"##,
        ));
    } else if let Some(letter) = configuration_mark(residue) {
        let ff = opts.font_family.as_str();
        let size = 16.0 * opts.scale;
        svg.push_str(&format!(
            r##"<text x="{cx}" y="{cy}" font-family="{ff}" font-size="{size}px" font-weight="bold" fill="#000" text-anchor="middle" dominant-baseline="central">{letter}</text>
"##,
        ));
    }
    svg.push_str("</svg>\n");
    Ok(svg)
}

fn svg_to_png(svg: &str) -> SnfgResult<Vec<u8>> {
    let mut options = resvg::usvg::Options::default();
    options.fontdb_mut().load_system_fonts();
    let tree = resvg::usvg::Tree::from_str(svg, &options)
        .map_err(|error| SnfgError::SvgParse(error.to_string()))?;
    let svg_size = tree.size().to_int_size();
    let width = svg_size
        .width()
        .checked_mul(PNG_SCALE)
        .ok_or(SnfgError::RasterDimensions)?;
    let height = svg_size
        .height()
        .checked_mul(PNG_SCALE)
        .ok_or(SnfgError::RasterDimensions)?;
    let mut pixmap =
        resvg::tiny_skia::Pixmap::new(width, height).ok_or(SnfgError::RasterAllocation)?;
    let transform = resvg::tiny_skia::Transform::from_scale(PNG_SCALE as f32, PNG_SCALE as f32);
    resvg::render(&tree, transform, &mut pixmap.as_mut());
    pixmap
        .encode_png()
        .map_err(|error| SnfgError::PngEncoding(error.to_string()))
}

fn highlight_class(selected: bool) -> &'static str {
    if selected {
        "motif-match"
    } else {
        "motif-dimmed"
    }
}

const MOTIF_DEEMPHASIS_CSS: &str = r##"  .motif-match { opacity: 1; }
  .motif-dimmed [fill="#0072BC"] { fill: #CDE7EF; }
  .motif-dimmed [fill="#00A651"] { fill: #CDE9DF; }
  .motif-dimmed [fill="#FFD400"] { fill: #FFF6DE; }
  .motif-dimmed [fill="#F47920"] { fill: #FDE7E0; }
  .motif-dimmed [fill="#F69EA1"] { fill: #FDF0F1; }
  .motif-dimmed [fill="#A54399"] { fill: #F1E6ED; }
  .motif-dimmed [fill="#8FCCE9"] { fill: #EEF8FB; }
  .motif-dimmed [fill="#A17A4D"] { fill: #F1E9E5; }
  .motif-dimmed [fill="#ED1C24"] { fill: #F7E0E0; }
  .motif-dimmed [stroke="#000000"],
  .motif-dimmed line { stroke: #D9D9D9; }
  .motif-dimmed text { fill: #D9D9D9; }
"##;

fn render_svg_internal(
    graph: &ResidueGraph,
    opts: &RenderOptions,
    highlights: Option<&HighlightSelection>,
) -> SnfgResult<String> {
    let inner = graph.inner();
    if inner.node_count() == 0 {
        return Ok(empty_svg(graph, opts));
    }
    if graph.is_composition() {
        return render_composition_svg(graph, opts, highlights);
    }

    let root = graph.root().unwrap_or_else(|| NodeIndex::from(0u32));
    let float = opts.floating_fragments
        && !(graph.undefined_linkages().is_empty() && graph.undefined_modifications().is_empty());
    let (layout, bands) = if float {
        compute_floating_layout(graph, root)
    } else {
        (compute_layout(graph, root), Vec::new())
    };

    let repeats = if opts.repeat_brackets {
        repeat_brackets(graph, &layout)
    } else {
        Vec::new()
    };
    let lanes = layout
        .values()
        .map(|li| li.y)
        .chain(bands.iter().map(|band| band.lane));
    let min_y = lanes.clone().fold(f64::MAX, f64::min);
    let max_y = lanes.fold(f64::MIN, f64::max);
    let label_reach = bands
        .iter()
        .filter(|band| {
            band.items
                .iter()
                .any(|item| matches!(item, FloatingItem::Substituent(..)))
        })
        .map(|band| band.stem + FLOATING_REACH + 0.45);
    let max_depth = layout
        .values()
        .map(|li| li.x)
        .chain(label_reach)
        .fold(0.0f64, f64::max);

    let s = opts.scale;
    let vertical = opts.orientation.is_vertical();
    let spacing = if opts.compact { COMPACT_SPACING } else { 1.0 };
    let has_modification_label = |residue: &Monosaccharide| {
        let labels = build_modification_labels(residue);
        !labels.above.is_empty() || !labels.below.is_empty()
    };
    let margin = if !opts.compact {
        70.0
    } else if inner.node_weights().any(has_modification_label) {
        52.0
    } else {
        36.0
    };
    let mut pad = if graph.undefined_modifications().is_empty() || float {
        margin * s
    } else {
        110.0 * s
    };
    if !repeats.is_empty() {
        pad = pad.max(90.0 * s);
    }
    if let Some(longest) = opts.cleavages.iter().map(|m| m.label.chars().count()).max() {
        let label = if vertical {
            longest as f64 * CLEAVAGE_LABEL_SIZE * 0.62
        } else {
            CLEAVAGE_LABEL_SIZE
        };
        pad = pad.max((CLEAVAGE_REACH + 14.0 + label + 6.0) * s);
    }

    // Past the root the canvas needs room for the reducing-end text.
    let reducing_end = opts
        .reducing_end
        .as_deref()
        .map(str::trim)
        .filter(|text| !text.is_empty());
    let text_width = reducing_end.map_or(0.0, |text| {
        text.chars().count() as f64 * LABEL_SIZE * s * 0.62
    });
    let mut root_extra = 0.0;
    if reducing_end.is_some() {
        let reach =
            REDUCING_END_STUB * s + 8.0 * s + if vertical { LABEL_SIZE * s } else { text_width };
        root_extra = (reach + 8.0 * s - pad).max(0.0);
        if vertical {
            pad = pad.max(text_width / 2.0 + 8.0 * s);
        }
    }

    let depth_extent = max_depth * (H_SPACING * spacing) * s;
    let lane_scale = spacing * s;
    let lane_extent = (max_y - min_y) * lane_scale;
    let (canvas_w, canvas_h) = if vertical {
        (
            lane_extent + 2.0 * pad,
            depth_extent + 2.0 * pad + root_extra,
        )
    } else {
        (
            depth_extent + 2.0 * pad + root_extra,
            lane_extent + 2.0 * pad,
        )
    };
    let frame = Frame {
        orientation: opts.orientation,
        width: canvas_w,
        height: canvas_h,
        root_pad: pad + root_extra,
        depth_step: H_SPACING * spacing,
        scale: s,
        lane_scale,
        lane_offset: -min_y * lane_scale + pad,
    };
    let to_canvas = |li: &LayoutInfo| frame.at(li);

    // ── SVG header ────────────────────────────────────────────────────
    let highlight_css = if highlights.is_some() {
        MOTIF_DEEMPHASIS_CSS
    } else {
        ""
    };
    let metadata = svg_metadata(graph, opts);
    let mut extra_css = String::new();
    if !repeats.is_empty() || !bands.is_empty() {
        extra_css.push_str(&format!(
            "  .bracket {{ fill: none; stroke-width: {}; stroke-linecap: square; }}\n",
            2.5 * s
        ));
    }
    if !repeats.is_empty() {
        extra_css.push_str(&format!(
            "  .repeat-count {{ font-family: {}; font-size: {}px; font-style: italic; fill: #000; text-anchor: middle; dominant-baseline: central; }}\n",
            opts.font_family,
            LABEL_SIZE * s
        ));
    }
    if reducing_end.is_some() {
        extra_css.push_str(&format!(
            "  .reducing-end {{ font-family: {}; font-size: {}px; fill: #000; dominant-baseline: central; }}\n",
            opts.font_family,
            LABEL_SIZE * s
        ));
    }
    if !opts.cleavages.is_empty() {
        extra_css.push_str(&format!(
            "  .cleavage {{ fill: none; stroke: {CLEAVAGE_COLOUR}; stroke-width: {}; stroke-linejoin: round; stroke-linecap: round; }}\n  .cleavage-label {{ font-family: {}; font-size: {}px; fill: {CLEAVAGE_COLOUR}; dominant-baseline: central; }}\n",
            2.5 * s,
            opts.font_family,
            CLEAVAGE_LABEL_SIZE * s
        ));
    }
    let mut svg = format!(
        r#"<svg xmlns="http://www.w3.org/2000/svg" role="img" aria-labelledby="snfg-title snfg-desc" viewBox="0 0 {w} {h}" width="{w}" height="{h}">
{metadata}<style>
  .bond {{ stroke: #000; stroke-width: {bw}; fill: none; stroke-linecap: round; }}
  .uncertain {{ stroke: #555; stroke-width: {ubw}; fill: none; stroke-linecap: round; stroke-dasharray: 8 7; }}
  .link {{ font-family: {ff}; font-size: {ls}px; fill: #000; text-anchor: middle; dominant-baseline: central; }}
  .mod-label {{ font-family: {ff}; font-size: {ls}px; fill: #000; text-anchor: middle; dominant-baseline: central; }}
  .node {{ fill: none; }}
  .res-label {{ font-family: {ff}; font-size: 11px; fill: #000; text-anchor: middle; dominant-baseline: central; }}
{config_css}{extra_css}{highlight_css}</style>
"#,
        config_css = configuration_label_css(graph, &opts.font_family, s),
        w = canvas_w,
        h = canvas_h,
        bw = BOND_W * s,
        ubw = 2.5 * s,
        ff = opts.font_family,
        ls = LABEL_SIZE * s,
        highlight_css = highlight_css,
        metadata = metadata,
    );

    let linkage_text = opts.show_linkages && !opts.compact;

    // ── Edges and linkage labels ──────────────────────────────────────
    for edge in inner.edge_references() {
        if opts.repeat_brackets && edge.weight().repeat.is_some() {
            continue;
        }
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
        if let Some(selection) = highlights {
            svg.push_str(&format!(
                r#"<g data-edge-index="{}" class="{}">
"#,
                edge.id().index(),
                highlight_class(selection.edge_indices.contains(&edge.id().index()))
            ));
        }
        svg.push_str(&format!(
            r#"<line x1="{px}" y1="{py}" x2="{cx}" y2="{cy}" class="{class}"/>
"#,
        ));
        if linkage_text {
            draw_bond_label(
                &mut svg,
                opts.orientation,
                px,
                py,
                cx,
                cy,
                &linkage_label_for(inner, edge.target(), edge.weight()),
                s,
            );
        }
        if highlights.is_some() {
            svg.push_str("</g>\n");
        }
    }

    for undefined in graph.undefined_linkages() {
        let floating = bands.iter().any(|band| {
            band.items.iter().any(|item| {
                matches!(item, FloatingItem::Fragment(child) if *child == undefined.child.index())
            })
        });
        if floating {
            continue;
        }
        let Some(child_layout) = layout.get(&undefined.child.index()) else {
            continue;
        };
        let (cx, cy) = to_canvas(child_layout);
        for (candidate_index, parent) in undefined.parents.iter().enumerate() {
            let Some(parent_layout) = layout.get(&parent.residue.index()) else {
                continue;
            };
            let (px, py) = to_canvas(parent_layout);
            if highlights.is_some() {
                svg.push_str(r#"<g class="motif-dimmed" data-undefined-linkage="true">"#);
                svg.push('\n');
            }
            svg.push_str(&format!(
                r#"<line x1="{px}" y1="{py}" x2="{cx}" y2="{cy}" class="uncertain"/>
"#,
            ));
            if linkage_text && candidate_index == 0 {
                draw_bond_label(&mut svg, opts.orientation, px, py, cx, cy, "?", s);
            }
            if highlights.is_some() {
                svg.push_str("</g>\n");
            }
        }
    }

    for band in &bands {
        draw_floating_band(
            &mut svg,
            band,
            graph,
            &layout,
            &frame,
            linkage_text,
            highlights.is_some(),
        );
    }

    let unbanded_modifications = if bands.is_empty() {
        graph.undefined_modifications()
    } else {
        &[]
    };
    for modification in unbanded_modifications {
        let candidates = modification
            .parents
            .iter()
            .filter_map(|parent| layout.get(&parent.residue.index()))
            .map(&to_canvas)
            .collect::<Vec<_>>();
        if candidates.is_empty() {
            continue;
        }
        let mean_x = || candidates.iter().map(|(x, _)| x).sum::<f64>() / candidates.len() as f64;
        let mean_y = || candidates.iter().map(|(_, y)| y).sum::<f64>() / candidates.len() as f64;
        let (label_x, label_y) = match opts.orientation {
            Orientation::RootRight => (
                candidates.iter().map(|(x, _)| *x).fold(f64::MIN, f64::max) + 65.0 * s,
                mean_y(),
            ),
            Orientation::RootLeft => (
                candidates.iter().map(|(x, _)| *x).fold(f64::MAX, f64::min) - 65.0 * s,
                mean_y(),
            ),
            Orientation::RootBottom => (
                mean_x(),
                candidates.iter().map(|(_, y)| *y).fold(f64::MIN, f64::max) + 65.0 * s,
            ),
            Orientation::RootTop => (
                mean_x(),
                candidates.iter().map(|(_, y)| *y).fold(f64::MAX, f64::min) - 65.0 * s,
            ),
        };
        if highlights.is_some() {
            svg.push_str(r#"<g class="motif-dimmed" data-undefined-modification="true">"#);
            svg.push('\n');
        }
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
        if highlights.is_some() {
            svg.push_str("</g>\n");
        }
    }

    for repeat in &repeats {
        draw_repeat_bracket(&mut svg, repeat, &frame, highlights.is_some());
    }

    if let (Some(text), Some(root_layout)) = (reducing_end, layout.get(&root.index())) {
        draw_reducing_end(
            &mut svg,
            text,
            root_layout,
            &frame,
            highlights.is_some_and(|selection| !selection.node_indices.contains(&root.index())),
        );
    }

    for mark in &opts.cleavages {
        let Some(edge) = inner
            .edge_references()
            .find(|edge| edge.id().index() == mark.edge)
        else {
            return Err(SnfgError::InvalidBondIndex {
                index: mark.edge,
                count: inner.edge_count(),
            });
        };
        if let (Some(parent), Some(child)) = (
            layout.get(&edge.source().index()),
            layout.get(&edge.target().index()),
        ) {
            draw_cleavage(
                &mut svg,
                mark,
                frame.at(parent),
                frame.at(child),
                linkage_text,
                s,
            );
        }
    }

    // ── Nodes ──────────────────────────────────────────────────────────
    for node_idx in inner.node_indices() {
        if let (Some(li), Some(residue)) =
            (layout.get(&node_idx.index()), inner.node_weight(node_idx))
        {
            let (cx, cy) = to_canvas(li);
            let symbol = symbol_for(residue)?;
            if let Some(selection) = highlights {
                svg.push_str(&format!(
                    r#"<g data-node-index="{}" class="{}">
"#,
                    node_idx.index(),
                    highlight_class(selection.node_indices.contains(&node_idx.index()))
                ));
            }
            draw_shape(&mut svg, &symbol, cx, cy, NODE_R * s, opts);

            // Put N- and lower-ring O-sulfates below the symbol and the
            // remaining O-sulfates above it.
            let mod_labels = build_modification_labels(residue);
            if vertical {
                // Bonds run through the top and bottom of a vertical symbol.
                let side = (NODE_R + 9.0) * s;
                for (text, x, anchor) in [
                    (&mod_labels.above, cx + side, "start"),
                    (&mod_labels.below, cx - side, "end"),
                ] {
                    if !text.is_empty() {
                        svg.push_str(&format!(
                            "<text x=\"{x}\" y=\"{cy}\" class=\"mod-label\" style=\"text-anchor: {anchor}\">{}</text>\n",
                            escape_xml_text(text),
                        ));
                    }
                }
            } else {
                if !mod_labels.above.is_empty() {
                    let mod_label = escape_xml_text(&mod_labels.above);
                    svg.push_str(&format!(
                        "<text x=\"{x}\" y=\"{y}\" class=\"mod-label\">{lbl}</text>\n",
                        x = cx,
                        y = cy - NODE_R * s - 14.0 * s,
                        lbl = mod_label,
                    ));
                }
                if !mod_labels.below.is_empty() {
                    let mod_label = escape_xml_text(&mod_labels.below);
                    svg.push_str(&format!(
                        "<text x=\"{x}\" y=\"{y}\" class=\"mod-label\">{lbl}</text>\n",
                        x = cx,
                        y = cy + NODE_R * s + 14.0 * s,
                        lbl = mod_label,
                    ));
                }
            }

            if opts.show_labels || is_assigned_symbol(residue, &symbol) {
                let label = escape_xml_text(&symbol.label);
                svg.push_str(&format!(
                    r#"<text x="{x}" y="{y}" class="res-label">{lbl}</text>
"#,
                    x = cx,
                    y = cy,
                    lbl = label,
                ));
            } else {
                draw_configuration_mark(&mut svg, residue, cx, cy);
            }
            if highlights.is_some() {
                svg.push_str("</g>\n");
            }
        }
    }

    svg.push_str("</svg>\n");
    Ok(svg)
}

fn render_composition_svg(
    graph: &ResidueGraph,
    opts: &RenderOptions,
    highlights: Option<&HighlightSelection>,
) -> SnfgResult<String> {
    type Group = (
        String,
        Symbol,
        ModificationLabels,
        bool,
        Option<char>,
        Vec<usize>,
    );
    let mut groups: Vec<Group> = Vec::new();
    for node in graph.inner().node_indices() {
        let residue = &graph.inner()[node];
        // Composition diagrams normally coalesce identical residues.  Keep
        // selected and unselected copies in separate groups so an explicit
        // selection remains exact even when the composition contains repeats.
        let selected =
            highlights.is_some_and(|selection| selection.node_indices.contains(&node.index()));
        let key = format!("{residue:?}|selected={selected}");
        if let Some((_, _, _, _, _, nodes)) = groups
            .iter_mut()
            .find(|(value, _, _, _, _, _)| *value == key)
        {
            nodes.push(node.index());
        } else {
            let symbol = symbol_for(residue)?;
            groups.push((
                key,
                symbol.clone(),
                build_modification_labels(residue),
                is_assigned_symbol(residue, &symbol),
                configuration_mark(residue),
                vec![node.index()],
            ));
        }
    }

    let scale = opts.scale;
    let spacing = 155.0 * scale;
    let width = (groups.len().max(1) as f64 * spacing) + 50.0 * scale;
    let height = 175.0 * scale;
    let highlight_css = if highlights.is_some() {
        MOTIF_DEEMPHASIS_CSS
    } else {
        ""
    };
    let metadata = svg_metadata(graph, opts);
    let mut svg = format!(
        r#"<svg xmlns="http://www.w3.org/2000/svg" role="img" aria-labelledby="snfg-title snfg-desc" viewBox="0 0 {width} {height}" width="{width}" height="{height}">
{metadata}<style>
  .res-label {{ font-family: {font}; font-size: 11px; fill: #000; text-anchor: middle; dominant-baseline: central; }}
{config_css}  .mod-label {{ font-family: {font}; font-size: {mod_size}px; fill: #000; text-anchor: middle; dominant-baseline: central; }}
  .count {{ font-family: {font}; font-size: {count_size}px; font-weight: 600; fill: #000; text-anchor: middle; }}
{highlight_css}</style>
"#,
        font = opts.font_family,
        mod_size = LABEL_SIZE * scale,
        config_css = configuration_label_css(graph, &opts.font_family, scale),
        count_size = 18.0 * scale,
        highlight_css = highlight_css,
        metadata = metadata,
    );
    for (index, (_, symbol, modification, assigned, mark, nodes)) in groups.iter().enumerate() {
        let x = 75.0 * scale + index as f64 * spacing;
        let y = 62.0 * scale;
        if let Some(selection) = highlights {
            let selected = nodes
                .iter()
                .any(|node| selection.node_indices.contains(node));
            let indices = nodes
                .iter()
                .map(usize::to_string)
                .collect::<Vec<_>>()
                .join(",");
            svg.push_str(&format!(
                r#"<g data-node-indices="{indices}" class="{}">
"#,
                highlight_class(selected)
            ));
        }
        draw_shape(&mut svg, symbol, x, y, NODE_R * scale, opts);
        if opts.show_labels || *assigned {
            let label = escape_xml_text(&symbol.label);
            svg.push_str(&format!(
                r#"<text x="{x}" y="{y}" class="res-label">{label}</text>
"#,
            ));
        } else if let Some(letter) = mark {
            svg.push_str(&format!(
                "<text x=\"{x}\" y=\"{y}\" class=\"config-label\">{letter}</text>\n"
            ));
        }
        if !modification.above.is_empty() {
            let modification = escape_xml_text(&modification.above);
            svg.push_str(&format!(
                r#"<text x="{x}" y="{}" class="mod-label">{modification}</text>
"#,
                y - 42.0 * scale
            ));
        }
        if !modification.below.is_empty() {
            let modification = escape_xml_text(&modification.below);
            svg.push_str(&format!(
                r#"<text x="{x}" y="{}" class="mod-label">{modification}</text>
"#,
                y + 42.0 * scale
            ));
        }
        svg.push_str(&format!(
            r#"<text x="{x}" y="{}" class="count">×{count}</text>
"#,
            y + 75.0 * scale,
            count = nodes.len()
        ));
        if highlights.is_some() {
            svg.push_str("</g>\n");
        }
    }
    svg.push_str("</svg>\n");
    Ok(svg)
}

fn is_assigned_symbol(residue: &Monosaccharide, symbol: &Symbol) -> bool {
    symbol.shape == Shape::Pentagon
        && classify_residue(residue).is_some_and(|kind| kind == ResidueKind::Assigned)
}

/// Maps layout positions onto the canvas for one orientation.
struct Frame {
    orientation: Orientation,
    width: f64,
    height: f64,
    /// Margin on the reducing-end side of the canvas.
    root_pad: f64,
    depth_step: f64,
    scale: f64,
    lane_scale: f64,
    lane_offset: f64,
}

impl Frame {
    fn at(&self, position: &LayoutInfo) -> (f64, f64) {
        self.offset(position, 0.0, 0.0)
    }

    /// The canvas point `depth` pixels further from the reducing end and
    /// `lane` pixels further down the lanes than `position`.
    fn offset(&self, position: &LayoutInfo, depth: f64, lane: f64) -> (f64, f64) {
        let along = position.x * self.depth_step * self.scale + depth;
        let across = position.y * self.lane_scale + self.lane_offset + lane;
        match self.orientation {
            Orientation::RootRight => (self.width - self.root_pad - along, across),
            Orientation::RootLeft => (self.root_pad + along, across),
            // The upper lanes of a horizontal drawing go to the right.
            Orientation::RootBottom => (self.width - across, self.height - self.root_pad - along),
            Orientation::RootTop => (self.width - across, self.root_pad + along),
        }
    }
}

/// Length of the bond drawn to the reducing-end text, from the root's centre.
const REDUCING_END_STUB: f64 = 52.0;
const CLEAVAGE_COLOUR: &str = "#C2185B";
const CLEAVAGE_LABEL_SIZE: f64 = 17.0;
/// Half the length of a cleavage stroke, clear of the symbols beside the bond.
const CLEAVAGE_REACH: f64 = 33.0;

fn draw_reducing_end(svg: &mut String, text: &str, root: &LayoutInfo, frame: &Frame, dimmed: bool) {
    let s = frame.scale;
    let (x1, y1) = frame.at(root);
    let (x2, y2) = frame.offset(root, -REDUCING_END_STUB * s, 0.0);
    if dimmed {
        svg.push_str("<g class=\"motif-dimmed\">\n");
    }
    svg.push_str(&format!(
        "<line x1=\"{x1}\" y1=\"{y1}\" x2=\"{x2}\" y2=\"{y2}\" class=\"bond\"/>\n"
    ));
    let (reach, anchor) = match frame.orientation {
        Orientation::RootRight => (8.0, "start"),
        Orientation::RootLeft => (8.0, "end"),
        Orientation::RootBottom | Orientation::RootTop => (8.0 + LABEL_SIZE / 2.0, "middle"),
    };
    let (x, y) = frame.offset(root, -(REDUCING_END_STUB + reach) * s, 0.0);
    svg.push_str(&format!(
        "<text x=\"{x}\" y=\"{y}\" class=\"reducing-end\" style=\"text-anchor: {anchor}\" data-reducing-end=\"true\">{}</text>\n",
        escape_xml_text(text)
    ));
    if dimmed {
        svg.push_str("</g>\n");
    }
}

/// A repeat unit as the layout positions its bracket is drawn around.
#[derive(Debug, Clone)]
struct RepeatBracket {
    near: f64,
    far: f64,
    top: f64,
    bottom: f64,
    count: String,
}

fn repeat_count_label(count: &RepeatCount) -> String {
    match count {
        RepeatCount::Unknown
        | RepeatCount::Range {
            min: None,
            max: None,
        } => "n".into(),
        RepeatCount::Exact(count) => format!("n = {count}"),
        RepeatCount::Range {
            min: Some(min),
            max: Some(max),
        } => format!("n = {min}\u{2013}{max}"),
        RepeatCount::Range {
            min: Some(min),
            max: None,
        } => format!("n \u{2265} {min}"),
        RepeatCount::Range {
            min: None,
            max: Some(max),
        } => format!("n \u{2264} {max}"),
    }
}

fn repeat_brackets(
    graph: &ResidueGraph,
    layout: &HashMap<usize, LayoutInfo>,
) -> Vec<RepeatBracket> {
    graph
        .inner()
        .edge_references()
        .filter_map(|edge| {
            let count = edge.weight().repeat.as_ref()?;
            let unit = graph.repeat_unit(edge.id());
            let positions = unit
                .iter()
                .filter_map(|node| layout.get(node))
                .collect::<Vec<_>>();
            let bound = |value: fn(&LayoutInfo) -> f64, start: f64, pick: fn(f64, f64) -> f64| {
                positions
                    .iter()
                    .map(|position| value(position))
                    .fold(start, pick)
            };
            (!positions.is_empty()).then(|| RepeatBracket {
                near: bound(|p| p.x, f64::MAX, f64::min),
                far: bound(|p| p.x, f64::MIN, f64::max),
                top: bound(|p| p.y, f64::MAX, f64::min),
                bottom: bound(|p| p.y, f64::MIN, f64::max),
                count: repeat_count_label(count),
            })
        })
        .collect()
}

fn path_through(points: &[(f64, f64)]) -> String {
    points
        .iter()
        .enumerate()
        .map(|(index, (x, y))| format!("{}{x} {y}", if index == 0 { "M" } else { "L" }))
        .collect::<Vec<_>>()
        .join(" ")
}

fn draw_repeat_bracket(svg: &mut String, repeat: &RepeatBracket, frame: &Frame, dimmed: bool) {
    let s = frame.scale;
    // Halfway to the neighbouring residue, clear of the largest symbols.
    let gap = (0.4 * frame.depth_step * s).max(33.0 * s);
    let reach = NODE_R * 1.2 * s + 8.0 * s;
    let tick = 12.0 * s;
    let near = LayoutInfo {
        x: repeat.near,
        y: repeat.top,
    };
    let far = LayoutInfo {
        x: repeat.far,
        y: repeat.top,
    };
    let span = (repeat.bottom - repeat.top) * frame.lane_scale;
    if dimmed {
        svg.push_str("<g class=\"motif-dimmed\">\n");
    }
    for (anchor, depth, inward) in [(&near, -gap, tick), (&far, gap, -tick)] {
        let corners = [
            frame.offset(anchor, depth + inward, -reach),
            frame.offset(anchor, depth, -reach),
            frame.offset(anchor, depth, span + reach),
            frame.offset(anchor, depth + inward, span + reach),
        ];
        svg.push_str(&format!(
            "<path d=\"{}\" class=\"bracket\" stroke=\"#000000\"/>\n",
            path_through(&corners)
        ));
    }
    let (x, y) = if frame.orientation.is_vertical() {
        frame.offset(&near, -gap - 18.0 * s, span / 2.0)
    } else {
        frame.offset(&near, -gap, span + reach + 15.0 * s)
    };
    svg.push_str(&format!(
        "<text x=\"{x}\" y=\"{y}\" class=\"repeat-count\">{}</text>\n",
        escape_xml_text(&repeat.count)
    ));
    if dimmed {
        svg.push_str("</g>\n");
    }
}

fn draw_floating_band(
    svg: &mut String,
    band: &FloatingBand,
    graph: &ResidueGraph,
    layout: &HashMap<usize, LayoutInfo>,
    frame: &Frame,
    linkage_text: bool,
    dimmed: bool,
) {
    let s = frame.scale;
    let inner = graph.inner();
    let point = |depth: f64, lane: f64| LayoutInfo { x: depth, y: lane };
    svg.push_str(if dimmed {
        "<g class=\"motif-dimmed\" data-floating-band=\"true\">\n"
    } else {
        "<g data-floating-band=\"true\">\n"
    });

    // The bracket spans the candidate parents with a tick above each of them.
    let (first, last) = (band.depths.first(), band.depths.last());
    let (from, to) = match (first, last) {
        (Some(first), Some(last)) if first != last => (*first, *last),
        _ => (band.stem - 0.3, band.stem + 0.3),
    };
    let mut bracket = path_through(&[
        frame.at(&point(from, band.lane)),
        frame.at(&point(to, band.lane)),
    ]);
    let tick = 12.0 * s;
    let mut ticks = band.depths.clone();
    if band.depths.len() < 2 {
        ticks = vec![from, to];
    }
    for depth in ticks {
        bracket.push(' ');
        bracket.push_str(&path_through(&[
            frame.at(&point(depth, band.lane)),
            frame.offset(&point(depth, band.lane), 0.0, tick),
        ]));
    }
    svg.push_str(&format!(
        "<path d=\"{bracket}\" class=\"bracket\" stroke=\"#000000\"/>\n"
    ));

    let lane_of = |item: &FloatingItem| match item {
        FloatingItem::Fragment(child) => layout.get(child).map_or(band.lane, |li| li.y),
        FloatingItem::Substituent(_, lane) => *lane,
    };
    let highest = band.items.iter().map(lane_of).fold(band.lane, f64::min);
    let (x1, y1) = frame.at(&point(band.stem, band.lane));
    let (x2, y2) = frame.at(&point(band.stem, highest));
    svg.push_str(&format!(
        "<line x1=\"{x1}\" y1=\"{y1}\" x2=\"{x2}\" y2=\"{y2}\" class=\"uncertain\"/>\n"
    ));

    for item in &band.items {
        let lane = lane_of(item);
        let stem_point = frame.at(&point(band.stem, lane));
        match item {
            FloatingItem::Fragment(child) => {
                let reach = layout.get(child).map_or(band.stem, |li| li.x);
                let (cx, cy) = frame.at(&point(reach, lane));
                svg.push_str(&format!(
                    "<line x1=\"{}\" y1=\"{}\" x2=\"{cx}\" y2=\"{cy}\" class=\"uncertain\"/>\n",
                    stem_point.0, stem_point.1
                ));
                if linkage_text {
                    let residue = &inner[NodeIndex::new(*child)];
                    let anomer = anomer_char(residue.anomeric_symbol, &residue.anomeric_prefix);
                    let (lx, ly) = frame.at(&point(band.stem + 0.5 * (reach - band.stem), lane));
                    draw_bond_label(
                        svg,
                        frame.orientation,
                        stem_point.0,
                        stem_point.1,
                        lx,
                        ly,
                        &format!("{anomer} ?"),
                        s,
                    );
                }
            }
            FloatingItem::Substituent(map_code, _) => {
                let (cx, cy) = frame.at(&point(band.stem + FLOATING_REACH, lane));
                svg.push_str(&format!(
                    "<line x1=\"{}\" y1=\"{}\" x2=\"{cx}\" y2=\"{cy}\" class=\"uncertain\"/>\n",
                    stem_point.0, stem_point.1
                ));
                let label = map_modification_label(map_code).unwrap_or("Sub");
                svg.push_str(&format!(
                    r##"<rect x="{}" y="{}" width="{}" height="{}" rx="{}" fill="#fff"/>
<text x="{cx}" y="{cy}" class="link" data-undefined-modification="true">{{{label}?}}</text>
"##,
                    cx - 42.0 * s,
                    cy - 18.0 * s,
                    84.0 * s,
                    36.0 * s,
                    6.0 * s,
                ));
            }
        }
    }
    svg.push_str("</g>\n");
}

/// A cleavage drawn across a bond: a stroke with a hook turned towards the
/// fragment that is kept, and the ion label past the end of the stroke.
fn draw_cleavage(
    svg: &mut String,
    mark: &CleavageMark,
    parent: (f64, f64),
    child: (f64, f64),
    linkage_text: bool,
    scale: f64,
) {
    let (dx, dy) = (child.0 - parent.0, child.1 - parent.1);
    let length = (dx * dx + dy * dy).sqrt();
    if length < 1.0 {
        return;
    }
    let along = (dx / length, dy / length);
    // Past the linkage text, which sits at the middle of the bond.
    let share = if linkage_text { 0.68 } else { 0.5 };
    let (cx, cy) = (parent.0 + dx * share, parent.1 + dy * share);
    // Across the bond, upwards or (on an upright bond) to the left.
    let mut across = (-along.1, along.0);
    if across.1 > 0.0 || (across.1 == 0.0 && across.0 > 0.0) {
        across = (-across.0, -across.1);
    }
    // The two sides hook at opposite ends, so a B and a Y on one bond read
    // as the usual pair.
    let (side, kept) = match mark.side {
        CleavageSide::NonReducing => ("non-reducing", 1.0),
        CleavageSide::Reducing => ("reducing", -1.0),
    };
    let across = (across.0 * kept, across.1 * kept);
    let reach = CLEAVAGE_REACH * scale;
    let hook = 12.0 * scale * kept;
    let end = (cx + across.0 * reach, cy + across.1 * reach);
    let points = [
        (cx - across.0 * reach, cy - across.1 * reach),
        end,
        (end.0 + along.0 * hook, end.1 + along.1 * hook),
    ];
    let (label_x, label_y) = (
        end.0 + across.0 * 14.0 * scale,
        end.1 + across.1 * 14.0 * scale,
    );
    let anchor = if across.0.abs() > 0.5 {
        if across.0 > 0.0 { "start" } else { "end" }
    } else {
        "middle"
    };
    svg.push_str(&format!(
        "<g data-cleavage=\"{}\" data-kept=\"{side}\">\n<path d=\"{}\" class=\"cleavage\"/>\n<text x=\"{label_x}\" y=\"{label_y}\" class=\"cleavage-label\" style=\"text-anchor: {anchor}\">{}</text>\n</g>\n",
        mark.edge,
        path_through(&points),
        escape_xml_text(&mark.label)
    ));
}

/// A bond label. Vertical drawings keep the text of an upright bond readable
/// and beside the bond instead of turning it with the bond.
#[allow(clippy::too_many_arguments)]
fn draw_bond_label(
    svg: &mut String,
    orientation: Orientation,
    px: f64,
    py: f64,
    cx: f64,
    cy: f64,
    label: &str,
    scale: f64,
) {
    if orientation.is_vertical() && (cx - px).abs() < (cy - py).abs() * 0.6 {
        svg.push_str(&format!(
            "<text x=\"{x}\" y=\"{y}\" class=\"link\" style=\"text-anchor: start\">{label}</text>\n",
            x = (px + cx) / 2.0 + 9.0 * scale,
            y = (py + cy) / 2.0,
            label = escape_xml_text(label),
        ));
    } else {
        draw_linkage_text(svg, px, py, cx, cy, label, scale);
    }
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
        if first.1 < 0.0 { first } else { second }
    } else {
        (0.0, -14.0 * scale)
    };
    let mut angle = dy.atan2(dx).to_degrees();
    if angle > 90.0 {
        angle -= 180.0;
    } else if angle < -90.0 {
        angle += 180.0;
    }
    let label = escape_xml_text(label);
    svg.push_str(&format!(
        r#"<text x="0" y="0" class="link" transform="translate({x},{y}) rotate({angle})">{label}</text>
"#,
        x = mx + ox,
        y = my + oy,
    ));
}

fn escape_xml_text(value: &str) -> String {
    value
        .replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
        .replace('\'', "&apos;")
}

fn retained_source_notation(graph: &ResidueGraph) -> Option<SourceNotation> {
    graph
        .source_wurcs()
        .map(|value| SourceNotation::new("wurcs", value))
        .or_else(|| {
            graph
                .source_iupac()
                .map(|value| SourceNotation::new("iupac-condensed", value))
        })
        .or_else(|| {
            graph
                .source_iupac_extended()
                .map(|value| SourceNotation::new("iupac-extended", value))
        })
        .or_else(|| {
            graph
                .source_glycam()
                .map(|value| SourceNotation::new("glycam", value))
        })
}

fn svg_metadata(graph: &ResidueGraph, opts: &RenderOptions) -> String {
    let iupac = crabwurcs_iupac::write_iupac_condensed_canonical(graph).ok();
    let wurcs = crabwurcs_core::write_wurcs_canonical(graph).ok();
    let source = opts
        .source_notation
        .clone()
        .or_else(|| retained_source_notation(graph));

    let title = iupac
        .as_deref()
        .filter(|value| !value.is_empty())
        .map(|value| format!("SNFG glycan: {value}"))
        .unwrap_or_else(|| "SNFG glycan".into());
    let desc = match (iupac.as_deref(), wurcs.as_deref()) {
        (Some(iupac), Some(wurcs)) => {
            format!("Canonical IUPAC condensed: {iupac}. Canonical WURCS: {wurcs}.")
        }
        (Some(iupac), None) => {
            format!("Canonical IUPAC condensed: {iupac}. Canonical WURCS is unavailable.")
        }
        (None, Some(wurcs)) => {
            format!("Canonical IUPAC condensed is unavailable. Canonical WURCS: {wurcs}.")
        }
        (None, None) => "Canonical IUPAC condensed and canonical WURCS are unavailable.".into(),
    };

    let mut metadata = format!(
        "<title id=\"snfg-title\">{}</title>\n\
<desc id=\"snfg-desc\">{}</desc>\n\
<metadata id=\"crabwurcs-notations\">\n\
  <crabwurcs:notations xmlns:crabwurcs=\"https://github.com/Ojas-Singh/crabWURCS/ns/metadata/1\">\n",
        escape_xml_text(&title),
        escape_xml_text(&desc),
    );
    match iupac {
        Some(value) => metadata.push_str(&format!(
            "    <crabwurcs:iupac-condensed canonical=\"true\" available=\"true\">{}</crabwurcs:iupac-condensed>\n",
            escape_xml_text(&value)
        )),
        None => metadata.push_str(
            "    <crabwurcs:iupac-condensed canonical=\"true\" available=\"false\"/>\n",
        ),
    }
    match wurcs {
        Some(value) => metadata.push_str(&format!(
            "    <crabwurcs:wurcs canonical=\"true\" available=\"true\">{}</crabwurcs:wurcs>\n",
            escape_xml_text(&value)
        )),
        None => {
            metadata.push_str("    <crabwurcs:wurcs canonical=\"true\" available=\"false\"/>\n")
        }
    }
    if let Some(source) = source {
        metadata.push_str(&format!(
            "    <crabwurcs:source format=\"{}\">{}</crabwurcs:source>\n",
            escape_xml_text(source.format.trim()),
            escape_xml_text(source.value.trim())
        ));
    }
    metadata.push_str("  </crabwurcs:notations>\n</metadata>\n");
    metadata
}

fn empty_svg(graph: &ResidueGraph, opts: &RenderOptions) -> String {
    format!(
        r##"<svg xmlns="http://www.w3.org/2000/svg" role="img" aria-labelledby="snfg-title snfg-desc" viewBox="0 0 120 40" width="120" height="40">
{}  <text x="10" y="25" font-family="sans-serif" font-size="11" fill="#999">(empty)</text>
</svg>
"##,
        svg_metadata(graph, opts)
    )
}

// ── Modification label ─────────────────────────────────────────────────────

#[derive(Clone, Debug, Default, PartialEq, Eq)]
struct ModificationLabels {
    above: String,
    below: String,
}

/// Build short O-sulfation labels around an SNFG symbol.
///
/// N-sulfation is written as `NS` below the diagonally divided NSquare.
/// O-sulfates at positions 2 and 3 are also placed below the symbol; positions
/// 4 and above are placed above it. For example, GlcNS6S yields `6S` above and
/// `NS` below, while IdoA2S yields `2S` below.
fn build_modification_labels(res: &Monosaccharide) -> ModificationLabels {
    let mut above_positions: Vec<u8> = Vec::new();
    let mut below_positions: Vec<u8> = Vec::new();
    let mut has_n_sulfate = false;

    for m in &res.modifications {
        let desc = &m.descriptor;
        if desc.contains("NSO") {
            has_n_sulfate = true;
        }
        if desc.contains("OSO") {
            if m.position.0 <= 3 {
                below_positions.push(m.position.0);
            } else {
                above_positions.push(m.position.0);
            }
        }
    }

    above_positions.sort();
    above_positions.dedup();
    below_positions.sort();
    below_positions.dedup();

    let format_positions = |positions: Vec<u8>| {
        positions
            .into_iter()
            .map(|position| format!("{position}S"))
            .collect::<String>()
    };
    ModificationLabels {
        above: format_positions(above_positions),
        below: format!(
            "{}{}",
            if has_n_sulfate { "NS" } else { "" },
            format_positions(below_positions)
        ),
    }
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
        Shape::DividedTriangle => {
            let h = r * 1.732;
            let top = cy - h * 0.667;
            let bottom = cy + h * 0.333;
            svg.push_str(&format!(
                r#"<polygon points="{cx},{top} {left},{bottom} {right},{bottom}" fill="white" stroke="{stroke}" stroke-width="{sw}"/>
"#,
                left = cx - r,
                right = cx + r,
            ));
            svg.push_str(&format!(
                r#"<polygon points="{cx},{top} {cx},{bottom} {right},{bottom}" fill="{fill}" stroke="none"/>
<line x1="{cx}" y1="{top}" x2="{cx}" y2="{bottom}" stroke="{stroke}" stroke-width="1.5"/>
"#,
                right = cx + r,
            ));
        }
        Shape::FlatRectangle => {
            svg.push_str(&format!(
                r#"<rect x="{x}" y="{y}" width="{w}" height="{h}" fill="{fill}" stroke="{stroke}" stroke-width="{sw}"/>
"#,
                x = cx - r,
                y = cy - r * 0.38,
                w = r * 2.0,
                h = r * 0.76,
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
        Shape::FlatDiamond => {
            let dx = r * 1.25;
            let dy = r * 0.62;
            svg.push_str(&format!(
                r#"<polygon points="{cx},{top} {right},{cy} {cx},{bottom} {left},{cy}" fill="{fill}" stroke="{stroke}" stroke-width="{sw}"/>
"#,
                top = cy - dy,
                right = cx + dx,
                bottom = cy + dy,
                left = cx - dx,
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
        Shape::FlatHexagon => {
            let dx = r * 1.15;
            let shoulder = r * 0.72;
            let dy = r * 0.58;
            svg.push_str(&format!(
                r#"<polygon points="{l0},{cy} {l1},{top} {r1},{top} {r0},{cy} {r1},{bottom} {l1},{bottom}" fill="{fill}" stroke="{stroke}" stroke-width="{sw}"/>
"#,
                l0 = cx - dx,
                l1 = cx - shoulder,
                r1 = cx + shoulder,
                r0 = cx + dx,
                top = cy - dy,
                bottom = cy + dy,
            ));
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
        let corpus_lines = include_str!("../data/glycoshape_notations.tsv")
            .lines()
            .chain(include_str!("../data/glycoshape_derived_notations.tsv").lines());
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
    fn n_and_o_sulfated_glucosamine_uses_nsquare_and_positioned_o_sulfates() {
        let g = parse(
            "WURCS=2.0/2,2,1/[u2122h_2*NSO/3=O/3=O_6*OSO/3=O/3=O][a2121A-1a_1-5_2*OSO/3=O/3=O]/1-2/a4-b1",
        );
        let root = g.residue(g.root().unwrap()).unwrap();
        let symbol = symbol_for(root).unwrap();
        assert_eq!(symbol.shape, Shape::NSquare);

        let svg = render_svg(&g).unwrap();
        assert!(svg.contains("SNFG glycan: IdoA2S(a1-4)GlcNS6S"));
        assert!(svg.contains(">6S</text>"));
        assert!(svg.contains(">2S</text>"));
        assert!(svg.contains(">NS</text>"));
        assert!(!svg.contains(">S6S</text>"));
        assert!(svg.contains("class=\"mod-label\""));
        assert!(
            svg.contains(
                ".mod-label { font-family: Arial, Helvetica, sans-serif; font-size: 20px;"
            )
        );
    }

    #[test]
    fn render_symbol_svg_draws_one_shape_on_a_transparent_canvas() {
        let glc = crabwurcs_core::residue_from_kind(ResidueKind::Glc).unwrap();
        let svg = render_symbol_svg(&glc, &RenderOptions::default()).unwrap();
        assert!(svg.starts_with("<svg xmlns=\"http://www.w3.org/2000/svg\""));
        assert!(svg.contains("<circle"));
        assert!(svg.contains("fill=\"#0072BC\""));
        // Single-symbol output must not leak the full-graph scaffolding.
        assert!(!svg.contains("class=\"bond\""));
        assert!(!svg.contains("<metadata"));
        assert!(!svg.contains("<title"));

        // Labels are opt-in.
        let labeled = render_symbol_svg(
            &glc,
            &RenderOptions {
                show_labels: true,
                ..RenderOptions::default()
            },
        )
        .unwrap();
        assert!(labeled.contains(">Glc</text>"));
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
            ("1121", "L-Gul", colour::GUL),
            ("2222", "All", colour::ALL),
            ("2221", "L-Tal", colour::TAL),
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
    fn enantiomers_use_the_snfg_symbol_of_their_name_with_a_configuration_letter() {
        // GlycoShape GS00025 (D-Fuc) and GS00657 (D-Rha) were drawn with the
        // Gal circle and the Man circle.
        for (wurcs, shape, fill, letter) in [
            (
                "WURCS=2.0/1,1,0/[a2112m-1a_1-5]/1/",
                Shape::Triangle,
                colour::RED,
                'D',
            ),
            (
                "WURCS=2.0/1,1,0/[a1122m-1a_1-5]/1/",
                Shape::Triangle,
                colour::GREEN,
                'D',
            ),
            (
                "WURCS=2.0/1,1,0/[a122h-1b_1-4]/1/",
                Shape::Star,
                colour::GREEN,
                'D',
            ),
            (
                "WURCS=2.0/1,1,0/[a2211h-1a_1-5]/1/",
                Shape::Circle,
                colour::GREEN,
                'L',
            ),
            (
                "WURCS=2.0/1,1,0/[a1121h-1a_1-5]/1/",
                Shape::Circle,
                colour::ORANGE,
                'L',
            ),
            (
                "WURCS=2.0/1,1,0/[a2221h-1a_1-5]/1/",
                Shape::Circle,
                colour::LIGHT_BLUE,
                'L',
            ),
            (
                "WURCS=2.0/1,1,0/[a1211m-1a_1-5]/1/",
                Shape::Triangle,
                colour::BLUE,
                'L',
            ),
            (
                "WURCS=2.0/1,1,0/[a2112m-1a_1-5_2*NCC/3=O]/1/",
                Shape::DividedTriangle,
                colour::RED,
                'D',
            ),
        ] {
            let g = parse(wurcs);
            let residue = g.residue(g.root().unwrap()).unwrap();
            let symbol = symbol_for(residue).unwrap();
            assert_eq!((symbol.shape, symbol.fill), (shape, fill), "{wurcs}");
            assert_eq!(configuration_mark(residue), Some(letter), "{wurcs}");
            let svg = render_svg(&g).unwrap();
            assert!(
                svg.contains(&format!("class=\"config-label\">{letter}</text>")),
                "{wurcs}"
            );
        }
        // Default configurations carry no letter.
        let fuc = parse("WURCS=2.0/1,1,0/[a1221m-1a_1-5]/1/");
        assert!(!render_svg(&fuc).unwrap().contains("config-label\">"));

        // Every concrete registry residue's enantiomer keeps its symbol.
        for kind in ResidueKind::ALL {
            if kind.is_generic() || kind.unique_residue().is_none() {
                continue;
            }
            let mut residue = crabwurcs_core::residue_from_kind(*kind).unwrap();
            let expected = registered_symbol(*kind, None);
            residue.skeleton_code = crabwurcs_core::mirror_skeleton_code(&residue.skeleton_code);
            residue.residue_kind = None;
            if crabwurcs_core::classify_residue(&residue).is_some() {
                // The mirror image is itself a registry entry (Abe/Col).
                continue;
            }
            let symbol = symbol_for(&residue).unwrap();
            assert_eq!(
                (symbol.shape, symbol.fill),
                (expected.shape, expected.fill),
                "{kind:?}"
            );
        }
    }

    #[test]
    fn test_symbol_xyl_uses_exact_snfg_orange() {
        let g = parse("WURCS=2.0/1,1,0/[u212h]/1/");
        let res = g.residue(g.root().unwrap()).unwrap();
        let sym = symbol_for(res).unwrap();
        assert_eq!(sym.shape, Shape::Star);
        assert_eq!(sym.fill, "#F47920");
    }

    #[test]
    fn lyx_absolute_configurations_use_yellow() {
        for (code, label) in [("112", "Lyx"), ("221", "L-Lyx")] {
            let g = parse(&format!("WURCS=2.0/1,1,0/[u{code}h]/1/"));
            let symbol = symbol_for(g.residue(g.root().unwrap()).unwrap()).unwrap();
            assert_eq!(symbol.shape, Shape::Star, "{code}");
            assert_eq!(symbol.label, label, "{code}");
            assert_eq!(symbol.fill, colour::YELLOW, "{code}");
        }
    }

    #[test]
    fn official_registry_covers_every_snfg_shape_and_colour() {
        for kind in ResidueKind::ALL {
            let symbol = registered_symbol(*kind, Some("Example"));
            assert!(!symbol.label.is_empty(), "{kind:?}");
        }

        for (kind, shape, fill) in [
            (ResidueKind::Hex, Shape::Circle, colour::WHITE),
            (ResidueKind::Glc, Shape::Circle, colour::BLUE),
            (ResidueKind::ManNAc, Shape::Square, colour::GREEN),
            (ResidueKind::GalN, Shape::NSquare, colour::YELLOW),
            (ResidueKind::IdoA, Shape::SplitDiamondBottom, colour::BROWN),
            (ResidueKind::Fuc, Shape::Triangle, colour::RED),
            (ResidueKind::QuiNAc, Shape::DividedTriangle, colour::BLUE),
            (ResidueKind::Dig, Shape::FlatRectangle, colour::PURPLE),
            (ResidueKind::Xyl, Shape::Star, colour::ORANGE),
            (ResidueKind::Neu5Gc, Shape::Diamond, colour::LIGHT_BLUE),
            (ResidueKind::Leg, Shape::FlatDiamond, colour::YELLOW),
            (ResidueKind::MurNAc, Shape::FlatHexagon, colour::PURPLE),
            (ResidueKind::Psi, Shape::Pentagon, colour::PINK),
        ] {
            let symbol = registered_symbol(kind, None);
            assert_eq!(symbol.shape, shape, "{kind:?}");
            assert_eq!(symbol.fill, fill, "{kind:?}");
        }
    }

    #[test]
    fn arbitrary_name_uses_assigned_white_pentagon_and_safe_label() {
        let mut residue = crabwurcs_core::residue_from_kind(ResidueKind::Hex).unwrap();
        residue.residue_kind = None;
        residue.display_name = Some("<foo&bar>".into());
        let symbol = symbol_for(&residue).unwrap();
        assert_eq!(symbol.shape, Shape::Pentagon);
        assert_eq!(symbol.fill, colour::WHITE);
        assert_eq!(symbol.label, "F");

        let mut graph = ResidueGraph::new();
        graph.add_residue(residue);
        let svg = render_svg(&graph).unwrap();
        assert!(svg.contains(">F</text>"));
        assert!(!svg.contains("<foo"));
    }

    #[test]
    fn every_registered_symbol_can_be_rendered() {
        for &kind in ResidueKind::ALL {
            let mut residue = crabwurcs_core::residue_from_kind(kind)
                .unwrap_or_else(|_| crabwurcs_core::residue_from_kind(ResidueKind::Hex).unwrap());
            residue.residue_kind = Some(kind);

            let mut graph = ResidueGraph::new();
            graph.add_residue(residue);
            let svg = render_svg(&graph).unwrap();
            assert!(
                svg.contains("<svg") && svg.contains("</svg>"),
                "failed to render {}",
                kind.canonical_name()
            );
        }
    }

    #[test]
    fn dynamic_svg_text_is_xml_escaped() {
        assert_eq!(escape_xml_text("<&>\"'"), "&lt;&amp;&gt;&quot;&apos;");
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
    fn official_sialic_acid_palette_and_fruf_symbol_are_exact() {
        let neu5gc = parse("WURCS=2.0/1,1,0/[AUd21122h_5*NCCO/3=O]/1/");
        let neu5gc = symbol_for(neu5gc.residue(neu5gc.root().unwrap()).unwrap()).unwrap();
        assert_eq!(neu5gc.shape, Shape::Diamond);
        assert_eq!(neu5gc.fill, colour::LIGHT_BLUE);

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
        assert_eq!(kdo_symbol.shape, Shape::FlatHexagon);
        assert_eq!(kdo_symbol.fill, colour::KDO);
        assert_eq!(kdo_symbol.label, "Kdo");

        let bac = parse("WURCS=2.0/1,1,0/[u2122m_2*NCC/3=O_4*NCC/3=O]/1/");
        let bac_symbol = symbol_for(bac.residue(bac.root().unwrap()).unwrap()).unwrap();
        assert_eq!(bac_symbol.shape, Shape::FlatHexagon);
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
        assert!(svg.contains("aria-labelledby=\"snfg-title snfg-desc\""));
        assert_eq!(svg.matches("class=\"count\"").count(), 4);
        assert!(svg.contains("×6"));
        assert!(svg.contains("×5"));
    }

    #[test]
    fn svg_metadata_contains_canonical_and_original_notations() {
        let graph = crabwurcs_iupac::parse_iupac_condensed("Glucose").unwrap();
        let svg = render_svg(&graph).unwrap();
        assert!(svg.contains("<title id=\"snfg-title\">SNFG glycan: Glc</title>"));
        assert!(svg.contains("<metadata id=\"crabwurcs-notations\">"));
        assert!(svg.contains(
            "<crabwurcs:iupac-condensed canonical=\"true\" available=\"true\">Glc</crabwurcs:iupac-condensed>"
        ));
        assert!(svg.contains("<crabwurcs:wurcs canonical=\"true\" available=\"true\">"));
        assert!(
            svg.contains("<crabwurcs:source format=\"iupac-condensed\">Glucose</crabwurcs:source>")
        );
        assert!(svg.contains("aria-labelledby=\"snfg-title snfg-desc\""));
    }

    #[test]
    fn metadata_is_best_effort_and_escapes_source_text() {
        let graph = crabwurcs_iupac::parse_iupac_condensed("Foo").unwrap();
        let svg = render_svg_with_options(
            &graph,
            &RenderOptions {
                source_notation: Some(SourceNotation::new("custom<&", "Foo<&")),
                ..RenderOptions::default()
            },
        )
        .unwrap();
        assert!(svg.contains(
            "<crabwurcs:iupac-condensed canonical=\"true\" available=\"true\">Foo</crabwurcs:iupac-condensed>"
        ));
        assert!(svg.contains("<crabwurcs:wurcs canonical=\"true\" available=\"false\"/>"));
        assert!(svg.contains(
            "<crabwurcs:source format=\"custom&lt;&amp;\">Foo&lt;&amp;</crabwurcs:source>"
        ));
    }

    #[test]
    fn png_is_transparent_rgba_at_twice_the_svg_dimensions() {
        let graph = crabwurcs_iupac::parse_iupac_condensed("Glc").unwrap();
        let png_bytes = render_png(&graph).unwrap();
        assert_eq!(&png_bytes[..8], b"\x89PNG\r\n\x1a\n");

        let decoder = png::Decoder::new(std::io::Cursor::new(&png_bytes));
        let mut reader = decoder.read_info().unwrap();
        let mut pixels = vec![0; reader.output_buffer_size().unwrap()];
        let info = reader.next_frame(&mut pixels).unwrap();
        assert_eq!(info.width, 280);
        assert_eq!(info.height, 280);
        assert_eq!(info.color_type, png::ColorType::Rgba);
        assert_eq!(
            pixels[3], 0,
            "top-left background pixel must be transparent"
        );
    }

    #[test]
    fn highlighted_png_uses_the_same_motif_rendering_path() {
        let graph = crabwurcs_iupac::parse_iupac_condensed("Fuc(a1-3)GlcNAc(b1-4)GlcNAc").unwrap();
        let motif = crabwurcs_iupac::parse_iupac_condensed("Fuc(a1-?)GlcNAc").unwrap();
        let png = render_png_with_motifs(&graph, &[motif], &RenderOptions::default()).unwrap();
        assert_eq!(&png[..8], b"\x89PNG\r\n\x1a\n");
    }

    #[test]
    fn explicit_selection_highlights_only_the_requested_graph_elements() {
        let graph =
            crabwurcs_iupac::parse_iupac_condensed("Fuc(a1-3)GlcNAc(b1-4)Fuc(a1-3)GlcNAc").unwrap();
        let motif = crabwurcs_iupac::parse_iupac_condensed("Fuc(a1-?)GlcNAc").unwrap();
        let all = render_svg_with_motifs(&graph, &[motif], &RenderOptions::default()).unwrap();
        assert_eq!(all.matches("class=\"motif-match\"").count(), 6);

        let mut selection = HighlightSelection::default();
        selection.node_indices.extend([0, 1]);
        selection.edge_indices.insert(0);
        let exact =
            render_svg_with_selection(&graph, &selection, &RenderOptions::default()).unwrap();
        assert_eq!(exact.matches("class=\"motif-match\"").count(), 3);
        assert_eq!(exact.matches("class=\"motif-dimmed\"").count(), 4);
        assert!(exact.contains("data-node-index=\"0\" class=\"motif-match\""));
        assert!(exact.contains("data-node-index=\"1\" class=\"motif-match\""));
        assert!(exact.contains("data-edge-index=\"0\" class=\"motif-match\""));
        let png = render_png_with_selection(&graph, &selection, &RenderOptions::default()).unwrap();
        assert_eq!(&png[..8], b"\x89PNG\r\n\x1a\n");
    }

    const N_GLYCAN: &str = "Gal(b1-4)GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)[Fuc(a1-6)]GlcNAc";
    const REPEAT: &str = "WURCS=2.0/4,5,5/[a2112h-1b_1-5_2*NCC/3=O][a2122A-1b_1-5][a2112h-1a_1-5][a2122h-1b_1-5]/1-2-1-3-4/a3-b1_b4-c1_c4-d1_d3-e1_a1-e4~2-5";
    const FRAGMENT: &str = "WURCS=2.0/6,11,10/[a2122h-1x_1-5_2*NCC/3=O][a2122h-1b_1-5_2*NCC/3=O][a1122h-1b_1-5][a1122h-1a_1-5][a2112h-1b_1-5][a1221m-1a_1-5]/1-2-3-4-2-5-4-2-6-2-5/a4-b1_a6-i1_b4-c1_c3-d1_c6-g1_d2-e1_e4-f1_g2-h1_j4-k1_j1-d4|d6|g4|g6}";
    const SUBSTITUENT: &str = "WURCS=2.0/2,2,1/[u2122h][u2112h]/1-2/a?|b?}*OCC/3=O";

    fn iupac(value: &str) -> ResidueGraph {
        crabwurcs_iupac::parse_iupac_condensed(value).unwrap()
    }

    fn with_orientation(orientation: Orientation) -> RenderOptions {
        RenderOptions {
            orientation,
            ..RenderOptions::default()
        }
    }

    /// Every number of the attributes in `names`, across the drawing.
    fn attributes(svg: &str, names: &[&str]) -> Vec<f64> {
        let mut values = Vec::new();
        for name in names {
            let pattern = format!(" {name}=\"");
            let mut rest = svg;
            while let Some(start) = rest.find(&pattern) {
                rest = &rest[start + pattern.len()..];
                let end = rest.find('"').unwrap();
                values.push(rest[..end].parse::<f64>().unwrap());
            }
        }
        values
    }

    fn canvas_size(svg: &str) -> (f64, f64) {
        let header = &svg[..svg.find('>').unwrap()];
        (
            attributes(header, &["width"])[0],
            attributes(header, &["height"])[0],
        )
    }

    /// Whether everything is placed on the canvas. Path corners count too.
    fn assert_on_canvas(svg: &str, context: &str) {
        let (width, height) = canvas_size(svg);
        let body = &svg[svg.find("</style>").unwrap()..];
        for x in attributes(body, &["cx", "x1", "x2"]) {
            assert!((0.0..=width).contains(&x), "x {x} of {width}: {context}");
        }
        for y in attributes(body, &["cy", "y1", "y2"]) {
            assert!((0.0..=height).contains(&y), "y {y} of {height}: {context}");
        }
        let mut rest = body;
        while let Some(start) = rest.find(" d=\"") {
            rest = &rest[start + 4..];
            let end = rest.find('"').unwrap();
            let numbers = rest[..end]
                .split([' ', 'M', 'L'])
                .filter(|part| !part.is_empty())
                .map(|part| part.parse::<f64>().unwrap())
                .collect::<Vec<_>>();
            for point in numbers.chunks(2) {
                assert!((0.0..=width).contains(&point[0]), "path x: {context}");
                assert!((0.0..=height).contains(&point[1]), "path y: {context}");
            }
        }
    }

    /// `(x, y)` of the first `<circle>` and the first `<rect>`.
    fn circle_and_square(svg: &str) -> ((f64, f64), (f64, f64)) {
        let circle = &svg[svg.find("<circle").unwrap()..];
        let circle = &circle[..circle.find("/>").unwrap()];
        let square = &svg[svg.find("<rect").unwrap()..];
        let square = &square[..square.find("/>").unwrap()];
        let size = attributes(square, &["width"])[0] / 2.0;
        (
            (
                attributes(circle, &["cx"])[0],
                attributes(circle, &["cy"])[0],
            ),
            (
                attributes(square, &["x"])[0] + size,
                attributes(square, &["y"])[0] + size,
            ),
        )
    }

    #[test]
    fn orientations_put_the_reducing_end_where_they_say() {
        // The GlcNAc square is the reducing end, the Gal circle the child.
        let graph = iupac("Gal(b1-4)GlcNAc");
        for (orientation, dx, dy) in [
            (Orientation::RootRight, -100.0, 0.0),
            (Orientation::RootLeft, 100.0, 0.0),
            (Orientation::RootBottom, 0.0, -100.0),
            (Orientation::RootTop, 0.0, 100.0),
        ] {
            let svg = render_svg_with_options(&graph, &with_orientation(orientation)).unwrap();
            let (gal, glcnac) = circle_and_square(&svg);
            assert_eq!(gal.0 - glcnac.0, dx, "{orientation:?}");
            assert_eq!(gal.1 - glcnac.1, dy, "{orientation:?}");
        }
    }

    #[test]
    fn the_right_to_left_orientation_is_the_default_drawing() {
        let graph = iupac(N_GLYCAN);
        assert_eq!(
            render_svg(&graph).unwrap(),
            render_svg_with_options(&graph, &with_orientation(Orientation::RootRight)).unwrap()
        );
    }

    #[test]
    fn the_higher_acceptor_arm_stays_above_or_to_the_right() {
        // Lanes come from the layout, whatever the orientation. Compare the
        // canvas positions of the two Man arms of the core.
        let graph = iupac("Man(a1-3)[Man(a1-6)]Man");
        let arms = |orientation| {
            let svg = render_svg_with_options(&graph, &with_orientation(orientation)).unwrap();
            let body = &svg[svg.find("</style>").unwrap()..];
            let xs = attributes(body, &["cx"]);
            let ys = attributes(body, &["cy"]);
            // Node order: core (root), then the arms as parsed: α6, α3.
            ((xs[1], ys[1]), (xs[2], ys[2]))
        };
        let (alpha6, alpha3) = arms(Orientation::RootRight);
        assert!(alpha6.1 < alpha3.1, "α6 above α3 in the default drawing");
        let (alpha6, alpha3) = arms(Orientation::RootLeft);
        assert!(alpha6.1 < alpha3.1, "α6 above α3 when mirrored");
        for orientation in [Orientation::RootBottom, Orientation::RootTop] {
            let (alpha6, alpha3) = arms(orientation);
            assert!(alpha6.0 > alpha3.0, "α6 right of α3: {orientation:?}");
        }
    }

    #[test]
    fn orientations_parse_their_names_and_the_common_spellings() {
        for orientation in [
            Orientation::RootRight,
            Orientation::RootLeft,
            Orientation::RootBottom,
            Orientation::RootTop,
        ] {
            assert_eq!(orientation.as_str().parse(), Ok(orientation));
        }
        assert_eq!("top-down".parse(), Ok(Orientation::RootTop));
        assert_eq!("left_to_right".parse(), Ok(Orientation::RootLeft));
        assert_eq!("VERTICAL".parse(), Ok(Orientation::RootBottom));
        assert!("diagonal".parse::<Orientation>().is_err());
    }

    #[test]
    fn every_corpus_glycan_stays_on_the_canvas_in_every_layout() {
        let corpus = include_str!("../data/glycoshape_notations.tsv")
            .lines()
            .chain(include_str!("../data/glycoshape_derived_notations.tsv").lines());
        for line in corpus {
            let wurcs = line.split('\t').next().unwrap();
            let graph = parse_wurcs(wurcs).unwrap();
            for orientation in [
                Orientation::RootRight,
                Orientation::RootLeft,
                Orientation::RootBottom,
                Orientation::RootTop,
            ] {
                for compact in [false, true] {
                    let opts = RenderOptions {
                        orientation,
                        compact,
                        reducing_end: Some("Ser/Thr".into()),
                        repeat_brackets: true,
                        floating_fragments: true,
                        ..RenderOptions::default()
                    };
                    let svg = render_svg_with_options(&graph, &opts).unwrap();
                    assert_on_canvas(&svg, &format!("{wurcs} {orientation:?} {compact}"));
                }
            }
        }
    }

    #[test]
    fn compact_drawings_are_smaller_and_carry_no_linkage_text() {
        let graph = iupac(N_GLYCAN);
        let full = render_svg(&graph).unwrap();
        let opts = RenderOptions {
            compact: true,
            ..RenderOptions::default()
        };
        let compact = render_svg_with_options(&graph, &opts).unwrap();
        let (full_w, full_h) = canvas_size(&full);
        let (compact_w, compact_h) = canvas_size(&compact);
        assert!(compact_w < full_w && compact_h < full_h);
        assert!(full.contains("class=\"link\""));
        assert!(!compact.contains("class=\"link\""));
        assert_on_canvas(&compact, "compact N-glycan");

        // Labels above and below a symbol keep their room.
        let sulfated = parse(
            "WURCS=2.0/2,2,1/[u2122h_2*NSO/3=O/3=O_6*OSO/3=O/3=O][a2121A-1a_1-5_2*OSO/3=O/3=O]/1-2/a4-b1",
        );
        let svg = render_svg_with_options(&sulfated, &opts).unwrap();
        let body = &svg[svg.find("</style>").unwrap()..];
        let lowest = attributes(body, &["y"])
            .into_iter()
            .fold(f64::MIN, f64::max);
        assert!(lowest + LABEL_SIZE / 2.0 <= canvas_size(&svg).1);
    }

    #[test]
    fn the_reducing_end_label_widens_the_canvas_and_is_escaped() {
        let graph = iupac("Gal(b1-4)GlcNAc");
        let plain = canvas_size(&render_svg(&graph).unwrap());
        for (orientation, grows_width) in [
            (Orientation::RootRight, true),
            (Orientation::RootLeft, true),
            (Orientation::RootBottom, false),
            (Orientation::RootTop, false),
        ] {
            let opts = RenderOptions {
                orientation,
                reducing_end: Some("Ser<Thr>".into()),
                ..RenderOptions::default()
            };
            let svg = render_svg_with_options(&graph, &opts).unwrap();
            assert!(svg.contains(">Ser&lt;Thr&gt;</text>"), "{orientation:?}");
            assert!(svg.contains("data-reducing-end=\"true\""));
            assert_on_canvas(&svg, "reducing end");
            let (width, height) = canvas_size(&svg);
            if grows_width {
                assert!(width > plain.0, "{orientation:?}");
            } else {
                // A vertical drawing is as tall as a horizontal one is wide.
                assert!(height > plain.0, "{orientation:?}");
            }
        }

        // A blank label draws nothing.
        let blank = RenderOptions {
            reducing_end: Some("  ".into()),
            ..RenderOptions::default()
        };
        assert_eq!(
            render_svg_with_options(&graph, &blank).unwrap(),
            render_svg(&graph).unwrap()
        );
    }

    #[test]
    fn the_reducing_end_label_dims_with_its_residue() {
        let graph = iupac("Gal(b1-4)GlcNAc");
        let opts = RenderOptions {
            reducing_end: Some("Asn".into()),
            ..RenderOptions::default()
        };
        // GlcNAc, the reducing end, is node 0 and Gal node 1.
        let group = "<g class=\"motif-dimmed\">
<line";
        let outside = HighlightSelection::from_residues(&graph, [1]).unwrap();
        let inside = HighlightSelection::from_residues(&graph, [0]).unwrap();
        let dimmed = render_svg_with_selection(&graph, &outside, &opts).unwrap();
        let kept = render_svg_with_selection(&graph, &inside, &opts).unwrap();
        assert!(dimmed.contains(group), "{dimmed}");
        assert!(!kept.contains(group), "{kept}");
    }

    #[test]
    fn repeat_brackets_replace_the_closing_bond_and_state_the_count() {
        let graph = parse(REPEAT);
        let plain = render_svg(&graph).unwrap();
        assert!(!plain.contains("class=\"bracket\""));
        assert_eq!(plain.matches("class=\"uncertain\"").count(), 1);

        let opts = RenderOptions {
            repeat_brackets: true,
            ..RenderOptions::default()
        };
        let svg = render_svg_with_options(&graph, &opts).unwrap();
        assert_eq!(svg.matches("class=\"bracket\"").count(), 2);
        assert_eq!(svg.matches("class=\"uncertain\"").count(), 0);
        assert!(svg.contains(">n = 2\u{2013}5</text>"), "{svg}");
        assert_eq!(svg.matches("class=\"bond\"").count(), 4);
        for orientation in [
            Orientation::RootLeft,
            Orientation::RootBottom,
            Orientation::RootTop,
        ] {
            let opts = RenderOptions {
                orientation,
                ..opts.clone()
            };
            assert_on_canvas(
                &render_svg_with_options(&graph, &opts).unwrap(),
                "repeat bracket",
            );
        }

        for (count, label) in [
            (RepeatCount::Unknown, "n"),
            (RepeatCount::Exact(3), "n = 3"),
            (
                RepeatCount::Range {
                    min: Some(2),
                    max: Some(5),
                },
                "n = 2\u{2013}5",
            ),
            (
                RepeatCount::Range {
                    min: Some(2),
                    max: None,
                },
                "n \u{2265} 2",
            ),
            (
                RepeatCount::Range {
                    min: None,
                    max: Some(5),
                },
                "n \u{2264} 5",
            ),
            (
                RepeatCount::Range {
                    min: None,
                    max: None,
                },
                "n",
            ),
        ] {
            assert_eq!(repeat_count_label(&count), label);
        }
    }

    #[test]
    fn a_repeat_bracket_leaves_a_cap_outside() {
        // Gal(b1-3)[GlcNAc(b1-4)Gal(b1-3)]nGalNAc: the Gal cap is outside.
        let mut graph = iupac("Gal(b1-3)GlcNAc(b1-4)Gal(b1-3)GalNAc");
        let mut closure = crabwurcs_core::Linkage::new(
            crabwurcs_core::CarbonPosition(3),
            crabwurcs_core::CarbonPosition(1),
        );
        closure.repeat = Some(RepeatCount::Exact(2));
        // Nodes from the reducing end: GalNAc 0, Gal 1, GlcNAc 2, Gal 3.
        graph.add_linkage(NodeIndex::new(2), NodeIndex::new(1), closure);
        let edge = graph
            .inner()
            .edge_references()
            .find(|edge| edge.weight().repeat.is_some())
            .unwrap()
            .id();
        assert_eq!(graph.repeat_unit(edge), BTreeSet::from([1, 2]));

        let layout = compute_layout(&graph, graph.root().unwrap());
        let [bracket] = repeat_brackets(&graph, &layout).try_into().unwrap();
        assert_eq!((bracket.near, bracket.far), (1.0, 2.0));
    }

    #[test]
    fn undefined_fragments_float_over_a_bracket_when_asked() {
        let graph = parse(FRAGMENT);
        let opts = RenderOptions {
            floating_fragments: true,
            show_labels: true,
            ..RenderOptions::default()
        };
        let svg = render_svg_with_options(&graph, &opts).unwrap();
        assert_eq!(svg.matches("data-floating-band=\"true\"").count(), 1);
        assert_eq!(svg.matches("class=\"bracket\"").count(), 1);
        assert_eq!(svg.matches("class=\"res-label\"").count(), 11);
        assert!(svg.contains(">\u{03B2} ?</text>"), "{svg}");
        // The stem and one stub: no dashed line reaches a candidate parent.
        assert_eq!(svg.matches("class=\"uncertain\"").count(), 2);
        assert_on_canvas(&svg, "floating fragment");

        let (layout, bands) = compute_floating_layout(&graph, graph.root().unwrap());
        assert_eq!(layout.len(), graph.node_count());
        let [band] = bands.try_into().unwrap();
        let main_top = graph
            .subtree(graph.root().unwrap())
            .iter()
            .map(|node| layout[node].y)
            .fold(f64::MAX, f64::min);
        assert!(band.lane < main_top, "bracket above the main structure");
        for node in graph.subtree(NodeIndex::new(9)) {
            assert!(layout[&node].y < band.lane, "fragment above its bracket");
        }
    }

    #[test]
    fn floating_substituents_join_the_band_of_their_candidates() {
        let graph = parse(SUBSTITUENT);
        let opts = RenderOptions {
            floating_fragments: true,
            ..RenderOptions::default()
        };
        let svg = render_svg_with_options(&graph, &opts).unwrap();
        assert_eq!(svg.matches("{Ac?}").count(), 1);
        assert_eq!(svg.matches("class=\"bracket\"").count(), 1);
        assert_on_canvas(&svg, "floating substituent");
        // The label's box stays on the canvas too.
        let body = &svg[svg.find("</style>").unwrap()..];
        assert!(attributes(body, &["x"]).into_iter().all(|x| x >= 0.0));
    }

    #[test]
    fn a_glycan_with_nothing_undefined_ignores_the_floating_option() {
        let graph = iupac(N_GLYCAN);
        let opts = RenderOptions {
            floating_fragments: true,
            repeat_brackets: true,
            ..RenderOptions::default()
        };
        assert_eq!(
            render_svg_with_options(&graph, &opts).unwrap(),
            render_svg(&graph).unwrap()
        );
    }

    #[test]
    fn highlight_selections_build_from_residues_and_matches() {
        let graph = iupac("Fuc(a1-3)GlcNAc(b1-4)GlcNAc");
        let selection = HighlightSelection::from_residues(&graph, [0, 1]).unwrap();
        assert_eq!(selection.node_indices, BTreeSet::from([0, 1]));
        assert_eq!(selection.edge_indices.len(), 1);

        // Residues that are not bonded select no bond.
        let apart = HighlightSelection::from_residues(&graph, [0, 2]).unwrap();
        assert!(apart.edge_indices.is_empty());

        assert!(matches!(
            HighlightSelection::from_residues(&graph, [3]),
            Err(SnfgError::InvalidResidueIndex { index: 3, count: 3 })
        ));

        let motif = iupac("Fuc(a1-?)GlcNAc");
        let matches = find_motif_matches(&graph, &motif).unwrap();
        let from_matches = HighlightSelection::from_matches(&matches);
        assert_eq!(from_matches.node_indices.len(), 2);
        assert_eq!(
            render_svg_with_selection(&graph, &from_matches, &RenderOptions::default()).unwrap(),
            render_svg_with_motifs(&graph, &[motif], &RenderOptions::default()).unwrap()
        );
    }

    #[test]
    fn cleavage_marks_hook_towards_the_kept_fragment() {
        let graph = iupac("Gal(b1-4)GlcNAc");
        let mark = |side, label: &str| CleavageMark {
            edge: 0,
            side,
            label: label.into(),
        };
        let opts = RenderOptions {
            cleavages: vec![
                mark(CleavageSide::NonReducing, "B1"),
                mark(CleavageSide::Reducing, "Y1"),
            ],
            ..RenderOptions::default()
        };
        let svg = render_svg_with_options(&graph, &opts).unwrap();
        assert_eq!(svg.matches("class=\"cleavage\"").count(), 2);
        assert!(svg.contains(">B1</text>") && svg.contains(">Y1</text>"));
        assert!(svg.contains("data-kept=\"non-reducing\""));
        assert!(svg.contains("data-kept=\"reducing\""));
        for orientation in [Orientation::RootBottom, Orientation::RootTop] {
            let opts = RenderOptions {
                orientation,
                compact: true,
                ..opts.clone()
            };
            let svg = render_svg_with_options(&graph, &opts).unwrap();
            assert_on_canvas(&svg, "cleavage");
            let body = &svg[svg.find("</style>").unwrap()..];
            let width = canvas_size(&svg).0;
            // Labels start or end at x; the room for them was reserved.
            for x in attributes(body, &["x"]) {
                let room = x.min(width - x);
                assert!(room >= CLEAVAGE_LABEL_SIZE, "label at {x} of {width}");
            }
        }

        let bad = RenderOptions {
            cleavages: vec![CleavageMark {
                edge: 7,
                side: CleavageSide::Reducing,
                label: "Y1".into(),
            }],
            ..RenderOptions::default()
        };
        assert!(matches!(
            render_svg_with_options(&graph, &bad),
            Err(SnfgError::InvalidBondIndex { index: 7, count: 1 })
        ));
    }

    #[test]
    fn dimmed_residue_stays_opaque_over_its_linkage() {
        let graph = crabwurcs_iupac::parse_iupac_condensed("Man(b1-4)GlcNAc").unwrap();
        let motif = crabwurcs_iupac::parse_iupac_condensed("GlcNAc").unwrap();
        let png = render_png_with_motifs(&graph, &[motif], &RenderOptions::default()).unwrap();

        let decoder = png::Decoder::new(std::io::Cursor::new(&png));
        let mut reader = decoder.read_info().unwrap();
        let mut pixels = vec![0; reader.output_buffer_size().unwrap()];
        let info = reader.next_frame(&mut pixels).unwrap();

        // The linkage terminates at the centre of the dimmed Man symbol. Its
        // muted green fill must fully occlude that line.
        let x = 140usize;
        let y = 140usize;
        let offset = (y * info.width as usize + x) * 4;
        let pixel = &pixels[offset..offset + 4];
        assert_eq!(pixel[3], 255, "dimmed symbols must remain fully opaque");
        assert_eq!(
            pixel,
            [0xCD, 0xE9, 0xDF, 0xFF],
            "dimmed Man should use GlycoDraw's exact muted green"
        );
    }
}
