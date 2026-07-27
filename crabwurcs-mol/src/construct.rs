use std::collections::{HashMap, HashSet};

use chematic::core::{
    Atom, AtomIdx, BondOrder, Chirality, Element, Molecule, MoleculeBuilder, STEREO_H_SENTINEL,
};
use crabwurcs_core::{AnomericSymbol, Monosaccharide, RepeatCount, ResidueGraph, RingClosure};
use petgraph::visit::EdgeRef;

use crate::map::MapGraph;
use crate::{MolError, MolResult, NonConcreteFeature};

#[derive(Debug)]
struct BuiltResidue {
    carbons: Vec<AtomIdx>,
    descriptors: Vec<char>,
    primary_modification: Vec<Option<AtomIdx>>,
    anomeric_oxygen: AtomIdx,
    ring_oxygen: Option<AtomIdx>,
    anomeric_position: usize,
}

fn carbon_descriptors(residue: &Monosaccharide) -> MolResult<Vec<char>> {
    let mut descriptors = Vec::new();
    let prefix = residue.anomeric_prefix.as_str();
    if prefix.starts_with("Aa") || prefix.starts_with("AU") {
        descriptors.extend(['A', prefix.chars().nth(1).unwrap_or('U'), 'd']);
    } else if prefix.starts_with("ha") || prefix.starts_with("hU") {
        descriptors.extend(['h', prefix.chars().nth(1).unwrap_or('U')]);
    } else {
        descriptors.push(prefix.chars().next().unwrap_or('u'));
    }
    descriptors.extend(residue.skeleton_code.chars());
    if descriptors.len() < 4 {
        return Err(MolError::UnsupportedChemistry(format!(
            "backbone `{prefix}{}` is too short",
            residue.skeleton_code
        )));
    }
    Ok(descriptors)
}

fn add_bond(
    builder: &mut MoleculeBuilder,
    left: AtomIdx,
    right: AtomIdx,
    order: BondOrder,
) -> MolResult<()> {
    builder
        .add_bond(left, right, order)
        .map(|_| ())
        .map_err(|error| MolError::UnsupportedChemistry(error.to_string()))
}

fn add_oxygen_with_order(
    builder: &mut MoleculeBuilder,
    atom: AtomIdx,
    order: BondOrder,
) -> MolResult<AtomIdx> {
    let oxygen = builder.add_atom(Atom::new(Element::O));
    add_bond(builder, atom, oxygen, order)?;
    Ok(oxygen)
}

fn backbone_bond_order(left: char, right: char) -> BondOrder {
    if left == 'K' || matches!(right, 'n' | 'N' | 'e' | 'z' | 'f' | 'E' | 'Z' | 'F' | 'K') {
        BondOrder::Double
    } else if matches!(right, 't' | 'T') {
        BondOrder::Triple
    } else {
        BondOrder::Single
    }
}

fn default_oxygen_bonds(descriptor: char, terminal: bool) -> MolResult<Vec<BondOrder>> {
    use BondOrder::{Double, Single, Triple};
    let bonds = match descriptor {
        'm' | 'd' | 'n' | 't' => vec![],
        'h' => vec![Single],
        'c' => vec![Single, Single],
        'C' => vec![Single; if terminal { 3 } else { 2 }],
        '1' | '2' | '3' | '4' | 'x' => vec![Single; if terminal { 2 } else { 1 }],
        '5' | '6' | '7' | '8' | 'X' => vec![Single; if terminal { 3 } else { 2 }],
        'o' | 'O' | 'K' => vec![Double],
        'A' => vec![Double, Single],
        'N' => vec![Single; if terminal { 2 } else { 1 }],
        'e' | 'z' | 'f' => terminal.then_some(Single).into_iter().collect(),
        'E' | 'Z' | 'F' => vec![Single; if terminal { 2 } else { 1 }],
        'T' => vec![Triple],
        // Anomeric descriptors receive their ring and anomeric oxygens
        // explicitly. `u`/`U` are retained for the common reducing-end
        // convention and conservatively use one unspecified oxygen.
        'a' => vec![Single, Single],
        'u' | 'U' => vec![Single],
        'Q' | '?' => {
            return Err(MolError::InvalidSkeleton {
                offset: 0,
                message: format!("descriptor `{descriptor}` does not define atom valence"),
            });
        }
        other => {
            return Err(MolError::InvalidSkeleton {
                offset: 0,
                message: format!("unknown carbon descriptor `{other}`"),
            });
        }
    };
    Ok(bonds)
}

fn add_modification(
    builder: &mut MoleculeBuilder,
    backbone: AtomIdx,
    descriptor: &str,
) -> MolResult<AtomIdx> {
    let graph = MapGraph::parse(&format!("*{descriptor}"))?;
    let built = graph.build(builder)?;
    let (attachment, order) =
        built
            .attachments
            .get(&1)
            .copied()
            .ok_or_else(|| MolError::InvalidMap {
                offset: 0,
                message: "residue MAP has no primary attachment".into(),
            })?;
    add_bond(builder, backbone, attachment, order)?;
    Ok(attachment)
}

fn inferred_ring_positions(residue: &Monosaccharide, carbon_count: usize) -> (usize, usize) {
    if let (Some(start), Some(end)) = (residue.ring_start, residue.ring_end) {
        return (start as usize, end as usize);
    }
    let start = effective_anomeric_position(residue);
    let span = match residue.ring {
        RingClosure::Furanose => 3,
        RingClosure::Pyranose => 4,
        RingClosure::Open => 0,
        RingClosure::Unknown => {
            if carbon_count.saturating_sub(start) >= 4 {
                4
            } else {
                3
            }
        }
    };
    (start, start + span)
}

fn effective_anomeric_position(residue: &Monosaccharide) -> usize {
    if residue.anomeric_position > 0 {
        usize::from(residue.anomeric_position)
    } else if residue.anomeric_prefix.starts_with("AU") || residue.anomeric_prefix.starts_with("hU")
    {
        2
    } else {
        1
    }
}

fn materialize_exact_repeat(graph: &ResidueGraph) -> MolResult<Option<ResidueGraph>> {
    let repeats = graph
        .inner()
        .edge_references()
        .filter(|edge| edge.weight().repeat.is_some())
        .collect::<Vec<_>>();
    if repeats.is_empty() {
        return Ok(None);
    }
    if repeats.iter().any(|edge| {
        matches!(
            edge.weight().repeat,
            Some(RepeatCount::Unknown | RepeatCount::Range { .. })
        )
    }) {
        return Err(MolError::NonConcrete(NonConcreteFeature::VariableRepeat));
    }
    if repeats.len() != 1 {
        return Err(MolError::UnsupportedChemistry(
            "multiple independent exact repeat units".into(),
        ));
    }
    let repeat = repeats[0];
    let count = match repeat.weight().repeat {
        Some(RepeatCount::Exact(count)) if count > 0 => count as usize,
        Some(RepeatCount::Exact(_)) => {
            return Err(MolError::InvalidSkeleton {
                offset: 0,
                message: "an exact repeat count must be positive".into(),
            });
        }
        _ => unreachable!(),
    };
    let repeat_edge = repeat.id();
    let repeat_source = repeat.source();
    let repeat_target = repeat.target();
    let mut expanded = ResidueGraph::new();
    let mut copies = Vec::with_capacity(count);
    for _ in 0..count {
        let mut nodes = HashMap::new();
        for node in graph.inner().node_indices() {
            let copied = expanded.add_residue(
                graph
                    .residue(node)
                    .expect("source graph node has a residue")
                    .clone(),
            );
            nodes.insert(node, copied);
        }
        copies.push(nodes);
    }
    for nodes in &copies {
        for edge in graph.inner().edge_references() {
            if edge.id() == repeat_edge {
                continue;
            }
            let mut linkage = edge.weight().clone();
            linkage.repeat = None;
            expanded.add_linkage(nodes[&edge.source()], nodes[&edge.target()], linkage);
        }
    }
    for index in 0..count.saturating_sub(1) {
        let mut linkage = repeat.weight().clone();
        linkage.repeat = None;
        linkage.cyclic = false;
        expanded.add_linkage(
            copies[index][&repeat_source],
            copies[index + 1][&repeat_target],
            linkage,
        );
    }
    if let Some(root) = graph.root() {
        expanded.set_root(copies[0][&root]);
    }
    Ok(Some(expanded))
}

pub(crate) fn construct_molecule(graph: &ResidueGraph) -> MolResult<Molecule> {
    if let Some(expanded) = materialize_exact_repeat(graph)? {
        return construct_molecule(&expanded);
    }
    if graph.is_composition() {
        return Err(MolError::NonConcrete(NonConcreteFeature::Composition));
    }
    if !graph.undefined_linkages().is_empty() {
        return Err(MolError::NonConcrete(NonConcreteFeature::UndefinedLinkage));
    }
    if !graph.undefined_modifications().is_empty() {
        return Err(MolError::NonConcrete(
            NonConcreteFeature::UndefinedModification,
        ));
    }

    let mut occupied_positions = HashSet::new();
    for edge in graph.inner().edge_references() {
        let linkage = edge.weight();
        if !linkage.parent_position_alternatives.is_empty()
            || !linkage.child_position_alternatives.is_empty()
        {
            return Err(MolError::NonConcrete(
                NonConcreteFeature::AlternativeLinkagePosition,
            ));
        }
        if linkage.parent_position.0 == 0 || linkage.child_position.0 == 0 {
            return Err(MolError::NonConcrete(
                NonConcreteFeature::UnknownLinkagePosition,
            ));
        }
        if linkage.parent_probability.is_some() || linkage.child_probability.is_some() {
            return Err(MolError::NonConcrete(
                NonConcreteFeature::LinkageProbability,
            ));
        }
        debug_assert!(linkage.repeat.is_none());
        occupied_positions.insert((edge.source().index(), linkage.parent_position.0 as usize));
        if linkage.map_code.is_some() {
            occupied_positions.insert((edge.target().index(), linkage.child_position.0 as usize));
        }
    }

    let mut builder = MoleculeBuilder::new();
    let mut built = HashMap::new();
    for node in graph.inner().node_indices() {
        let residue = graph.residue(node).expect("node has residue");
        let descriptors = carbon_descriptors(residue)?;
        let carbon_count = descriptors.len();
        let anomeric_position = effective_anomeric_position(residue);
        let reference = descriptors
            .iter()
            .copied()
            .filter(|descriptor| matches!(descriptor, '1' | '2'))
            .nth(3)
            .or_else(|| {
                descriptors
                    .iter()
                    .rev()
                    .copied()
                    .find(|descriptor| matches!(descriptor, '1' | '2'))
            });
        let anomeric_absolute = match (residue.anomeric_symbol, reference) {
            (AnomericSymbol::Beta, reference) => reference,
            (AnomericSymbol::Alpha, Some('1')) => Some('2'),
            (AnomericSymbol::Alpha, Some('2')) => Some('1'),
            _ => None,
        };

        let mut carbons = Vec::with_capacity(carbon_count);
        for (offset, descriptor) in descriptors.iter().copied().enumerate() {
            let position = offset + 1;
            let stereo = if position == anomeric_position {
                anomeric_absolute
            } else {
                matches!(descriptor, '1' | '2').then_some(descriptor)
            };
            let chirality = match stereo {
                Some('1') => Chirality::Clockwise,
                Some('2') => Chirality::CounterClockwise,
                _ => Chirality::None,
            };
            let atom = if chirality == Chirality::None {
                Atom::new(Element::C)
            } else {
                Atom::bracket(
                    Element::C,
                    None,
                    chirality,
                    u8::from(position != anomeric_position || anomeric_position == 1),
                    0,
                    None,
                )
            };
            carbons.push(builder.add_atom(atom));
        }
        for (offset, pair) in carbons.windows(2).enumerate() {
            add_bond(
                &mut builder,
                pair[0],
                pair[1],
                backbone_bond_order(descriptors[offset], descriptors[offset + 1]),
            )?;
        }

        let mut primary_modification = vec![None; carbon_count];
        let open_chain = residue.ring == RingClosure::Open;
        let (ring_start, ring_end) = inferred_ring_positions(residue, carbon_count);
        if anomeric_position == 0 || anomeric_position > carbon_count {
            return Err(MolError::InvalidSkeleton {
                offset: 0,
                message: format!(
                    "anomeric position {anomeric_position} on {carbon_count}-carbon backbone"
                ),
            });
        }
        let ring_oxygen = if open_chain {
            None
        } else {
            if ring_start == 0 || ring_end > carbon_count || ring_start == ring_end {
                return Err(MolError::InvalidSkeleton {
                    offset: 0,
                    message: format!(
                        "ring {ring_start}-{ring_end} on {carbon_count}-carbon backbone"
                    ),
                });
            }
            let oxygen = builder.add_atom(Atom::new(Element::O));
            add_bond(
                &mut builder,
                carbons[ring_start - 1],
                oxygen,
                BondOrder::Single,
            )?;
            add_bond(
                &mut builder,
                carbons[ring_end - 1],
                oxygen,
                BondOrder::Single,
            )?;
            Some(oxygen)
        };

        let anomeric_oxygen = add_oxygen_with_order(
            &mut builder,
            carbons[anomeric_position - 1],
            if open_chain {
                BondOrder::Double
            } else {
                BondOrder::Single
            },
        )?;
        primary_modification[anomeric_position - 1] = Some(anomeric_oxygen);
        if let Some(ring_oxygen) = ring_oxygen {
            primary_modification[ring_end - 1] = Some(ring_oxygen);
        }
        let mut consumed_modifications = vec![0usize; carbon_count];
        if !open_chain {
            consumed_modifications[ring_start - 1] += 1;
            consumed_modifications[ring_end - 1] += 1;
        }
        consumed_modifications[anomeric_position - 1] += 1;

        let modifications = residue
            .modifications
            .iter()
            .map(|modification| {
                if modification.position.0 == 0 || modification.probability.is_some() {
                    return Err(MolError::NonConcrete(
                        NonConcreteFeature::UndefinedModification,
                    ));
                }
                Ok((modification.position.0 as usize, &modification.descriptor))
            })
            .collect::<MolResult<HashMap<_, _>>>()?;
        for (offset, descriptor) in descriptors.iter().copied().enumerate() {
            let position = offset + 1;
            if let Some(map) = modifications.get(&position) {
                let atom = add_modification(&mut builder, carbons[offset], map)?;
                primary_modification[offset].get_or_insert(atom);
                consumed_modifications[offset] += 1;
            }
            if occupied_positions.contains(&(node.index(), position)) {
                consumed_modifications[offset] += 1;
            }
            let terminal = offset == 0 || offset + 1 == carbon_count;
            let default_bonds = default_oxygen_bonds(descriptor, terminal)?;
            for order in default_bonds
                .into_iter()
                .skip(consumed_modifications[offset])
            {
                let oxygen = add_oxygen_with_order(&mut builder, carbons[offset], order)?;
                primary_modification[offset].get_or_insert(oxygen);
            }
        }
        built.insert(
            node.index(),
            BuiltResidue {
                carbons,
                descriptors,
                primary_modification,
                anomeric_oxygen,
                ring_oxygen,
                anomeric_position,
            },
        );
    }

    for edge in graph.inner().edge_references() {
        let linkage = edge.weight();
        if let Some(map_code) = linkage.map_code.as_deref() {
            let map = MapGraph::parse(map_code)?;
            let emitted = map.build(&mut builder)?;
            let parent_star = linkage.parent_modification_position.unwrap_or(1);
            let child_star = linkage.child_modification_position.unwrap_or_else(|| {
                if emitted.attachments.contains_key(&2) {
                    2
                } else {
                    1
                }
            });
            let (parent_atom, parent_order) = emitted
                .attachments
                .get(&parent_star)
                .copied()
                .ok_or_else(|| MolError::InvalidMap {
                    offset: 0,
                    message: format!("MAP has no parent attachment star {parent_star}"),
                })?;
            let (child_atom, child_order) = emitted
                .attachments
                .get(&child_star)
                .copied()
                .ok_or_else(|| MolError::InvalidMap {
                    offset: 0,
                    message: format!("MAP has no child attachment star {child_star}"),
                })?;
            let parent =
                built
                    .get(&edge.source().index())
                    .ok_or_else(|| MolError::InvalidSkeleton {
                        offset: 0,
                        message: "missing parent residue".into(),
                    })?;
            let child =
                built
                    .get(&edge.target().index())
                    .ok_or_else(|| MolError::InvalidSkeleton {
                        offset: 0,
                        message: "missing child residue".into(),
                    })?;
            let parent_carbon = *parent
                .carbons
                .get(linkage.parent_position.0 as usize - 1)
                .ok_or_else(|| MolError::InvalidSkeleton {
                    offset: 0,
                    message: "parent linkage position is outside its backbone".into(),
                })?;
            let child_carbon = *child
                .carbons
                .get(linkage.child_position.0 as usize - 1)
                .ok_or_else(|| MolError::InvalidSkeleton {
                    offset: 0,
                    message: "child linkage position is outside its backbone".into(),
                })?;
            add_bond(&mut builder, parent_carbon, parent_atom, parent_order)?;
            add_bond(&mut builder, child_carbon, child_atom, child_order)?;
            continue;
        }
        let child_position = linkage.child_position.0 as usize;
        let donor = built
            .get(&edge.target().index())
            .and_then(|residue| {
                residue
                    .primary_modification
                    .get(child_position.wrapping_sub(1))
                    .copied()
                    .flatten()
            })
            .ok_or_else(|| MolError::InvalidSkeleton {
                offset: 0,
                message: format!("missing donor modification at position {child_position}"),
            })?;
        let parent =
            built
                .get_mut(&edge.source().index())
                .ok_or_else(|| MolError::InvalidSkeleton {
                    offset: 0,
                    message: "missing acceptor residue".into(),
                })?;
        let position = linkage.parent_position.0 as usize;
        let carbon = *parent
            .carbons
            .get(position.wrapping_sub(1))
            .ok_or_else(|| MolError::InvalidSkeleton {
                offset: 0,
                message: format!("acceptor position {position}"),
            })?;
        add_bond(&mut builder, carbon, donor, BondOrder::Single)?;
        parent.primary_modification[position - 1] = Some(donor);
    }

    for residue in built.values() {
        for position in 1..residue.carbons.len().saturating_sub(1) {
            if position + 1 == residue.anomeric_position
                || !matches!(residue.descriptors[position], '1' | '2')
            {
                continue;
            }
            let Some(modification) = residue.primary_modification[position] else {
                continue;
            };
            builder.set_stereo_neighbor_order(
                residue.carbons[position],
                vec![
                    residue.carbons[position - 1].0,
                    modification.0,
                    residue.carbons[position + 1].0,
                    STEREO_H_SENTINEL,
                ],
            );
        }
        let anomer = residue.anomeric_position - 1;
        if builder.atom_at(residue.carbons[anomer]).chirality == Chirality::None {
            continue;
        }
        let Some(ring_oxygen) = residue.ring_oxygen else {
            continue;
        };
        let previous = if anomer == 0 {
            STEREO_H_SENTINEL
        } else {
            residue.carbons[anomer - 1].0
        };
        builder.set_stereo_neighbor_order(
            residue.carbons[anomer],
            vec![
                previous,
                residue.anomeric_oxygen.0,
                residue.carbons[anomer + 1].0,
                ring_oxygen.0,
            ],
        );
    }
    Ok(builder.build())
}
