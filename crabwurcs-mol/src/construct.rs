use std::collections::{HashMap, HashSet};

use chematic::core::{
    Atom, AtomIdx, BondOrder, Chirality, Element, Molecule, MoleculeBuilder, STEREO_H_SENTINEL,
};
use crabwurcs_core::{AnomericSymbol, Monosaccharide, ResidueGraph, RingClosure};
use petgraph::visit::EdgeRef;

use crate::{MolError, MolResult};

#[derive(Debug)]
struct BuiltResidue {
    carbons: Vec<AtomIdx>,
    descriptors: Vec<char>,
    primary_modification: Vec<Option<AtomIdx>>,
    anomeric_oxygen: AtomIdx,
    ring_oxygen: AtomIdx,
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

fn add_oxygen(builder: &mut MoleculeBuilder, atom: AtomIdx) -> MolResult<AtomIdx> {
    let oxygen = builder.add_atom(Atom::new(Element::O));
    add_bond(builder, atom, oxygen, BondOrder::Single)?;
    Ok(oxygen)
}

fn add_carboxyl(builder: &mut MoleculeBuilder, carbon: AtomIdx) -> MolResult<AtomIdx> {
    let hydroxyl = add_oxygen(builder, carbon)?;
    let carbonyl = builder.add_atom(Atom::new(Element::O));
    add_bond(builder, carbon, carbonyl, BondOrder::Double)?;
    Ok(hydroxyl)
}

fn add_sulfate(
    builder: &mut MoleculeBuilder,
    backbone: AtomIdx,
    through_nitrogen: bool,
) -> MolResult<AtomIdx> {
    let primary = builder.add_atom(Atom::new(if through_nitrogen {
        Element::N
    } else {
        Element::O
    }));
    add_bond(builder, backbone, primary, BondOrder::Single)?;
    let sulfur = builder.add_atom(Atom::new(Element::S));
    add_bond(builder, primary, sulfur, BondOrder::Single)?;
    for order in [BondOrder::Double, BondOrder::Double, BondOrder::Single] {
        let oxygen = builder.add_atom(Atom::new(Element::O));
        add_bond(builder, sulfur, oxygen, order)?;
    }
    Ok(primary)
}

fn add_n_acyl(
    builder: &mut MoleculeBuilder,
    backbone: AtomIdx,
    glycolyl: bool,
) -> MolResult<AtomIdx> {
    let nitrogen = builder.add_atom(Atom::new(Element::N));
    add_bond(builder, backbone, nitrogen, BondOrder::Single)?;
    let carbonyl = builder.add_atom(Atom::new(Element::C));
    add_bond(builder, nitrogen, carbonyl, BondOrder::Single)?;
    let oxygen = builder.add_atom(Atom::new(Element::O));
    add_bond(builder, carbonyl, oxygen, BondOrder::Double)?;
    let side = builder.add_atom(Atom::new(Element::C));
    add_bond(builder, carbonyl, side, BondOrder::Single)?;
    if glycolyl {
        add_oxygen(builder, side)?;
    }
    Ok(nitrogen)
}

fn add_oxygen_modification(
    builder: &mut MoleculeBuilder,
    backbone: AtomIdx,
    descriptor: &str,
) -> MolResult<AtomIdx> {
    let oxygen = add_oxygen(builder, backbone)?;
    match descriptor {
        "OC" => {
            let methyl = builder.add_atom(Atom::new(Element::C));
            add_bond(builder, oxygen, methyl, BondOrder::Single)?;
        }
        "OCC/3=O" => {
            let carbonyl = builder.add_atom(Atom::new(Element::C));
            add_bond(builder, oxygen, carbonyl, BondOrder::Single)?;
            let carbonyl_oxygen = builder.add_atom(Atom::new(Element::O));
            add_bond(builder, carbonyl, carbonyl_oxygen, BondOrder::Double)?;
            let methyl = builder.add_atom(Atom::new(Element::C));
            add_bond(builder, carbonyl, methyl, BondOrder::Single)?;
        }
        "OP^XOCCNC/7C/7C/3O/3=O" => {
            let phosphorus = builder.add_atom(Atom::new(Element::P));
            add_bond(builder, oxygen, phosphorus, BondOrder::Single)?;
            let phosphoryl = builder.add_atom(Atom::new(Element::O));
            add_bond(builder, phosphorus, phosphoryl, BondOrder::Double)?;
            add_oxygen(builder, phosphorus)?;
            let bridge = builder.add_atom(Atom::new(Element::O));
            add_bond(builder, phosphorus, bridge, BondOrder::Single)?;
            let first = builder.add_atom(Atom::new(Element::C));
            let second = builder.add_atom(Atom::new(Element::C));
            let mut nitrogen = Atom::new(Element::N);
            nitrogen.charge = 1;
            let nitrogen = builder.add_atom(nitrogen);
            add_bond(builder, bridge, first, BondOrder::Single)?;
            add_bond(builder, first, second, BondOrder::Single)?;
            add_bond(builder, second, nitrogen, BondOrder::Single)?;
            for _ in 0..3 {
                let methyl = builder.add_atom(Atom::new(Element::C));
                add_bond(builder, nitrogen, methyl, BondOrder::Single)?;
            }
        }
        _ => {
            return Err(MolError::UnsupportedChemistry(format!(
                "MAP `{descriptor}`"
            )));
        }
    }
    Ok(oxygen)
}

fn add_modification(
    builder: &mut MoleculeBuilder,
    backbone: AtomIdx,
    descriptor: &str,
) -> MolResult<AtomIdx> {
    match descriptor {
        "NCC/3=O" => add_n_acyl(builder, backbone, false),
        "NCCO/3=O" => add_n_acyl(builder, backbone, true),
        "NSO/3=O/3=O" => add_sulfate(builder, backbone, true),
        "OSO/3=O/3=O" => add_sulfate(builder, backbone, false),
        value if value.starts_with('O') => add_oxygen_modification(builder, backbone, value),
        _ => Err(MolError::UnsupportedChemistry(format!(
            "MAP `{descriptor}`"
        ))),
    }
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

pub(crate) fn construct_molecule(graph: &ResidueGraph) -> MolResult<Molecule> {
    if graph.is_composition()
        || !graph.undefined_linkages().is_empty()
        || !graph.undefined_modifications().is_empty()
    {
        return Err(MolError::UnsupportedChemistry(
            "compositions or undefined fragments".into(),
        ));
    }

    let mut occupied_acceptors = HashSet::new();
    for edge in graph.inner().edge_references() {
        let linkage = edge.weight();
        if !linkage.parent_position_alternatives.is_empty()
            || !linkage.child_position_alternatives.is_empty()
            || linkage.repeat.is_some()
            || linkage.map_code.is_some()
        {
            return Err(MolError::UnsupportedChemistry(
                "ambiguous, repeating, or MAP-bridged linkages".into(),
            ));
        }
        occupied_acceptors.insert((edge.source().index(), linkage.parent_position.0 as usize));
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
        for pair in carbons.windows(2) {
            add_bond(&mut builder, pair[0], pair[1], BondOrder::Single)?;
        }

        let (ring_start, ring_end) = inferred_ring_positions(residue, carbon_count);
        if ring_start == 0 || ring_end > carbon_count || ring_start == ring_end {
            return Err(MolError::UnsupportedChemistry(format!(
                "ring {ring_start}-{ring_end} on {carbon_count}-carbon backbone"
            )));
        }
        let ring_oxygen = builder.add_atom(Atom::new(Element::O));
        add_bond(
            &mut builder,
            carbons[ring_start - 1],
            ring_oxygen,
            BondOrder::Single,
        )?;
        add_bond(
            &mut builder,
            carbons[ring_end - 1],
            ring_oxygen,
            BondOrder::Single,
        )?;

        let mut primary_modification = vec![None; carbon_count];
        let anomeric_oxygen = add_oxygen(&mut builder, carbons[anomeric_position - 1])?;
        primary_modification[anomeric_position - 1] = Some(anomeric_oxygen);
        primary_modification[ring_end - 1] = Some(ring_oxygen);

        let modifications = residue
            .modifications
            .iter()
            .map(|modification| (modification.position.0 as usize, &modification.descriptor))
            .collect::<HashMap<_, _>>();
        for (offset, descriptor) in descriptors.iter().copied().enumerate() {
            let position = offset + 1;
            if position == anomeric_position || position == ring_end {
                continue;
            }
            if let Some(map) = modifications.get(&position) {
                primary_modification[offset] =
                    Some(add_modification(&mut builder, carbons[offset], map)?);
                continue;
            }
            if occupied_acceptors.contains(&(node.index(), position)) {
                continue;
            }
            primary_modification[offset] = match descriptor {
                '1' | '2' | 'x' | 'h' => Some(add_oxygen(&mut builder, carbons[offset])?),
                'A' => Some(add_carboxyl(&mut builder, carbons[offset])?),
                'm' | 'd' => None,
                _ => None,
            };
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
        let donor = built
            .get(&edge.target().index())
            .ok_or_else(|| MolError::UnsupportedChemistry("missing donor residue".into()))?
            .anomeric_oxygen;
        let parent = built
            .get_mut(&edge.source().index())
            .ok_or_else(|| MolError::UnsupportedChemistry("missing acceptor residue".into()))?;
        let position = linkage.parent_position.0 as usize;
        let carbon = *parent
            .carbons
            .get(position.wrapping_sub(1))
            .ok_or_else(|| {
                MolError::UnsupportedChemistry(format!("acceptor position {position}"))
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
                residue.ring_oxygen.0,
            ],
        );
    }
    Ok(builder.build())
}
