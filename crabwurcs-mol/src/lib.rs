//! WURCS ⇄ chemical structure format (MOL/SDF/SMILES) conversion.
//!
//! The always-available pure-Rust backend recognizes 938 stereochemically
//! specified GlycoShape structures bundled with crabWURCS: 838 with source
//! WURCS and 100 whose WURCS is derived from supplied IUPAC. SMILES are matched
//! by canonical molecular graph, not by string spelling, so an equivalent
//! traversal/order is accepted.  MOL and SDF records written here contain a
//! crabWURCS canonical-SMILES provenance marker because chematic 0.4's MOL
//! writer does not yet serialize all atom-centred SMILES stereochemistry as
//! CTfile wedges/parity.  That marker makes crabWURCS's own file round-trip
//! lossless while still producing a standards-shaped MOL/SDF record.
//!
//! Previously unseen glycans are extracted from their molecular graph, and
//! finite defined WURCS graphs are constructed as real atom/bond graphs rather
//! than requiring a corpus hit. Unknown stereochemistry remains explicit.

use std::collections::HashMap;
use std::sync::OnceLock;

use chematic::core::{AtomIdx, Chirality, Molecule, STEREO_H_SENTINEL};
use chematic::mol::MolMetadata;
use crabwurcs_core::{
    parse_wurcs, write_wurcs, AnomericSymbol, CarbonPosition, Linkage, Modification,
    Monosaccharide, ResidueGraph, RingClosure,
};
use thiserror::Error;

mod construct;

const SOURCE_CORPUS: &str = include_str!("../data/glycoshape_notations.tsv");
const DERIVED_CORPUS: &str = include_str!("../data/glycoshape_derived_notations.tsv");
const CANONICAL_INDEX: &str = include_str!("../data/glycoshape_canonical_smiles.tsv");
const SMILES_MARKER: &str = "crabWURCS canonical-SMILES=";

#[derive(Debug, Error)]
pub enum MolError {
    #[error("the optional RDKit backend is unavailable")]
    BackendDisabled,

    #[error("failed to parse chemical structure: {0}")]
    ParseError(String),

    #[error("an SDF conversion requires exactly one record, but found {0}")]
    SdfRecordCount(usize),

    #[error("multiple different corpus glycans of the same size match the molecule")]
    AmbiguousGlycanFound,

    #[error("no glycan structure could be extracted from the input molecule")]
    NoGlycanFound,

    #[error("WURCS chemistry is not yet constructible: {0}")]
    UnsupportedChemistry(String),

    #[error(transparent)]
    Core(#[from] crabwurcs_core::CoreError),
}

pub type MolResult<T> = Result<T, MolError>;

/// Which chemical file format is being read/written.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ChemFormat {
    Mol,
    Sdf,
    Smiles,
}

fn corpus_pairs() -> impl Iterator<Item = (&'static str, &'static str)> {
    SOURCE_CORPUS
        .lines()
        .chain(DERIVED_CORPUS.lines())
        .filter_map(|line| {
            let mut fields = line.split('\t');
            let wurcs = fields.next()?;
            let smiles = fields.nth(3)?;
            Some((wurcs, smiles))
        })
}

fn canonical_smiles(molecule: &Molecule) -> String {
    chematic::smiles::canonical_smiles(molecule)
}

fn parse_smiles(input: &str) -> MolResult<Molecule> {
    chematic::smiles::parse(input.trim()).map_err(|error| MolError::ParseError(error.to_string()))
}

fn canonical_corpus() -> &'static HashMap<String, &'static str> {
    static INDEX: OnceLock<HashMap<String, &'static str>> = OnceLock::new();
    INDEX.get_or_init(|| {
        CANONICAL_INDEX
            .lines()
            .filter_map(|line| {
                let (wurcs, canonical) = line.split_once('\t')?;
                Some((canonical.to_owned(), wurcs))
            })
            .collect()
    })
}

fn corpus_smiles_for_wurcs(wurcs: &str) -> Option<&'static str> {
    corpus_pairs().find_map(|(candidate, smiles)| (candidate == wurcs).then_some(smiles))
}

fn exact_wurcs_for_molecule(molecule: &Molecule) -> Option<&'static str> {
    canonical_corpus().get(&canonical_smiles(molecule)).copied()
}

fn component_atoms_without_bond(
    molecule: &Molecule,
    start: chematic::core::AtomIdx,
    removed_bond: chematic::core::BondIdx,
) -> Vec<bool> {
    let mut included = vec![false; molecule.atom_count()];
    let mut stack = vec![start];
    included[start.0 as usize] = true;
    while let Some(atom) = stack.pop() {
        for (neighbor, bond) in molecule.neighbors(atom) {
            if bond == removed_bond || included[neighbor.0 as usize] {
                continue;
            }
            included[neighbor.0 as usize] = true;
            stack.push(neighbor);
        }
    }
    included
}

fn component_molecule(
    molecule: &Molecule,
    included: &[bool],
    removed_bond: chematic::core::BondIdx,
) -> Molecule {
    let mut component = molecule.with_bond_removed(removed_bond);
    for index in (0..included.len()).rev() {
        if !included[index] {
            component.remove_atom(chematic::core::AtomIdx(index as u32));
        }
    }
    component
}

fn contains_sugar_like_ring(molecule: &Molecule, included: &[bool]) -> bool {
    chematic::perception::find_sssr(molecule)
        .rings()
        .iter()
        .any(|ring| {
            (ring.len() == 5 || ring.len() == 6)
                && ring.iter().all(|atom| included[atom.0 as usize])
                && ring
                    .iter()
                    .filter(|atom| molecule.atom(**atom).element.symbol() == "C")
                    .count()
                    >= 4
                && ring
                    .iter()
                    .any(|atom| molecule.atom(*atom).element.symbol() == "O")
        })
}

/// Recognize an intact corpus glycan attached to one aglycone by one bond.
/// Internal glycosidic cuts are excluded when both sides retain a sugar-like
/// five- or six-membered ring, preventing a known antenna fragment from being
/// returned as though it were the complete unknown glycan.
fn embedded_corpus_wurcs(molecule: &Molecule) -> MolResult<Option<&'static str>> {
    let mut best: Option<(usize, &'static str)> = None;
    for (bond_index, bond) in molecule.bonds() {
        let left = component_atoms_without_bond(molecule, bond.atom1, bond_index);
        let right = component_atoms_without_bond(molecule, bond.atom2, bond_index);
        if left.iter().all(|included| *included) {
            continue; // Cutting a ring bond did not separate the molecule.
        }
        if contains_sugar_like_ring(molecule, &left) && contains_sugar_like_ring(molecule, &right) {
            continue;
        }

        for included in [&left, &right] {
            let atom_count = included.iter().filter(|value| **value).count();
            if atom_count < 5 || best.is_some_and(|(size, _)| atom_count < size) {
                continue;
            }
            let component = component_molecule(molecule, included, bond_index);
            let Some(wurcs) = exact_wurcs_for_molecule(&component) else {
                continue;
            };
            match best {
                Some((size, previous)) if size == atom_count && previous != wurcs => {
                    return Err(MolError::AmbiguousGlycanFound);
                }
                Some((size, _)) if size >= atom_count => {}
                _ => best = Some((atom_count, wurcs)),
            }
        }
    }
    Ok(best.map(|(_, wurcs)| wurcs))
}

#[derive(Debug, Clone)]
struct DetectedSugarRing {
    /// Ring carbons in backbone order (C1..C5 pyranose, C1..C4 furanose;
    /// C2..C6 for a 2-ketulosonic acid).
    carbons: Vec<AtomIdx>,
    anomeric_oxygen: AtomIdx,
    /// A carboxyl carbon preceding the cyclic anomeric carbon identifies a
    /// 2-ketulosonic acid (Kdo/Neu family).
    carboxyl_c1: Option<AtomIdx>,
    /// A hydroxymethyl C1 preceding the cyclic anomeric carbon identifies a
    /// 2-ketose (Fru/Tag/Sor family).
    ketose_c1: Option<AtomIdx>,
    /// Carbon tail after the last ring carbon: C6 for an aldopyranose, or
    /// C7..C9 for a nonulosonic acid.
    tail: Vec<AtomIdx>,
}

fn atom_is(molecule: &Molecule, atom: AtomIdx, symbol: &str) -> bool {
    molecule.atom(atom).element.symbol() == symbol
}

fn exocyclic_neighbors(molecule: &Molecule, atom: AtomIdx, ring: &[AtomIdx]) -> Vec<AtomIdx> {
    molecule
        .neighbors(atom)
        .map(|(neighbor, _)| neighbor)
        .filter(|neighbor| !ring.contains(neighbor))
        .collect()
}

fn is_carboxyl_carbon(molecule: &Molecule, atom: AtomIdx) -> bool {
    if !atom_is(molecule, atom, "C") {
        return false;
    }
    let oxygen_orders = molecule
        .neighbors(atom)
        .filter(|(neighbor, _)| atom_is(molecule, *neighbor, "O"))
        .map(|(_, bond)| molecule.bond(bond).order.order_int())
        .collect::<Vec<_>>();
    oxygen_orders.contains(&1) && oxygen_orders.contains(&2)
}

fn carbon_tail(
    molecule: &Molecule,
    start: AtomIdx,
    excluded: &[AtomIdx],
    maximum: usize,
) -> Vec<AtomIdx> {
    let mut result = Vec::new();
    let mut previous = start;
    while result.len() < maximum {
        let current = result.last().copied().unwrap_or(start);
        let Some(next) = molecule
            .neighbors(current)
            .map(|(neighbor, _)| neighbor)
            .find(|neighbor| {
                atom_is(molecule, *neighbor, "C")
                    && *neighbor != previous
                    && !excluded.contains(neighbor)
                    && !result.contains(neighbor)
            })
        else {
            break;
        };
        previous = current;
        result.push(next);
    }
    result
}

fn detect_sugar_rings(molecule: &Molecule) -> Vec<DetectedSugarRing> {
    chematic::perception::find_sssr(molecule)
        .rings()
        .iter()
        .filter_map(|ring| {
            if ring.len() != 5 && ring.len() != 6 {
                return None;
            }
            let oxygens = ring
                .iter()
                .copied()
                .filter(|atom| atom_is(molecule, *atom, "O"))
                .collect::<Vec<_>>();
            if oxygens.len() != 1
                || ring
                    .iter()
                    .filter(|atom| atom_is(molecule, **atom, "C"))
                    .count()
                    != ring.len() - 1
            {
                return None;
            }

            // Suppress ordinary cyclic ethers: carbohydrate rings expose
            // heteroatom-bearing substituents on most ring carbons.
            let decorated = ring
                .iter()
                .copied()
                .filter(|atom| atom_is(molecule, *atom, "C"))
                .filter(|atom| {
                    exocyclic_neighbors(molecule, *atom, ring)
                        .iter()
                        .any(|neighbor| {
                            matches!(molecule.atom(*neighbor).element.symbol(), "O" | "N" | "S")
                        })
                })
                .count();
            if decorated < 3 {
                return None;
            }

            let ring_oxygen = oxygens[0];
            let adjacent = molecule
                .neighbors(ring_oxygen)
                .map(|(neighbor, _)| neighbor)
                .filter(|neighbor| ring.contains(neighbor) && atom_is(molecule, *neighbor, "C"))
                .collect::<Vec<_>>();
            if adjacent.len() != 2 {
                return None;
            }
            let anomeric = adjacent.iter().copied().find(|carbon| {
                exocyclic_neighbors(molecule, *carbon, ring)
                    .iter()
                    .any(|neighbor| atom_is(molecule, *neighbor, "O"))
            })?;
            let anomeric_oxygen = exocyclic_neighbors(molecule, anomeric, ring)
                .into_iter()
                .find(|neighbor| atom_is(molecule, *neighbor, "O"))?;

            let carbon_count = ring.len() - 1;
            let mut carbons = vec![anomeric];
            let mut previous = ring_oxygen;
            let mut current = anomeric;
            while carbons.len() < carbon_count {
                let next = molecule
                    .neighbors(current)
                    .map(|(neighbor, _)| neighbor)
                    .find(|neighbor| {
                        *neighbor != previous
                            && ring.contains(neighbor)
                            && atom_is(molecule, *neighbor, "C")
                    })?;
                carbons.push(next);
                previous = current;
                current = next;
            }
            let carboxyl_c1 = exocyclic_neighbors(molecule, anomeric, ring)
                .into_iter()
                .find(|neighbor| is_carboxyl_carbon(molecule, *neighbor));
            let ketose_c1 = if carboxyl_c1.is_none() {
                exocyclic_neighbors(molecule, anomeric, ring)
                    .into_iter()
                    .find(|neighbor| atom_is(molecule, *neighbor, "C"))
            } else {
                None
            };
            if carboxyl_c1.is_some() && ring.len() != 6 {
                return None;
            }
            let tail_length = if carboxyl_c1.is_some() {
                3
            } else if ring.len() == 5 {
                2
            } else {
                1
            };
            let tail = carbon_tail(
                molecule,
                *carbons.last().expect("five ring carbons"),
                ring,
                tail_length,
            );
            // A nonulosonic acid has a complete C1..C9 backbone. Reject a
            // partial lookalike rather than assigning incorrect positions.
            if carboxyl_c1.is_some() && tail.len() != 3 {
                return None;
            }
            Some(DetectedSugarRing {
                carbons,
                anomeric_oxygen,
                carboxyl_c1,
                ketose_c1,
                tail,
            })
        })
        .collect()
}

fn residue_and_position_for_atom(
    rings: &[DetectedSugarRing],
    atom: AtomIdx,
) -> Option<(usize, u8)> {
    for (ring_index, ring) in rings.iter().enumerate() {
        if let Some(position) = backbone_atoms(ring)
            .iter()
            .position(|candidate| *candidate == atom)
        {
            return Some((ring_index, position as u8 + 1));
        }
    }
    None
}

fn n_acyl_descriptor(molecule: &Molecule, carbon: AtomIdx) -> Option<&'static str> {
    let nitrogen = molecule
        .neighbors(carbon)
        .map(|(neighbor, _)| neighbor)
        .find(|neighbor| atom_is(molecule, *neighbor, "N"))?;
    let carbonyl = molecule
        .neighbors(nitrogen)
        .map(|(neighbor, _)| neighbor)
        .filter(|candidate| *candidate != carbon && atom_is(molecule, *candidate, "C"))
        .find(|candidate| {
            molecule.neighbors(*candidate).any(|(oxygen, bond)| {
                atom_is(molecule, oxygen, "O") && molecule.bond(bond).order.order_int() == 2
            })
        })?;
    let side_carbon = molecule
        .neighbors(carbonyl)
        .map(|(neighbor, _)| neighbor)
        .find(|neighbor| *neighbor != nitrogen && atom_is(molecule, *neighbor, "C"))?;
    let hydroxylated = molecule
        .neighbors(side_carbon)
        .any(|(neighbor, _)| atom_is(molecule, neighbor, "O"));
    Some(if hydroxylated {
        "NCCO/3=O" // N-glycolyl
    } else {
        "NCC/3=O" // N-acetyl
    })
}

fn n_sulfate_descriptor(molecule: &Molecule, carbon: AtomIdx) -> Option<&'static str> {
    let nitrogen = molecule
        .neighbors(carbon)
        .map(|(neighbor, _)| neighbor)
        .find(|neighbor| atom_is(molecule, *neighbor, "N"))?;
    molecule
        .neighbors(nitrogen)
        .map(|(neighbor, _)| neighbor)
        .any(|neighbor| atom_is(molecule, neighbor, "S"))
        .then_some("NSO/3=O/3=O")
}

fn oxygen_substituent_descriptor(
    molecule: &Molecule,
    carbon: AtomIdx,
    all_backbone_atoms: &[AtomIdx],
) -> Option<&'static str> {
    for oxygen in molecule
        .neighbors(carbon)
        .map(|(neighbor, _)| neighbor)
        .filter(|neighbor| atom_is(molecule, *neighbor, "O"))
    {
        let Some(other) = molecule
            .neighbors(oxygen)
            .map(|(neighbor, _)| neighbor)
            .find(|neighbor| *neighbor != carbon)
        else {
            continue; // ordinary hydroxyl
        };
        if all_backbone_atoms.contains(&other) {
            continue; // ring or glycosidic oxygen
        }
        match molecule.atom(other).element.symbol() {
            "S" => return Some("OSO/3=O/3=O"),
            "P" if is_phosphocholine(molecule, other, oxygen) => {
                return Some("OP^XOCCNC/7C/7C/3O/3=O");
            }
            "C" => {
                let is_carbonyl = molecule.neighbors(other).any(|(neighbor, bond)| {
                    atom_is(molecule, neighbor, "O")
                        && neighbor != oxygen
                        && molecule.bond(bond).order.order_int() == 2
                });
                if is_carbonyl {
                    return Some("OCC/3=O");
                }
                let heavy_neighbors = molecule
                    .neighbors(other)
                    .filter(|(neighbor, _)| !atom_is(molecule, *neighbor, "H"))
                    .count();
                if heavy_neighbors == 1 {
                    return Some("OC");
                }
            }
            _ => {}
        }
    }
    None
}

fn is_phosphocholine(molecule: &Molecule, phosphorus: AtomIdx, backbone_oxygen: AtomIdx) -> bool {
    molecule
        .neighbors(phosphorus)
        .map(|(neighbor, _)| neighbor)
        .filter(|neighbor| *neighbor != backbone_oxygen && atom_is(molecule, *neighbor, "O"))
        .any(|bridge_oxygen| {
            let Some(first_carbon) = molecule
                .neighbors(bridge_oxygen)
                .map(|(neighbor, _)| neighbor)
                .find(|neighbor| *neighbor != phosphorus && atom_is(molecule, *neighbor, "C"))
            else {
                return false;
            };
            let Some(second_carbon) = molecule
                .neighbors(first_carbon)
                .map(|(neighbor, _)| neighbor)
                .find(|neighbor| *neighbor != bridge_oxygen && atom_is(molecule, *neighbor, "C"))
            else {
                return false;
            };
            let Some(nitrogen) = molecule
                .neighbors(second_carbon)
                .map(|(neighbor, _)| neighbor)
                .find(|neighbor| *neighbor != first_carbon && atom_is(molecule, *neighbor, "N"))
            else {
                return false;
            };
            molecule
                .neighbors(nitrogen)
                .filter(|(neighbor, _)| {
                    *neighbor != second_carbon && atom_is(molecule, *neighbor, "C")
                })
                .count()
                == 3
        })
}

fn backbone_atoms(ring: &DetectedSugarRing) -> Vec<AtomIdx> {
    let mut atoms = Vec::with_capacity(1 + ring.carbons.len() + ring.tail.len());
    if let Some(carboxyl) = ring.carboxyl_c1 {
        atoms.push(carboxyl);
    } else if let Some(c1) = ring.ketose_c1 {
        atoms.push(c1);
    }
    atoms.extend(ring.carbons.iter().copied());
    atoms.extend(ring.tail.iter().copied());
    atoms
}

fn stereo_for_ligand_order(
    molecule: &Molecule,
    center: AtomIdx,
    desired: &[u32; 4],
) -> Option<char> {
    let atom = molecule.atom(center);
    if atom.chirality == Chirality::None {
        return None;
    }
    let original = molecule.stereo_neighbor_order(center)?;
    if original.len() != 4 {
        return None;
    }
    let permutation = original
        .iter()
        .map(|ligand| desired.iter().position(|candidate| candidate == ligand))
        .collect::<Option<Vec<_>>>()?;
    let inversions = (0..4)
        .flat_map(|left| (left + 1..4).map(move |right| (left, right)))
        .filter(|(left, right)| permutation[*left] > permutation[*right])
        .count();
    let clockwise = atom.chirality == Chirality::Clockwise;
    Some(if clockwise ^ (inversions % 2 == 1) {
        '1'
    } else {
        '2'
    })
}

fn exocyclic_c6(molecule: &Molecule, ring: &DetectedSugarRing) -> Option<AtomIdx> {
    let _ = molecule;
    (!ring.is_ulosonic())
        .then(|| ring.tail.last().copied())
        .flatten()
}

impl DetectedSugarRing {
    fn is_ulosonic(&self) -> bool {
        self.carboxyl_c1.is_some()
    }

    fn is_ketose(&self) -> bool {
        self.ketose_c1.is_some()
    }

    fn anomeric_position(&self) -> u8 {
        if self.is_ulosonic() || self.is_ketose() {
            2
        } else {
            1
        }
    }

    fn ring_end(&self) -> u8 {
        self.anomeric_position() + self.carbons.len() as u8 - 1
    }

    fn ring_closure(&self) -> RingClosure {
        if self.carbons.len() == 4 {
            RingClosure::Furanose
        } else {
            RingClosure::Pyranose
        }
    }
}

fn chain_stereo_descriptor(molecule: &Molecule, backbone: &[AtomIdx], position: usize) -> char {
    if position == 0 || position + 1 >= backbone.len() {
        return 'x';
    }
    let center = backbone[position];
    let previous = backbone[position - 1];
    let next = backbone[position + 1];
    let modifications = molecule
        .neighbors(center)
        .map(|(neighbor, _)| neighbor)
        .filter(|neighbor| *neighbor != previous && *neighbor != next)
        .collect::<Vec<_>>();
    if modifications.len() != 1 {
        return 'x';
    }
    stereo_for_ligand_order(
        molecule,
        center,
        &[previous.0, modifications[0].0, next.0, STEREO_H_SENTINEL],
    )
    .unwrap_or('x')
}

fn anomeric_stereo_descriptor(molecule: &Molecule, ring: &DetectedSugarRing) -> Option<char> {
    let center = ring.carbons[0];
    let next = ring.carbons[1];
    let previous = ring.carboxyl_c1.or(ring.ketose_c1);
    let modifications = molecule
        .neighbors(center)
        .map(|(neighbor, _)| neighbor)
        .filter(|neighbor| *neighbor != next && Some(*neighbor) != previous)
        .collect::<Vec<_>>();
    if modifications.len() != 2 {
        return None;
    }
    // MolWURCS orders the carbon-unit modification ligands independently of
    // SMILES traversal. Anchor that order on the exocyclic anomeric oxygen;
    // the other ligand is the ring oxygen. Global CIP order is unsuitable
    // here because it can flip when a remote glycosyl substituent changes.
    let first = ring.anomeric_oxygen;
    let second = modifications
        .iter()
        .copied()
        .find(|candidate| *candidate != first)?;
    let desired = [
        previous.map_or(STEREO_H_SENTINEL, |atom| atom.0),
        first.0,
        next.0,
        second.0,
    ];
    stereo_for_ligand_order(molecule, center, &desired)
}

fn terminal_descriptor(molecule: &Molecule, ring: &DetectedSugarRing) -> char {
    if ring.is_ulosonic() {
        return 'h';
    }
    let Some(c6) = exocyclic_c6(molecule, ring) else {
        return 'h';
    };
    let oxygen_bonds = molecule
        .neighbors(c6)
        .filter(|(neighbor, _)| atom_is(molecule, *neighbor, "O"))
        .map(|(_, bond)| molecule.bond(bond).order.order_int())
        .collect::<Vec<_>>();
    if oxygen_bonds.contains(&2) && oxygen_bonds.contains(&1) {
        'A'
    } else if oxygen_bonds.is_empty() {
        'm'
    } else {
        'h'
    }
}

fn de_novo_wurcs_for_molecule(molecule: &Molecule) -> MolResult<ResidueGraph> {
    let rings = detect_sugar_rings(molecule);
    if rings.is_empty() {
        return Err(MolError::NoGlycanFound);
    }

    let mut graph = ResidueGraph::new();
    let all_backbone_atoms = rings.iter().flat_map(backbone_atoms).collect::<Vec<_>>();
    let nodes = rings
        .iter()
        .map(|ring| {
            let backbone = backbone_atoms(ring);
            let (anomeric_prefix, skeleton, reference_stereo) = if ring.is_ulosonic() {
                // Nonulosonic acids have C1 carboxyl, C2 anomeric ketal, C3
                // deoxy, stereocentres C4..C8, and terminal C9 alcohol.
                let stereos = (3..=7)
                    .map(|position| chain_stereo_descriptor(molecule, &backbone, position))
                    .collect::<String>();
                let reference = stereos.chars().nth(3).or_else(|| stereos.chars().last());
                let mut skeleton = stereos;
                skeleton.push('h');
                ("Aad", skeleton, reference)
            } else if ring.is_ketose() {
                // 2-ketoses encode terminal C1 (`h`) and anomeric C2 (`a`)
                // in the prefix; the skeleton begins with stereogenic C3.
                let stereos = (2..backbone.len().saturating_sub(1))
                    .map(|position| chain_stereo_descriptor(molecule, &backbone, position))
                    .collect::<String>();
                let reference = stereos.chars().nth(3).or_else(|| stereos.chars().last());
                let mut skeleton = stereos;
                skeleton.push(terminal_descriptor(molecule, ring));
                ("ha", skeleton, reference)
            } else {
                let stereos = (1..backbone.len().saturating_sub(1))
                    .map(|position| chain_stereo_descriptor(molecule, &backbone, position))
                    .collect::<String>();
                let reference = stereos.chars().nth(3).or_else(|| stereos.chars().last());
                let mut skeleton = stereos;
                skeleton.push(terminal_descriptor(molecule, ring));
                ("a", skeleton, reference)
            };
            let anomeric_stereo = anomeric_stereo_descriptor(molecule, ring);
            let anomeric_symbol = match (anomeric_stereo, reference_stereo) {
                (Some(anomeric), Some(reference)) if anomeric == reference => AnomericSymbol::Beta,
                (Some(_), Some(_)) => AnomericSymbol::Alpha,
                _ => AnomericSymbol::Unknown,
            };
            let mut modifications = Vec::new();
            for (offset, carbon) in backbone.iter().enumerate() {
                if let Some(descriptor) = n_acyl_descriptor(molecule, *carbon) {
                    modifications.push(Modification {
                        position: CarbonPosition(offset as u8 + 1),
                        descriptor: descriptor.into(),
                        probability: None,
                    });
                }
                if let Some(descriptor) = n_sulfate_descriptor(molecule, *carbon) {
                    modifications.push(Modification {
                        position: CarbonPosition(offset as u8 + 1),
                        descriptor: descriptor.into(),
                        probability: None,
                    });
                }
                if let Some(descriptor) =
                    oxygen_substituent_descriptor(molecule, *carbon, &all_backbone_atoms)
                {
                    modifications.push(Modification {
                        position: CarbonPosition(offset as u8 + 1),
                        descriptor: descriptor.into(),
                        probability: None,
                    });
                }
            }
            graph.add_residue(Monosaccharide::new(
                skeleton.chars().filter(char::is_ascii_digit).count() as u8,
                skeleton,
                vec![],
                ring.ring_closure(),
                Some(ring.anomeric_position()),
                Some(ring.ring_end()),
                ring.anomeric_position(),
                anomeric_symbol,
                anomeric_prefix.into(),
                modifications,
            ))
        })
        .collect::<Vec<_>>();

    let mut has_parent = vec![false; rings.len()];
    for (child_index, ring) in rings.iter().enumerate() {
        let other = molecule
            .neighbors(ring.anomeric_oxygen)
            .map(|(neighbor, _)| neighbor)
            .find(|neighbor| *neighbor != ring.carbons[0]);
        let Some(other) = other else {
            continue;
        };
        let Some((parent_index, parent_position)) = residue_and_position_for_atom(&rings, other)
        else {
            continue; // reducing hydroxyl or aglycone
        };
        if parent_index == child_index {
            continue;
        }
        graph.add_linkage(
            nodes[parent_index],
            nodes[child_index],
            Linkage::new(
                CarbonPosition(parent_position),
                CarbonPosition(ring.anomeric_position()),
            ),
        );
        has_parent[child_index] = true;
    }

    let roots = has_parent
        .iter()
        .enumerate()
        .filter_map(|(index, has_parent)| (!has_parent).then_some(index))
        .collect::<Vec<_>>();
    if roots.len() != 1 || graph.edge_count() + 1 != graph.node_count() {
        return Err(MolError::AmbiguousGlycanFound);
    }
    let root = nodes[roots[0]];
    graph.set_root(root);
    if let Some(residue) = graph.residue_mut(root) {
        if residue.anomeric_symbol == AnomericSymbol::Unknown {
            residue.anomeric_position = 0;
            residue.anomeric_prefix = match residue.anomeric_prefix.as_str() {
                "Aad" => "AUd".into(),
                "ha" => "hU".into(),
                _ => "u".into(),
            };
        }
    }
    Ok(graph)
}

fn wurcs_for_molecule(molecule: &Molecule) -> MolResult<ResidueGraph> {
    if let Some(wurcs) = exact_wurcs_for_molecule(molecule) {
        return parse_wurcs(wurcs).map_err(MolError::from);
    }
    if let Some(wurcs) = embedded_corpus_wurcs(molecule)? {
        return parse_wurcs(wurcs).map_err(MolError::from);
    }
    de_novo_wurcs_for_molecule(molecule)
}

fn marked_smiles(metadata: &MolMetadata) -> Option<&str> {
    metadata.comment.strip_prefix(SMILES_MARKER)
}

fn graph_from_molecule_and_metadata(
    molecule: &Molecule,
    metadata: &MolMetadata,
) -> MolResult<ResidueGraph> {
    if let Some(smiles) = marked_smiles(metadata) {
        return wurcs_for_molecule(&parse_smiles(smiles)?);
    }
    wurcs_for_molecule(molecule)
}

fn parse_mol_record(input: &str) -> MolResult<(Molecule, MolMetadata)> {
    let result = if input
        .lines()
        .nth(3)
        .is_some_and(|line| line.contains("V3000"))
    {
        chematic::mol::parse_mol_v3000(input)
    } else {
        chematic::mol::parse_mol(input)
    };
    result.map_err(|error| MolError::ParseError(error.to_string()))
}

fn split_sdf_records(input: &str) -> Vec<&str> {
    input
        .split("$$$$")
        .map(str::trim)
        .filter(|record| !record.is_empty())
        .collect()
}

/// Parse a chemical structure and identify it against the bundled glycan
/// molecular corpus. Equivalent non-canonical SMILES spellings are accepted.
pub fn wurcs_from_molecule(input: &str, format: ChemFormat) -> MolResult<ResidueGraph> {
    match format {
        ChemFormat::Smiles => wurcs_for_molecule(&parse_smiles(input)?),
        ChemFormat::Mol => {
            let (molecule, metadata) = parse_mol_record(input)?;
            graph_from_molecule_and_metadata(&molecule, &metadata)
        }
        ChemFormat::Sdf => {
            let records = split_sdf_records(input);
            if records.len() != 1 {
                return Err(MolError::SdfRecordCount(records.len()));
            }
            let (molecule, metadata) = parse_mol_record(records[0])?;
            graph_from_molecule_and_metadata(&molecule, &metadata)
        }
    }
}

/// Parse every chemical record in an input stream.
///
/// MOL and SMILES contain one record. SDF may contain any positive number of
/// records and preserves their input order. Use [`wurcs_from_molecule`] when a
/// caller deliberately requires exactly one result.
pub fn wurcs_from_molecules(input: &str, format: ChemFormat) -> MolResult<Vec<ResidueGraph>> {
    if format != ChemFormat::Sdf {
        return wurcs_from_molecule(input, format).map(|graph| vec![graph]);
    }
    let records = split_sdf_records(input);
    if records.is_empty() {
        return Err(MolError::SdfRecordCount(0));
    }
    records
        .into_iter()
        .map(|record| {
            let (molecule, metadata) = parse_mol_record(record)?;
            graph_from_molecule_and_metadata(&molecule, &metadata)
        })
        .collect()
}

/// Serialize a WURCS graph to a chemical structure.
///
/// MOL/SDF output embeds canonical SMILES in the comment header so this crate
/// can recover complete atom-centred stereochemistry on a later read.
pub fn molecule_from_wurcs(graph: &ResidueGraph, format: ChemFormat) -> MolResult<String> {
    let wurcs = write_wurcs(graph)?;
    let molecule = if let Some(smiles) = corpus_smiles_for_wurcs(&wurcs) {
        parse_smiles(smiles)?
    } else {
        construct::construct_molecule(graph)?
    };
    let smiles = canonical_smiles(&molecule);
    if format == ChemFormat::Smiles {
        return Ok(smiles);
    }

    let canonical = canonical_smiles(&molecule);
    let metadata = MolMetadata::default()
        .with_name("crabWURCS glycan")
        .with_comment(&format!("{SMILES_MARKER}{canonical}"));
    let mol = chematic::mol::write_mol(&molecule, &metadata);
    Ok(match format {
        ChemFormat::Mol => mol,
        ChemFormat::Sdf => format!("{mol}$$$$\n"),
        ChemFormat::Smiles => unreachable!(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use petgraph::visit::EdgeRef;
    use petgraph::Direction;

    const WURCS: &str = "WURCS=2.0/4,4,3/[u2112h_2*NCC/3=O][a2112h-1b_1-5][a2112h-1a_1-5][a1221m-1a_1-5]/1-2-3-4/a3-b1_b3-c1_c2-d1";

    fn canonical_subtree(graph: &ResidueGraph, node: petgraph::graph::NodeIndex) -> String {
        let mut children = graph
            .inner()
            .edges_directed(node, Direction::Outgoing)
            .map(|edge| {
                format!(
                    "{:?}{}",
                    edge.weight(),
                    canonical_subtree(graph, edge.target())
                )
            })
            .collect::<Vec<_>>();
        children.sort();
        let residue = graph.residue(node).expect("graph node has a residue");
        let mut modifications = residue.modifications.clone();
        modifications.sort_by(|left, right| {
            (left.position.0, &left.descriptor).cmp(&(right.position.0, &right.descriptor))
        });
        // Reducing-end WURCS commonly omits the ring declaration even when
        // the molecular input is explicitly cyclic. Compare the chemical
        // residue identity while retaining declared ring data for donors.
        let ring = (residue.anomeric_position > 0).then_some((
            residue.ring,
            residue.ring_start,
            residue.ring_end,
        ));
        format!(
            "({:?}|{}|{:?}|{}{:?}|{:?}[{}])",
            residue.anomeric_prefix,
            residue.skeleton_code,
            ring,
            residue.anomeric_position,
            residue.anomeric_symbol,
            modifications,
            children.join(",")
        )
    }

    fn semantically_equal(left: &ResidueGraph, right: &ResidueGraph) -> bool {
        if left.node_count() != right.node_count() || left.edge_count() != right.edge_count() {
            return false;
        }
        match (left.root(), right.root()) {
            (Some(left_root), Some(right_root)) => {
                canonical_subtree(left, left_root) == canonical_subtree(right, right_root)
            }
            (None, None) => true,
            _ => false,
        }
    }

    fn assert_de_novo_semantic(expected: &str) {
        let (_, smiles) = corpus_pairs()
            .find(|(wurcs, _)| *wurcs == expected)
            .expect("test WURCS exists in the molecular corpus");
        let actual = de_novo_wurcs_for_molecule(&parse_smiles(smiles).unwrap()).unwrap();
        assert!(semantically_equal(&actual, &parse_wurcs(expected).unwrap()));
    }

    #[test]
    fn unknown_valid_smiles_is_not_a_disabled_backend_error() {
        assert!(matches!(
            wurcs_from_molecule("C1CC1", ChemFormat::Smiles),
            Err(MolError::NoGlycanFound)
        ));
    }

    #[test]
    fn unseen_stereounspecified_glucose_gets_conservative_de_novo_wurcs() {
        let graph = wurcs_from_molecule("OCC1OC(O)C(O)C(O)C1O", ChemFormat::Smiles).unwrap();
        assert_eq!(graph.node_count(), 1);
        assert_eq!(write_wurcs(&graph).unwrap(), "WURCS=2.0/1,1,0/[uxxxxh]/1/");
    }

    #[test]
    fn unseen_stereounspecified_disaccharide_preserves_topology() {
        let graph =
            wurcs_from_molecule("OC1OC(O)C(O)C(O)C1OC2OC(O)C(O)C(O)C2O", ChemFormat::Smiles)
                .unwrap();
        assert_eq!(graph.node_count(), 2);
        assert_eq!(graph.edge_count(), 1);
        assert!(write_wurcs(&graph).unwrap().contains("/a2-b1"));
    }

    #[test]
    fn canonical_index_has_no_conflicting_structures() {
        let mut seen = HashMap::new();
        let mut records = 0usize;
        for (wurcs, smiles) in corpus_pairs() {
            records += 1;
            let key = canonical_smiles(&parse_smiles(smiles).unwrap());
            if let Some(previous) = seen.insert(key, wurcs) {
                assert_eq!(previous, wurcs, "one molecular graph mapped to two WURCS");
            }
        }
        assert_eq!(records, 938);
        assert_eq!(seen.len(), canonical_corpus().len());
    }

    #[test]
    fn equivalent_smiles_atom_order_reaches_same_wurcs() {
        let (_, smiles) = corpus_pairs().find(|(wurcs, _)| *wurcs == WURCS).unwrap();
        let molecule = parse_smiles(smiles).unwrap();
        let reordered = chematic::smiles::write(&molecule);
        assert_eq!(
            write_wurcs(&wurcs_from_molecule(&reordered, ChemFormat::Smiles).unwrap()).unwrap(),
            WURCS
        );
    }

    #[test]
    fn corpus_glycan_is_extracted_from_a_methyl_glycoside() {
        const GLYCAN: &str =
            "O1C(O)[C@H](O)[C@@H](O)[C@H](O[C@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@H]2O)[C@H]1CO";
        const EXPECTED: &str = "WURCS=2.0/2,2,1/[u2122h][a2122h-1a_1-5]/1-2/a4-b1";
        let glycoside = GLYCAN.replacen("O1C(O)", "O1C(OC)", 1);
        let extracted = wurcs_from_molecule(&glycoside, ChemFormat::Smiles).unwrap();
        assert_eq!(write_wurcs(&extracted).unwrap(), EXPECTED);
    }

    #[test]
    fn de_novo_stereo_descriptors_match_known_glucose_disaccharide() {
        const GLYCAN: &str =
            "O1C(O)[C@H](O)[C@@H](O)[C@H](O[C@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@H]2O)[C@H]1CO";
        let molecule = parse_smiles(GLYCAN).unwrap();
        let graph = de_novo_wurcs_for_molecule(&molecule).unwrap();
        assert_eq!(
            write_wurcs(&graph).unwrap(),
            "WURCS=2.0/2,2,1/[u2122h][a2122h-1a_1-5]/1-2/a4-b1"
        );
    }

    #[test]
    fn de_novo_handles_nonulosonic_furanose_and_phosphocholine_backbones() {
        for expected in [
            "WURCS=2.0/3,3,2/[u2112h_2*NCC/3=O][a2112h-1a_1-5_2*NCC/3=O][Aad21122h-2a_2-6_5*NCC/3=O]/1-2-3/a3-b1_a6-c2",
            "WURCS=2.0/2,2,1/[u211h][a211h-1b_1-4]/1-2/a2-b1",
            "WURCS=2.0/7,9,8/[u2122h_2*NCC/3=O][a2122h-1b_1-5_2*NCC/3=O][a1122h-1b_1-5][a1122h-1a_1-5][a2122h-1b_1-5_2*NCC/3=O_6*OP^XOCCNC/7C/7C/3O/3=O][a1221m-1a_1-5][a2112h-1b_1-5_2*NCC/3=O]/1-2-3-4-5-6-7-4-6/a4-b1_a6-i1_b4-c1_c3-d1_c6-h1_d2-e1_e3-f1_e4-g1",
        ] {
            assert_de_novo_semantic(expected);
        }
    }

    #[test]
    fn constructed_wurcs_molecules_roundtrip_without_corpus_lookup() {
        for wurcs in [
            "WURCS=2.0/1,1,0/[a2122h-1a_1-5]/1/",
            "WURCS=2.0/2,2,1/[u2122h][a2112h-1b_1-5_4*OSO/3=O/3=O]/1-2/a4-b1",
            "WURCS=2.0/2,2,1/[u211h][Aad21122h-2a_2-6_5*NCC/3=O]/1-2/a2-b2",
        ] {
            let expected = parse_wurcs(wurcs).unwrap();
            let molecule = construct::construct_molecule(&expected).unwrap();
            let actual = de_novo_wurcs_for_molecule(&molecule).unwrap();
            assert!(
                semantically_equal(&actual, &expected),
                "constructed molecule changed {wurcs} into {}",
                write_wurcs(&actual).unwrap()
            );
        }
    }

    #[test]
    fn public_molecule_writer_constructs_an_edited_wurcs_graph() {
        let expected = parse_wurcs("WURCS=2.0/1,1,0/[a2122h-1a_1-5_3*OC]/1/").unwrap();
        let smiles = molecule_from_wurcs(&expected, ChemFormat::Smiles).unwrap();
        let recovered = wurcs_from_molecule(&smiles, ChemFormat::Smiles).unwrap();
        assert!(
            semantically_equal(&expected, &recovered),
            "{smiles} -> {}",
            write_wurcs(&recovered).unwrap()
        );
    }

    #[test]
    fn corpus_graph_roundtrips_through_mol_and_sdf() {
        let graph = parse_wurcs(WURCS).unwrap();
        for format in [ChemFormat::Mol, ChemFormat::Sdf] {
            let chemical = molecule_from_wurcs(&graph, format).unwrap();
            assert!(chemical.contains(SMILES_MARKER));
            let recovered = wurcs_from_molecule(&chemical, format).unwrap();
            assert_eq!(write_wurcs(&recovered).unwrap(), WURCS);
        }
    }

    #[test]
    fn multi_record_sdf_returns_every_glycan_in_order() {
        let first = parse_wurcs(WURCS).unwrap();
        let second_wurcs = "WURCS=2.0/1,1,0/[u2122h_2*NCC/3=O]/1/";
        let second = parse_wurcs(second_wurcs).unwrap();
        let sdf = format!(
            "{}{}",
            molecule_from_wurcs(&first, ChemFormat::Sdf).unwrap(),
            molecule_from_wurcs(&second, ChemFormat::Sdf).unwrap()
        );
        let recovered = wurcs_from_molecules(&sdf, ChemFormat::Sdf).unwrap();
        assert_eq!(recovered.len(), 2);
        assert_eq!(write_wurcs(&recovered[0]).unwrap(), WURCS);
        assert_eq!(write_wurcs(&recovered[1]).unwrap(), second_wurcs);
    }

    #[test]
    fn malformed_and_multi_record_inputs_are_explicit_errors() {
        assert!(matches!(
            wurcs_from_molecule("not a mol", ChemFormat::Mol),
            Err(MolError::ParseError(_))
        ));
        assert!(matches!(
            wurcs_from_molecule("$$$$\n$$$$\n", ChemFormat::Sdf),
            Err(MolError::SdfRecordCount(0))
        ));
    }

    #[test]
    #[ignore = "full 938-structure de-novo coverage audit"]
    fn de_novo_corpus_audit() {
        let mut extracted = 0usize;
        let mut exact = 0usize;
        let mut semantic = 0usize;
        let mut mismatches = 0usize;
        for (expected, smiles) in corpus_pairs() {
            let molecule = parse_smiles(smiles).unwrap();
            let Ok(graph) = de_novo_wurcs_for_molecule(&molecule) else {
                continue;
            };
            extracted += 1;
            let actual = write_wurcs(&graph).unwrap();
            if actual == expected {
                exact += 1;
            }
            let expected_graph = parse_wurcs(expected).unwrap();
            if semantically_equal(&graph, &expected_graph) {
                semantic += 1;
            } else if mismatches < 10 {
                eprintln!("expected: {expected}\nactual:   {actual}\n");
                mismatches += 1;
            }
        }
        eprintln!(
            "de-novo corpus audit: extracted={extracted}/938 semantic={semantic}/938 exact={exact}/938"
        );
        assert_eq!(extracted, 938);
        assert_eq!(semantic, 938);
        assert!(exact >= 226);
    }

    #[test]
    #[ignore = "full 938-structure de-novo construction audit"]
    fn constructed_corpus_audit() {
        let mut constructed = 0usize;
        let mut semantic = 0usize;
        let mut failures = 0usize;
        for (expected, _) in corpus_pairs() {
            let graph = parse_wurcs(expected).unwrap();
            let molecule = match construct::construct_molecule(&graph) {
                Ok(molecule) => molecule,
                Err(error) => {
                    if failures < 10 {
                        eprintln!("construction failed: {expected}\n{error}\n");
                    }
                    failures += 1;
                    continue;
                }
            };
            constructed += 1;
            let serialized = canonical_smiles(&molecule);
            let reparsed = parse_smiles(&serialized).unwrap();
            let Ok(recovered) = de_novo_wurcs_for_molecule(&reparsed) else {
                if failures < 10 {
                    eprintln!("re-extraction failed: {expected}\n");
                }
                failures += 1;
                continue;
            };
            if semantically_equal(&graph, &recovered) {
                semantic += 1;
            } else {
                if failures < 10 {
                    eprintln!(
                        "construction mismatch: {expected}\nactual: {}\n",
                        write_wurcs(&recovered).unwrap()
                    );
                }
                failures += 1;
            }
        }
        eprintln!(
            "construction corpus audit: constructed={constructed}/938 semantic={semantic}/938"
        );
        assert_eq!(constructed, 938);
        assert_eq!(semantic, 938);
    }
}
