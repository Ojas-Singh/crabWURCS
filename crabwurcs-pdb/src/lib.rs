use crabwurcs_core::{
    AnomericSymbol, CarbonPosition, Linkage, Modification, Monosaccharide, ResidueGraph,
    ResidueKind, RingClosure, residue_from_kind,
};
use pdbtbx::{
    ContainsAtomConformer, ContainsAtomConformerResidue, ContainsAtomConformerResidueChain,
};
use std::collections::{HashMap, HashSet};
use thiserror::Error;

#[derive(Debug, Error)]
pub enum PdbError {
    #[error("failed to parse structure file: {0}")]
    ParseError(String),

    #[error("no carbohydrate residues found in structure")]
    NoGlycansFound,

    #[error(transparent)]
    Core(#[from] crabwurcs_core::CoreError),
}

pub type PdbResult<T> = Result<T, PdbError>;

/// A graph-node provenance record retained from the source PDB/mmCIF file.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PdbResidueReference {
    pub node_index: usize,
    pub chain: String,
    pub sequence_number: isize,
    pub insertion_code: Option<String>,
}

#[derive(Debug, Clone)]
pub struct ExtractedGlycan {
    pub attachment_site: Option<String>,
    pub graph: ResidueGraph,
}

/// An extracted glycan together with the source residue identity for every
/// graph node. The legacy [`ExtractedGlycan`] API remains available for
/// callers that only need the graph.
#[derive(Debug, Clone)]
pub struct ExtractedGlycanWithProvenance {
    pub attachment_site: Option<String>,
    pub graph: ResidueGraph,
    pub residues: Vec<PdbResidueReference>,
}

const PDB_COMPONENTS: &str = include_str!("../data/pdb_carbohydrate_components.tsv");

fn is_sugar_residue(name: &str) -> bool {
    pdb_component_residue(name).is_some()
        || glycam_residue_kind(name).is_some()
        || registry_named_residue(name).is_some()
}

fn registry_named_residue(name: &str) -> Option<Monosaccharide> {
    let kind = ResidueKind::from_name(name)?;
    (!kind.is_generic())
        .then(|| residue_from_kind(kind).ok())
        .flatten()
}

fn apply_declared_form(
    mut residue: Monosaccharide,
    anomer: AnomericSymbol,
    ring: Option<RingClosure>,
) -> Monosaccharide {
    residue.anomeric_symbol = anomer;
    if let Some(ring) = ring {
        residue.ring = ring;
        let start = residue
            .ring_start
            .unwrap_or(residue.anomeric_position.max(1));
        residue.ring_start = Some(start);
        residue.ring_end = Some(match ring {
            RingClosure::Furanose => start.saturating_add(3),
            RingClosure::Pyranose => start.saturating_add(4),
            RingClosure::Open => start,
            RingClosure::Unknown => residue.ring_end.unwrap_or(start),
        });
    }
    residue
}

fn pdb_component_residue(name: &str) -> Option<Monosaccharide> {
    let fields = PDB_COMPONENTS
        .lines()
        .filter(|line| !line.starts_with('#') && !line.trim().is_empty())
        .find_map(|line| {
            let mut fields = line.split('\t');
            (fields.next()? == name).then(|| {
                (
                    fields.next().unwrap_or_default(),
                    fields.next().unwrap_or("x"),
                    fields.next().unwrap_or("x"),
                    fields.next().unwrap_or("n"),
                )
            })
        })?;
    let kind = ResidueKind::from_name(fields.0)?;
    let mut residue = residue_from_kind(kind).ok()?;
    if fields.3 == "y" {
        mirror_stereochemistry(&mut residue);
    }
    let anomer = match fields.1 {
        "a" => AnomericSymbol::Alpha,
        "b" => AnomericSymbol::Beta,
        _ => AnomericSymbol::Unknown,
    };
    let ring = match fields.2 {
        "p" => Some(RingClosure::Pyranose),
        "f" => Some(RingClosure::Furanose),
        _ => None,
    };
    Some(apply_declared_form(residue, anomer, ring))
}

fn mirror_stereochemistry(residue: &mut Monosaccharide) {
    residue.skeleton_code = residue
        .skeleton_code
        .chars()
        .map(|descriptor| match descriptor {
            '1' => '2',
            '2' => '1',
            '3' => '4',
            '4' => '3',
            '5' => '6',
            '6' => '5',
            '7' => '8',
            '8' => '7',
            other => other,
        })
        .collect();
}

fn glycam_residue_kind(name: &str) -> Option<Monosaccharide> {
    let bytes = name.as_bytes();
    if bytes.len() != 3 || !b"0123456789ZYXWVUTSRQP".contains(&bytes[0]) {
        return None;
    }
    let middle = bytes[1] as char;
    let last = bytes[2] as char;
    let kind = match (middle, last) {
        ('G' | 'g', 'L' | 'l') => ResidueKind::Neu5Gc,
        ('Y' | 'y', 'N' | 'n' | 'S' | 's') => ResidueKind::GlcN,
        ('K' | 'k', 'N' | 'n') => ResidueKind::Kdn,
        ('K' | 'k', 'O' | 'o') => ResidueKind::Kdo,
        ('B' | 'b', 'C' | 'c') => ResidueKind::Bac,
        ('Z' | 'z', _) => ResidueKind::GlcA,
        ('O' | 'o', _) => ResidueKind::GalA,
        ('U' | 'u', _) => ResidueKind::IdoA,
        ('C' | 'c', _) => ResidueKind::Fru,
        ('A' | 'a', _) => ResidueKind::Ara,
        ('H' | 'h', _) => ResidueKind::Rha,
        ('Q' | 'q', _) => ResidueKind::Qui,
        ('W' | 'w', _) => ResidueKind::ManNAc,
        ('K' | 'k', _) => ResidueKind::Gul,
        ('E' | 'e', _) => ResidueKind::Alt,
        ('N' | 'n', _) => ResidueKind::All,
        ('T' | 't', _) => ResidueKind::Tal,
        ('P' | 'p', 'D' | 'U') => ResidueKind::Psi,
        ('R' | 'r', 'D' | 'U') => ResidueKind::Rib,
        ('D' | 'd', 'A' | 'B') => ResidueKind::Lyx,
        ('J' | 'j', 'A' | 'B') => ResidueKind::Tag,
        ('b', 'A' | 'B') => ResidueKind::Sor,
        _ => match middle {
            'G' | 'g' => ResidueKind::Glc,
            'Y' | 'y' => ResidueKind::GlcNAc,
            'L' | 'l' => ResidueKind::Gal,
            'V' | 'v' => ResidueKind::GalNAc,
            'M' | 'm' => ResidueKind::Man,
            'F' | 'f' => ResidueKind::Fuc,
            'X' | 'x' => ResidueKind::Xyl,
            'S' | 's' => ResidueKind::Neu5Ac,
            _ => return None,
        },
    };
    let anomer = match last {
        'A' | 'D' => AnomericSymbol::Alpha,
        'B' | 'U' => AnomericSymbol::Beta,
        'N' | 'S' | 'L' => AnomericSymbol::Alpha,
        'n' | 's' | 'l' => AnomericSymbol::Beta,
        _ if matches!(kind, ResidueKind::Bac | ResidueKind::Kdo | ResidueKind::Kdn) => {
            AnomericSymbol::Alpha
        }
        _ => return None,
    };
    let mut residue = residue_from_kind(kind).ok()?;
    let registry_template_is_lowercase_glycam_form =
        matches!(kind, ResidueKind::Ara | ResidueKind::Fuc | ResidueKind::Rha);
    if (registry_template_is_lowercase_glycam_form && middle.is_ascii_uppercase())
        || (!registry_template_is_lowercase_glycam_form
            && middle.is_ascii_lowercase()
            && !matches!(kind, ResidueKind::Neu5Ac | ResidueKind::Neu5Gc))
    {
        mirror_stereochemistry(&mut residue);
    }
    let ring = matches!(last, 'D' | 'U').then_some(RingClosure::Furanose);
    residue = apply_declared_form(residue, anomer, ring);
    if matches!(last, 'S' | 's')
        && let Some(modification) = residue
            .modifications
            .iter_mut()
            .find(|modification| modification.descriptor == "N")
    {
        modification.descriptor = "NSO/3=O/3=O".into();
    }
    Some(residue)
}

fn modification(position: u8, descriptor: &str) -> Modification {
    Modification {
        position: CarbonPosition(position),
        descriptor: descriptor.to_string(),
        probability: None,
    }
}

fn residue_to_monosaccharide(name: &str) -> Monosaccharide {
    pdb_component_residue(name)
        .or_else(|| registry_named_residue(name))
        .or_else(|| glycam_residue_kind(name))
        .unwrap_or_else(|| {
            residue_from_kind(ResidueKind::Hex).expect("the generic hexose registry entry is valid")
        })
}

pub fn extract_glycans_from_file(path: &std::path::Path) -> PdbResult<Vec<ExtractedGlycan>> {
    use pdbtbx::ReadOptions;

    let contents = std::fs::read_to_string(path)
        .map_err(|error| PdbError::ParseError(format!("cannot read file: {error}")))?;
    let path_str = path.to_string_lossy();
    let (pdb, _errors) = ReadOptions::default()
        .set_level(pdbtbx::StrictnessLevel::Loose)
        .read(path_str.as_ref())
        .map_err(|e| PdbError::ParseError(format!("cannot parse file: {:?}", e)))?;

    Ok(extract_glycans_from_pdb_with_provenance(
        &pdb,
        &raw_pdb_residue_names(&contents),
        &raw_pdb_bonds(&contents),
    )?
    .into_iter()
    .map(|glycan| ExtractedGlycan {
        attachment_site: glycan.attachment_site,
        graph: glycan.graph,
    })
    .collect())
}

pub fn extract_glycans_from_str(contents: &str, is_mmcif: bool) -> PdbResult<Vec<ExtractedGlycan>> {
    use pdbtbx::{Format, ReadOptions, StrictnessLevel};
    use std::io::BufReader;

    let format = if is_mmcif { Format::Mmcif } else { Format::Pdb };
    let reader = BufReader::new(contents.as_bytes());
    let (pdb, _errors) = ReadOptions::default()
        .set_level(StrictnessLevel::Loose)
        .set_format(format)
        .read_raw(reader)
        .map_err(|e| PdbError::ParseError(format!("cannot parse: {:?}", e)))?;

    Ok(extract_glycans_from_pdb_with_provenance(
        &pdb,
        &raw_pdb_residue_names(contents),
        &raw_pdb_bonds(contents),
    )?
    .into_iter()
    .map(|glycan| ExtractedGlycan {
        attachment_site: glycan.attachment_site,
        graph: glycan.graph,
    })
    .collect())
}

/// Extract glycans from a PDB/mmCIF file while retaining source residue
/// identities for every graph node.
pub fn extract_glycans_with_provenance_from_file(
    path: &std::path::Path,
) -> PdbResult<Vec<ExtractedGlycanWithProvenance>> {
    use pdbtbx::ReadOptions;

    let contents = std::fs::read_to_string(path)
        .map_err(|error| PdbError::ParseError(format!("cannot read file: {error}")))?;
    let path_str = path.to_string_lossy();
    let (pdb, _errors) = ReadOptions::default()
        .set_level(pdbtbx::StrictnessLevel::Loose)
        .read(path_str.as_ref())
        .map_err(|e| PdbError::ParseError(format!("cannot parse file: {:?}", e)))?;
    extract_glycans_from_pdb_with_provenance(
        &pdb,
        &raw_pdb_residue_names(&contents),
        &raw_pdb_bonds(&contents),
    )
}

/// Extract glycans from a PDB or mmCIF string while retaining source residue
/// identities for every graph node.
pub fn extract_glycans_with_provenance_from_str(
    contents: &str,
    is_mmcif: bool,
) -> PdbResult<Vec<ExtractedGlycanWithProvenance>> {
    use pdbtbx::{Format, ReadOptions, StrictnessLevel};
    use std::io::BufReader;

    let format = if is_mmcif { Format::Mmcif } else { Format::Pdb };
    let reader = BufReader::new(contents.as_bytes());
    let (pdb, _errors) = ReadOptions::default()
        .set_level(StrictnessLevel::Loose)
        .set_format(format)
        .read_raw(reader)
        .map_err(|e| PdbError::ParseError(format!("cannot parse: {:?}", e)))?;
    extract_glycans_from_pdb_with_provenance(
        &pdb,
        &raw_pdb_residue_names(contents),
        &raw_pdb_bonds(contents),
    )
}

fn raw_pdb_bonds(contents: &str) -> Vec<(usize, usize)> {
    let mut bonds = HashSet::new();
    for line in contents.lines().filter(|line| line.starts_with("CONECT")) {
        let serials = line
            .as_bytes()
            .get(6..)
            .into_iter()
            .flat_map(|bytes| bytes.chunks(5))
            .filter_map(|field| {
                std::str::from_utf8(field)
                    .ok()?
                    .trim()
                    .parse::<usize>()
                    .ok()
            })
            .collect::<Vec<_>>();
        let Some(&first) = serials.first() else {
            continue;
        };
        for &second in &serials[1..] {
            bonds.insert(if first < second {
                (first, second)
            } else {
                (second, first)
            });
        }
    }
    bonds.into_iter().collect()
}

fn raw_pdb_residue_names(contents: &str) -> HashMap<(isize, String), String> {
    let mut names = HashMap::new();
    for line in contents.lines() {
        if !(line.starts_with("ATOM  ") || line.starts_with("HETATM")) || line.len() < 26 {
            continue;
        }
        let raw = line[17..20].trim();
        let Ok(sequence) = line[22..26].trim().parse::<isize>() else {
            continue;
        };
        names
            .entry((sequence, raw.to_ascii_uppercase()))
            .or_insert_with(|| raw.to_string());
    }
    names
}

fn extract_glycans_from_pdb_with_provenance(
    pdb: &pdbtbx::PDB,
    raw_names: &HashMap<(isize, String), String>,
    raw_bonds: &[(usize, usize)],
) -> PdbResult<Vec<ExtractedGlycanWithProvenance>> {
    #[derive(Debug, Clone)]
    struct ResidueMeta {
        name: String,
        sequence: isize,
        insertion_code: Option<String>,
        chain: String,
        order: usize,
    }

    #[derive(Debug, Clone)]
    struct AtomMeta {
        residue: usize,
        name: String,
        element: Option<String>,
        charge: i8,
        position: [f64; 3],
    }

    let mut residues: HashMap<usize, ResidueMeta> = HashMap::new();
    let mut atoms: HashMap<usize, AtomMeta> = HashMap::new();
    let mut atoms_by_serial = HashMap::new();
    for hierarchy in pdb.atoms_with_hierarchy() {
        let residue_key = std::ptr::from_ref(hierarchy.residue()) as usize;
        let atom_key = std::ptr::from_ref(hierarchy.atom()) as usize;
        let next_order = residues.len();
        residues.entry(residue_key).or_insert_with(|| {
            let parsed_name = hierarchy.residue().name().unwrap_or("UNK");
            let (sequence, insertion_code) = hierarchy.residue().id();
            ResidueMeta {
                name: raw_names
                    .get(&(sequence, parsed_name.to_string()))
                    .cloned()
                    .unwrap_or_else(|| parsed_name.to_string()),
                sequence,
                insertion_code: insertion_code.map(str::to_string),
                chain: hierarchy.chain().id().to_string(),
                order: next_order,
            }
        });
        atoms.insert(
            atom_key,
            AtomMeta {
                residue: residue_key,
                name: hierarchy.atom().name().to_string(),
                element: hierarchy
                    .atom()
                    .element()
                    .map(|element| element.symbol().to_string()),
                charge: hierarchy.atom().charge().clamp(-8, 8) as i8,
                position: [
                    hierarchy.atom().x(),
                    hierarchy.atom().y(),
                    hierarchy.atom().z(),
                ],
            },
        );
        atoms_by_serial.insert(hierarchy.atom().serial_number(), atom_key);
    }

    let mut sugar_keys = residues
        .iter()
        .filter_map(|(key, residue)| is_sugar_residue(&residue.name).then_some(*key))
        .collect::<HashSet<_>>();
    let mut inferred_residues = HashMap::<usize, Monosaccharide>::new();

    // A component ID is authoritative when it is in the bundled CCD table.
    // For renamed/private components, fall back to the actual atom graph. PDB
    // files do not always carry intra-residue CONECT records, so conservative
    // covalent-radius perception is used only inside one residue.
    for key in residues
        .keys()
        .copied()
        .filter(|key| !sugar_keys.contains(key))
        .collect::<Vec<_>>()
    {
        let residue_atoms = atoms
            .iter()
            .filter(|(_, atom)| atom.residue == key)
            .collect::<Vec<_>>();
        let carbon_or_oxygen = residue_atoms
            .iter()
            .filter(|(_, atom)| matches!(atom.element.as_deref(), Some("C" | "O")))
            .count();
        if residue_atoms.len() < 5 || carbon_or_oxygen < 4 {
            continue;
        }
        let mut input = crabwurcs_mol::MolecularGraphInput::default();
        let mut atom_indices = HashMap::new();
        for (index, (atom_key, atom)) in residue_atoms.iter().enumerate() {
            let Some(element) = atom.element.as_deref() else {
                continue;
            };
            atom_indices.insert(**atom_key, input.atoms.len());
            input.atoms.push(crabwurcs_mol::MolecularAtom {
                element: element.to_string(),
                formal_charge: atom.charge,
                coordinates: Some(atom.position),
                stereo: crabwurcs_mol::InputStereo::Unspecified,
            });
            let _ = index;
        }
        if input.atoms.len() < 5 {
            continue;
        }
        let radius = |symbol: &str| -> f64 {
            match symbol {
                "H" => 0.31,
                "C" => 0.76,
                "N" => 0.71,
                "O" => 0.66,
                "P" => 1.07,
                "S" => 1.05,
                "F" => 0.57,
                "Cl" => 1.02,
                _ => 0.8,
            }
        };
        for left in 0..input.atoms.len() {
            for right in left + 1..input.atoms.len() {
                let left_source = residue_atoms
                    .iter()
                    .find(|(key, _)| atom_indices.get(key) == Some(&left))
                    .map(|(_, atom)| atom);
                let right_source = residue_atoms
                    .iter()
                    .find(|(key, _)| atom_indices.get(key) == Some(&right))
                    .map(|(_, atom)| atom);
                let (Some(left_source), Some(right_source)) = (left_source, right_source) else {
                    continue;
                };
                let distance_squared = left_source
                    .position
                    .iter()
                    .zip(right_source.position)
                    .map(|(a, b)| (a - b).powi(2))
                    .sum::<f64>();
                let maximum =
                    radius(&input.atoms[left].element) + radius(&input.atoms[right].element) + 0.35;
                if !(0.45f64.powi(2)..=maximum.powi(2)).contains(&distance_squared) {
                    continue;
                }
                let distance = distance_squared.sqrt();
                let order = if matches!(
                    (
                        input.atoms[left].element.as_str(),
                        input.atoms[right].element.as_str()
                    ),
                    ("C", "O") | ("O", "C")
                ) && distance < 1.32
                {
                    crabwurcs_mol::MolecularBondOrder::Double
                } else {
                    crabwurcs_mol::MolecularBondOrder::Single
                };
                input.bonds.push(crabwurcs_mol::MolecularBond {
                    atom1: left,
                    atom2: right,
                    order,
                });
            }
        }
        let Ok(graph) = crabwurcs_mol::wurcs_from_atom_graph(&input) else {
            continue;
        };
        if graph.node_count() != 1 {
            continue;
        }
        let Some(residue) = graph.root().and_then(|root| graph.residue(root)).cloned() else {
            continue;
        };
        inferred_residues.insert(key, residue);
        sugar_keys.insert(key);
    }
    if sugar_keys.is_empty() {
        return Ok(vec![]);
    }

    // A glycosidic bond is C(anomeric)–O(acceptor).  CONECT/LINK records are
    // authoritative: residue order in a coordinate file has no topological
    // meaning and must never be used to invent a linear 1→4 chain.
    let mut glycosidic = HashSet::new();
    let mut external_attachments: HashMap<usize, String> = HashMap::new();
    let mut bonded_modifications = HashSet::new();
    {
        let mut record_bond = |first_meta: &AtomMeta, second_meta: &AtomMeta| {
            let first_residue = &first_meta.residue;
            let second_residue = &second_meta.residue;
            if first_residue == second_residue {
                return;
            }

            let first_is_sugar = sugar_keys.contains(first_residue);
            let second_is_sugar = sugar_keys.contains(second_residue);
            if first_is_sugar && second_is_sugar {
                if let Some((parent, child, parent_position, child_position)) =
                    orient_glycosidic_bond(
                        *first_residue,
                        &first_meta.name,
                        *second_residue,
                        &second_meta.name,
                    )
                {
                    glycosidic.insert((parent, child, parent_position, child_position));
                }
            } else if first_is_sugar || second_is_sugar {
                let (sugar, sugar_atom, other) = if first_is_sugar {
                    (*first_residue, first_meta, *second_residue)
                } else {
                    (*second_residue, second_meta, *first_residue)
                };
                if let Some(other) = residues.get(&other) {
                    let descriptor = match other.name.as_str() {
                        "MEX" => Some("OC"),
                        "ACX" => Some("OCC/3=O"),
                        "PCX" => Some("OP^XOCCNC/7C/7C/3O/3=O"),
                        _ => None,
                    };
                    if let (Some(position), Some(descriptor)) =
                        (atom_position(&sugar_atom.name, 'O'), descriptor)
                    {
                        bonded_modifications.insert((sugar, position, descriptor.to_string()));
                    } else {
                        external_attachments.entry(sugar).or_insert_with(|| {
                            format!("{}/{}/{}", other.chain, other.name, other.sequence)
                        });
                    }
                }
            }
        };

        for (first, second, _) in pdb.bonds() {
            let first_key = std::ptr::from_ref(first) as usize;
            let second_key = std::ptr::from_ref(second) as usize;
            let (Some(first_meta), Some(second_meta)) =
                (atoms.get(&first_key), atoms.get(&second_key))
            else {
                continue;
            };
            record_bond(first_meta, second_meta);
        }
        for &(first_serial, second_serial) in raw_bonds {
            if let (Some(first), Some(second)) = (
                atoms_by_serial
                    .get(&first_serial)
                    .and_then(|key| atoms.get(key)),
                atoms_by_serial
                    .get(&second_serial)
                    .and_then(|key| atoms.get(key)),
            ) {
                record_bond(first, second);
            }
        }
    }

    // pdbtbx imports LINK/SSBOND records but not the PDB CONECT records used
    // by GlycoShape. Recover those covalent links from geometry. Restricting
    // the carbon endpoint to the child's declared anomeric carbon prevents
    // ordinary close contacts from becoming false glycosidic linkages.
    let sugar_atoms = atoms
        .values()
        .filter(|atom| sugar_keys.contains(&atom.residue))
        .collect::<Vec<_>>();
    if glycosidic.is_empty() {
        for (index, first) in sugar_atoms.iter().enumerate() {
            for second in sugar_atoms.iter().skip(index + 1) {
                if first.residue == second.residue {
                    continue;
                }
                let Some((parent, child, parent_position, child_position)) = orient_glycosidic_bond(
                    first.residue,
                    &first.name,
                    second.residue,
                    &second.name,
                ) else {
                    continue;
                };
                if child_position
                    != residue_to_monosaccharide(&residues[&child].name).anomeric_position
                {
                    continue;
                }
                let distance_squared = first
                    .position
                    .iter()
                    .zip(second.position)
                    .map(|(a, b)| (a - b).powi(2))
                    .sum::<f64>();
                if (0.8f64.powi(2)..=1.8f64.powi(2)).contains(&distance_squared) {
                    glycosidic.insert((parent, child, parent_position, child_position));
                }
            }
        }
    }

    let sulfate_atoms = atoms
        .values()
        .filter(|atom| {
            residues
                .get(&atom.residue)
                .is_some_and(|residue| residue.name == "SO3")
                && atom.name.starts_with('S')
        })
        .collect::<Vec<_>>();
    let mut sulfate_positions: HashMap<usize, HashSet<u8>> = HashMap::new();
    for sugar_atom in &sugar_atoms {
        let Some(position) = atom_position(&sugar_atom.name, 'O') else {
            continue;
        };
        for sulfate in &sulfate_atoms {
            let distance_squared = sugar_atom
                .position
                .iter()
                .zip(sulfate.position)
                .map(|(a, b)| (a - b).powi(2))
                .sum::<f64>();
            if (0.8f64.powi(2)..=1.9f64.powi(2)).contains(&distance_squared) {
                sulfate_positions
                    .entry(sugar_atom.residue)
                    .or_default()
                    .insert(position);
            }
        }
    }

    let mut adjacency: HashMap<usize, Vec<usize>> = HashMap::new();
    for &(parent, child, _, _) in &glycosidic {
        adjacency.entry(parent).or_default().push(child);
        adjacency.entry(child).or_default().push(parent);
    }

    let mut remaining = sugar_keys.clone();
    let mut components = Vec::new();
    while let Some(start) = remaining.iter().next().copied() {
        let mut stack = vec![start];
        let mut component = HashSet::new();
        while let Some(current) = stack.pop() {
            if !component.insert(current) {
                continue;
            }
            remaining.remove(&current);
            stack.extend(adjacency.get(&current).into_iter().flatten().copied());
        }
        components.push(component);
    }

    components.sort_by_key(|component| {
        component
            .iter()
            .filter_map(|key| residues.get(key).map(|residue| residue.order))
            .min()
            .unwrap_or(usize::MAX)
    });

    let mut extracted = Vec::with_capacity(components.len());
    for component in components {
        let mut ordered = component.iter().copied().collect::<Vec<_>>();
        ordered.sort_by_key(|key| residues[key].order);

        let mut graph = ResidueGraph::new();
        let mut nodes = HashMap::new();
        let mut source_residues = Vec::with_capacity(ordered.len());
        for key in &ordered {
            let mut residue = inferred_residues
                .get(key)
                .cloned()
                .unwrap_or_else(|| residue_to_monosaccharide(&residues[key].name));
            if let Some(positions) = sulfate_positions.get(key) {
                let mut positions = positions.iter().copied().collect::<Vec<_>>();
                positions.sort_unstable();
                residue.modifications.extend(
                    positions
                        .into_iter()
                        .map(|position| modification(position, "OSO/3=O/3=O")),
                );
            }
            residue.modifications.extend(
                bonded_modifications
                    .iter()
                    .filter(|(residue, _, _)| residue == key)
                    .map(|(_, position, descriptor)| modification(*position, descriptor)),
            );
            residue.modifications.sort_by(|first, second| {
                (first.position.0, first.descriptor.as_str())
                    .cmp(&(second.position.0, second.descriptor.as_str()))
            });
            residue.modifications.dedup_by(|first, second| {
                first.position == second.position && first.descriptor == second.descriptor
            });
            let node = graph.add_residue(residue);
            nodes.insert(*key, node);
            source_residues.push(PdbResidueReference {
                node_index: node.index(),
                chain: residues[key].chain.clone(),
                sequence_number: residues[key].sequence,
                insertion_code: residues[key].insertion_code.clone(),
            });
        }
        for &(parent, child, parent_position, child_position) in &glycosidic {
            if component.contains(&parent) && component.contains(&child) {
                graph.add_linkage(
                    nodes[&parent],
                    nodes[&child],
                    Linkage::new(
                        CarbonPosition(parent_position),
                        CarbonPosition(child_position),
                    ),
                );
            }
        }

        let children = glycosidic
            .iter()
            .filter(|(parent, child, _, _)| component.contains(parent) && component.contains(child))
            .map(|(_, child, _, _)| *child)
            .collect::<HashSet<_>>();
        let root = ordered
            .iter()
            .copied()
            .find(|key| !children.contains(key))
            .unwrap_or(ordered[0]);
        graph.set_root(nodes[&root]);

        let attachment_site = component
            .iter()
            .filter_map(|key| external_attachments.get(key))
            .next()
            .cloned()
            .or_else(|| {
                let residue = &residues[&root];
                Some(format!(
                    "{}/{}/{}",
                    residue.chain, residue.name, residue.sequence
                ))
            });
        extracted.push(ExtractedGlycanWithProvenance {
            attachment_site,
            graph,
            residues: source_residues,
        });
    }

    Ok(extracted)
}

fn atom_position(name: &str, element: char) -> Option<u8> {
    let suffix = name.strip_prefix(element)?;
    let digits = suffix
        .chars()
        .take_while(char::is_ascii_digit)
        .collect::<String>();
    (!digits.is_empty()).then(|| digits.parse().ok()).flatten()
}

fn orient_glycosidic_bond(
    first_residue: usize,
    first_atom: &str,
    second_residue: usize,
    second_atom: &str,
) -> Option<(usize, usize, u8, u8)> {
    if let (Some(parent_position), Some(child_position)) = (
        atom_position(first_atom, 'O'),
        atom_position(second_atom, 'C'),
    ) {
        return Some((
            first_residue,
            second_residue,
            parent_position,
            child_position,
        ));
    }
    let parent_position = atom_position(second_atom, 'O')?;
    let child_position = atom_position(first_atom, 'C')?;
    Some((
        second_residue,
        first_residue,
        parent_position,
        child_position,
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[allow(clippy::too_many_arguments)]
    fn make_hetatm_line(
        serial: u32,
        atom_name: &str,
        res_name: &str,
        chain: &str,
        seq: u32,
        x: f64,
        y: f64,
        z: f64,
    ) -> String {
        let element = atom_name
            .chars()
            .find(char::is_ascii_alphabetic)
            .unwrap_or('C')
            .to_ascii_uppercase()
            .to_string();
        format!(
            "HETATM{:5} {:<4}{:1}{:3} {}{:4}{:1}   {:8.3}{:8.3}{:8.3}{:6.2}{:6.2}          {:>2}{:2}",
            serial, atom_name, "", res_name, chain, seq, "", x, y, z, 1.0f64, 0.0f64, element, ""
        )
    }

    #[test]
    fn test_sugar_residue_detection() {
        assert!(is_sugar_residue("NAG"));
        assert!(is_sugar_residue("MAN"));
        assert!(is_sugar_residue("FUC"));
        assert!(!is_sugar_residue("HOH"));
        assert!(!is_sugar_residue("ALA"));
    }

    #[test]
    fn test_residue_to_monosaccharide() {
        let nag = residue_to_monosaccharide("NAG");
        assert_eq!(nag.backbone_length, 4);

        let sia = residue_to_monosaccharide("SIA");
        assert_eq!(sia.backbone_length, 5);
    }

    #[test]
    fn glycam_and_pdb_special_residues_keep_their_chemistry() {
        let fructose = residue_to_monosaccharide("0CU");
        assert_eq!(fructose.ring, RingClosure::Furanose);
        assert_eq!(fructose.anomeric_position, 2);
        assert_eq!(fructose.anomeric_symbol, AnomericSymbol::Beta);

        for code in ["0GL", "0gL"] {
            let neu5gc = residue_to_monosaccharide(code);
            assert_eq!(neu5gc.anomeric_position, 2);
            assert!(
                neu5gc
                    .modifications
                    .iter()
                    .any(|modification| modification.position.0 == 5
                        && modification.descriptor == "NCCO/3=O")
            );
        }

        let arabinofuranose = residue_to_monosaccharide("0aU");
        assert_eq!(arabinofuranose.ring, RingClosure::Furanose);
        assert_eq!(arabinofuranose.skeleton_code, "211h");

        let glucosamine_sulfate = residue_to_monosaccharide("UYS");
        assert!(
            glucosamine_sulfate
                .modifications
                .iter()
                .any(|modification| modification.descriptor == "NSO/3=O/3=O")
        );

        let iduronate = residue_to_monosaccharide("IDR");
        assert_eq!(iduronate.skeleton_code, "2121A");

        assert_eq!(residue_to_monosaccharide("3hA").skeleton_code, "2211m");
        assert_eq!(residue_to_monosaccharide("3HA").skeleton_code, "1122m");
        assert_eq!(residue_to_monosaccharide("FUC").skeleton_code, "1221m");
        assert!(!is_sugar_residue("RHM"));
        assert_eq!(residue_to_monosaccharide("XYS").skeleton_code, "212h");
        assert_eq!(
            residue_to_monosaccharide("XYS").anomeric_symbol,
            AnomericSymbol::Alpha
        );
        assert_eq!(residue_to_monosaccharide("8SA").anomeric_position, 2);
    }

    #[test]
    fn every_bundled_ccd_component_has_a_concrete_registry_mapping() {
        for line in PDB_COMPONENTS
            .lines()
            .filter(|line| !line.starts_with('#') && !line.trim().is_empty())
        {
            let component = line.split('\t').next().unwrap();
            let residue = pdb_component_residue(component)
                .unwrap_or_else(|| panic!("missing bundled CCD mapping for {component}"));
            assert!(
                residue.residue_kind.is_some(),
                "{component} lost its registry identity"
            );
            assert!(!residue.skeleton_code.contains('x'), "{component}");
        }
    }

    #[test]
    fn every_concrete_registry_name_is_recognized_for_mmcif_and_private_components() {
        let concrete = ResidueKind::ALL
            .iter()
            .copied()
            .filter(|kind| !kind.is_generic())
            .collect::<Vec<_>>();
        assert_eq!(concrete.len(), 74);
        for kind in concrete {
            let name = kind.canonical_name();
            assert!(is_sugar_residue(name), "{name}");
            assert_eq!(
                residue_to_monosaccharide(name).residue_kind,
                Some(kind),
                "{name}"
            );
        }
    }

    #[test]
    fn test_parse_minimal_pdb_with_sugar() {
        let lines = [
            "HEADER    TEST                                                            END"
                .to_string(),
            make_hetatm_line(1, "C1", "NAG", "A", 1, -1.0, 0.0, 0.0),
            make_hetatm_line(2, "O4", "NAG", "A", 1, 0.0, 0.0, 0.0),
            make_hetatm_line(3, "C1", "BMA", "A", 2, 1.4, 0.0, 0.0),
            "CONECT    2    3".to_string(),
            "END".to_string(),
        ];
        let pdb_str = lines.join("\n") + "\n";

        let result = extract_glycans_from_str(&pdb_str, false);
        assert!(result.is_ok(), "Failed: {:?}", result.err());
        let glycans = result.unwrap();
        assert_eq!(glycans.len(), 1);
        assert!(glycans[0].graph.node_count() >= 2);
    }

    #[test]
    fn mmcif_fixture_extracts_registry_sugar_with_provenance() {
        let mmcif = r#"data_crabwurcs
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.auth_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.auth_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.pdbx_formal_charge
_atom_site.pdbx_PDB_model_num
HETATM 1 C C1 . NAG A B 1 1 42 . -1.0 0.0 0.0 1.0 0.0 0 1
HETATM 2 O O4 . NAG A B 1 1 42 .  0.0 0.0 0.0 1.0 0.0 0 1
"#;
        let glycans = extract_glycans_with_provenance_from_str(mmcif, true).unwrap();
        assert_eq!(glycans.len(), 1);
        assert_eq!(glycans[0].residues.len(), 1);
        assert_eq!(glycans[0].residues[0].chain, "B");
        assert_eq!(glycans[0].residues[0].sequence_number, 42);
        let root = glycans[0].graph.root().unwrap();
        assert_eq!(
            crabwurcs_core::classify_residue(glycans[0].graph.residue(root).unwrap()),
            Some(ResidueKind::GlcNAc)
        );
    }

    #[test]
    fn provenance_tracks_source_chain_and_sequence_for_each_node() {
        let lines = [
            "HEADER    TEST                                                            END"
                .to_string(),
            make_hetatm_line(1, "C1", "NAG", "B", 17, -1.0, 0.0, 0.0),
            make_hetatm_line(2, "O4", "NAG", "B", 17, 0.0, 0.0, 0.0),
            make_hetatm_line(3, "C1", "BMA", "B", 18, 1.4, 0.0, 0.0),
            "CONECT    2    3".to_string(),
            "END".to_string(),
        ];
        let pdb_str = lines.join("\n") + "\n";
        let glycans = extract_glycans_with_provenance_from_str(&pdb_str, false).unwrap();
        assert_eq!(glycans.len(), 1);
        assert_eq!(glycans[0].residues.len(), glycans[0].graph.node_count());
        assert_eq!(
            glycans[0]
                .residues
                .iter()
                .map(|residue| (residue.chain.as_str(), residue.sequence_number))
                .collect::<Vec<_>>(),
            vec![("B", 17), ("B", 18)]
        );
        assert_eq!(
            glycans[0]
                .residues
                .iter()
                .map(|residue| residue.node_index)
                .collect::<Vec<_>>(),
            vec![0, 1]
        );
    }

    #[test]
    fn renamed_component_is_recognized_from_its_coordinate_graph() {
        let atoms = [
            ("O5", 1.400, 0.000, 0.000),
            ("C1", 0.700, 1.212, 0.150),
            ("C2", -0.700, 1.212, -0.150),
            ("C3", -1.400, 0.000, 0.150),
            ("C4", -0.700, -1.212, -0.150),
            ("C5", 0.700, -1.212, 0.150),
            ("O1", 1.400, 2.424, 0.650),
            ("O2", -1.400, 2.424, -0.650),
            ("O3", -2.800, 0.000, 0.650),
            ("O4", -1.400, -2.424, -0.650),
            ("C6", 1.450, -2.511, 0.350),
            ("O6", 2.150, -3.723, 0.850),
        ];
        let mut lines = vec!["HEADER    PRIVATE CARBOHYDRATE".to_string()];
        lines.extend(atoms.iter().enumerate().map(|(index, (name, x, y, z))| {
            make_hetatm_line(index as u32 + 1, name, "ZZZ", "A", 1, *x, *y, *z)
        }));
        lines.push("END".into());

        let glycans = extract_glycans_from_str(&(lines.join("\n") + "\n"), false).unwrap();
        assert_eq!(glycans.len(), 1);
        assert_eq!(glycans[0].graph.node_count(), 1);
        assert_ne!(
            glycans[0]
                .graph
                .residue(glycans[0].graph.root().unwrap())
                .unwrap()
                .display_name
                .as_deref(),
            Some("ZZZ")
        );
    }

    #[test]
    fn test_no_sugar_pdb() {
        let lines = [
            "HEADER    TEST                                                            END"
                .to_string(),
            make_hetatm_line(1, "CA", "ALA", "A", 1, -1.0, 0.0, 0.0),
            "END".to_string(),
        ];
        let pdb_str = lines.join("\n") + "\n";

        let result = extract_glycans_from_str(&pdb_str, false);
        assert!(result.is_ok(), "Failed: {:?}", result.err());
        assert!(result.unwrap().is_empty());
    }
}
