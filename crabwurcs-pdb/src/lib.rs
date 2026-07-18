use crabwurcs_core::{
    AnomericSymbol, CarbonPosition, Linkage, Modification, Monosaccharide, ResidueGraph,
    RingClosure,
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

#[derive(Debug, Clone)]
pub struct ExtractedGlycan {
    pub attachment_site: Option<String>,
    pub graph: ResidueGraph,
}

const SUGAR_RESIDUE_NAMES: &[&str] = &[
    "NAG", "NDG", "BMA", "MAN", "FUC", "GAL", "GLC", "SIA", "NGN", "XYS", "XYL", "RIP", "RIB",
    "ARA", "LYX", "ALL", "ALT", "GUL", "IDO", "TAL", "GLA", "FRU", "SOR", "PSI", "TAG", "XYU",
    "BGC", "BMX", "FCA", "FCB", "FUL", "GCU", "GCV", "G6D", "GIV", "GL0", "GNS", "GTR", "GXL",
    "GYC", "GYG", "GYP", "GYS", "HSG", "HSQ", "HSY", "KDM", "KDO", "KDN", "L6S", "LAT", "LFR",
    "MHS", "NGA", "NGC", "NGE", "NGS", "NM6", "OAA", "RAM", "SGA", "SGC", "SGD", "SGN", "SGS",
    "SIB", "SID", "SLB", "SUS", "T6S", "XYP", "XXR", "Z0E", "ZB0", "ZB1", "A2G", "GCS", "BDP",
    "IDR", "IDS", "RHM", "ARB", "GZL",
];

fn is_sugar_residue(name: &str) -> bool {
    SUGAR_RESIDUE_NAMES.contains(&name) || glycam_residue_kind(name).is_some()
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum SugarKind {
    Glc,
    GlcNAc,
    GlcN,
    GlcA,
    Qui,
    Gal,
    GalNAc,
    GalA,
    Man,
    ManNAc,
    Fuc,
    Rha,
    Ara,
    Xyl,
    Fru,
    IdoA,
    Gul,
    Alt,
    All,
    Tal,
    Rib,
    Lyx,
    Psi,
    Sor,
    Tag,
    GalF,
    Bac,
    Kdo,
    Neu5Ac,
    Neu5Gc,
    Kdn,
    Unknown,
}

fn glycam_residue_kind(name: &str) -> Option<(SugarKind, AnomericSymbol)> {
    let bytes = name.as_bytes();
    if bytes.len() != 3 || !b"0123456789ZYXWVUTSRQP".contains(&bytes[0]) {
        return None;
    }
    let middle = bytes[1] as char;
    let last = bytes[2] as char;
    let kind = match (middle, last) {
        ('G' | 'g', 'L' | 'l') => SugarKind::Neu5Gc,
        ('Y' | 'y', 'N' | 'n' | 'S' | 's') => SugarKind::GlcN,
        ('K' | 'k', 'N' | 'n') => SugarKind::Kdn,
        ('K' | 'k', 'O' | 'o') => SugarKind::Kdo,
        ('B' | 'b', 'C' | 'c') => SugarKind::Bac,
        ('Z' | 'z', _) => SugarKind::GlcA,
        ('O' | 'o', _) => SugarKind::GalA,
        ('U' | 'u', _) => SugarKind::IdoA,
        ('C' | 'c', _) => SugarKind::Fru,
        ('A' | 'a', _) => SugarKind::Ara,
        ('H' | 'h', _) => SugarKind::Rha,
        ('Q' | 'q', _) => SugarKind::Qui,
        ('W' | 'w', _) => SugarKind::ManNAc,
        ('K' | 'k', _) => SugarKind::Gul,
        ('E' | 'e', _) => SugarKind::Alt,
        ('N' | 'n', _) => SugarKind::All,
        ('T' | 't', _) => SugarKind::Tal,
        ('P' | 'p', 'D' | 'U') => SugarKind::Psi,
        ('R' | 'r', 'D' | 'U') => SugarKind::Rib,
        ('D' | 'd', 'A' | 'B') => SugarKind::Lyx,
        ('J' | 'j', 'A' | 'B') => SugarKind::Tag,
        ('b', 'A' | 'B') => SugarKind::Sor,
        _ => match middle {
            'G' | 'g' => SugarKind::Glc,
            'Y' | 'y' => SugarKind::GlcNAc,
            'L' | 'l' => SugarKind::Gal,
            'V' | 'v' => SugarKind::GalNAc,
            'M' | 'm' => SugarKind::Man,
            'F' | 'f' => SugarKind::Fuc,
            'X' | 'x' => SugarKind::Xyl,
            'S' | 's' => SugarKind::Neu5Ac,
            _ => return None,
        },
    };
    let anomer = match last {
        'A' | 'D' => AnomericSymbol::Alpha,
        'B' | 'U' => AnomericSymbol::Beta,
        'N' | 'S' | 'L' => AnomericSymbol::Alpha,
        'n' | 's' | 'l' => AnomericSymbol::Beta,
        _ if matches!(kind, SugarKind::Bac | SugarKind::Kdo | SugarKind::Kdn) => {
            AnomericSymbol::Alpha
        }
        _ => return None,
    };
    Some((kind, anomer))
}

fn pdb_residue_kind(name: &str) -> (SugarKind, AnomericSymbol) {
    match name {
        "NAG" => (SugarKind::GlcNAc, AnomericSymbol::Beta),
        "NDG" => (SugarKind::GlcNAc, AnomericSymbol::Alpha),
        "NGA" => (SugarKind::GalNAc, AnomericSymbol::Beta),
        "A2G" => (SugarKind::GalNAc, AnomericSymbol::Alpha),
        "BMA" => (SugarKind::Man, AnomericSymbol::Beta),
        "MAN" => (SugarKind::Man, AnomericSymbol::Alpha),
        "FUC" | "FCA" => (SugarKind::Fuc, AnomericSymbol::Alpha),
        "FUL" | "FCB" => (SugarKind::Fuc, AnomericSymbol::Beta),
        "GAL" => (SugarKind::Gal, AnomericSymbol::Beta),
        "GLA" => (SugarKind::Gal, AnomericSymbol::Alpha),
        "BGC" => (SugarKind::Glc, AnomericSymbol::Beta),
        "GLC" => (SugarKind::Glc, AnomericSymbol::Alpha),
        "BDP" | "GCU" => (SugarKind::GlcA, AnomericSymbol::Beta),
        "IDR" | "IDS" => (SugarKind::IdoA, AnomericSymbol::Alpha),
        "RAM" | "RHM" => (SugarKind::Rha, AnomericSymbol::Alpha),
        "ARB" => (SugarKind::Ara, AnomericSymbol::Beta),
        "GZL" => (SugarKind::GalF, AnomericSymbol::Beta),
        "GCS" | "GNS" | "QYS" | "UYS" | "VYS" => (SugarKind::GlcN, AnomericSymbol::Alpha),
        "FRU" => (SugarKind::Fru, AnomericSymbol::Beta),
        "ARA" => (SugarKind::Ara, AnomericSymbol::Beta),
        "SIA" | "SLB" | "SUS" => (SugarKind::Neu5Ac, AnomericSymbol::Alpha),
        "SIB" | "SID" => (SugarKind::Neu5Ac, AnomericSymbol::Beta),
        "NGN" | "NGC" => (SugarKind::Neu5Gc, AnomericSymbol::Alpha),
        "KDN" => (SugarKind::Kdn, AnomericSymbol::Alpha),
        "KDO" => (SugarKind::Kdo, AnomericSymbol::Alpha),
        "XYS" => (SugarKind::Xyl, AnomericSymbol::Alpha),
        "XYL" | "XYP" => (SugarKind::Xyl, AnomericSymbol::Beta),
        _ => (SugarKind::Unknown, AnomericSymbol::Unknown),
    }
}

fn modification(position: u8, descriptor: &str) -> Modification {
    Modification {
        position: CarbonPosition(position),
        descriptor: descriptor.to_string(),
        probability: None,
    }
}

fn residue_to_monosaccharide(name: &str) -> Monosaccharide {
    let pdb_kind = pdb_residue_kind(name);
    let (kind, anomer) = if pdb_kind.0 != SugarKind::Unknown {
        pdb_kind
    } else {
        glycam_residue_kind(name).unwrap_or(pdb_kind)
    };
    let l_isomer = if pdb_kind.0 != SugarKind::Unknown {
        matches!(name, "FUC" | "FCA" | "FUL" | "FCB" | "RAM")
    } else {
        name.as_bytes().get(1).is_some_and(u8::is_ascii_lowercase)
    };
    let furanose_code = name
        .as_bytes()
        .get(2)
        .is_some_and(|code| matches!(*code as char, 'D' | 'U'));
    let ring = if kind == SugarKind::Bac {
        RingClosure::Unknown
    } else if matches!(
        kind,
        SugarKind::Fru | SugarKind::GalF | SugarKind::Rib | SugarKind::Psi
    ) || (kind == SugarKind::Ara && furanose_code)
    {
        RingClosure::Furanose
    } else {
        RingClosure::Pyranose
    };
    let (prefix, skeleton, ring_start, ring_end, position, modifications) = match kind {
        SugarKind::Glc => ("a", "2122h", 1, 5, 1, vec![]),
        SugarKind::GlcNAc => ("a", "2122h", 1, 5, 1, vec![modification(2, "NCC/3=O")]),
        SugarKind::GlcN => (
            "a",
            "2122h",
            1,
            5,
            1,
            vec![modification(
                2,
                if name.ends_with(['S', 's']) {
                    "NSO/3=O/3=O"
                } else {
                    "N"
                },
            )],
        ),
        SugarKind::GlcA => ("a", "2122Ah", 1, 5, 1, vec![]),
        SugarKind::Qui => ("a", "2122m", 1, 5, 1, vec![]),
        SugarKind::Gal => ("a", "2112h", 1, 5, 1, vec![]),
        SugarKind::GalNAc => ("a", "2112h", 1, 5, 1, vec![modification(2, "NCC/3=O")]),
        SugarKind::GalA => ("a", "2112Ah", 1, 5, 1, vec![]),
        SugarKind::Man => (
            "a",
            if l_isomer { "2211h" } else { "1122h" },
            1,
            5,
            1,
            vec![],
        ),
        SugarKind::ManNAc => ("a", "1122h", 1, 5, 1, vec![modification(2, "NCC/3=O")]),
        SugarKind::Fuc => (
            "a",
            if l_isomer { "1221m" } else { "2112m" },
            1,
            5,
            1,
            vec![],
        ),
        SugarKind::Rha => (
            "a",
            if l_isomer { "2211m" } else { "1122m" },
            1,
            5,
            1,
            vec![],
        ),
        SugarKind::Ara => (
            "a",
            if name.as_bytes().get(1).is_some_and(u8::is_ascii_lowercase) {
                "211h"
            } else {
                "122h"
            },
            1,
            if furanose_code { 4 } else { 5 },
            1,
            vec![],
        ),
        SugarKind::Xyl => ("a", "212h", 1, 5, 1, vec![]),
        SugarKind::Fru => ("ha", "122h", 2, 5, 2, vec![]),
        SugarKind::IdoA => ("a", "2121Ah", 1, 5, 1, vec![]),
        SugarKind::Gul => ("a", "1121h", 1, 5, 1, vec![]),
        SugarKind::Alt => ("a", "2111h", 1, 5, 1, vec![]),
        SugarKind::All => ("a", "2222h", 1, 5, 1, vec![]),
        SugarKind::Tal => ("a", "2221h", 1, 5, 1, vec![]),
        SugarKind::Rib => ("a", "222h", 1, 4, 1, vec![]),
        SugarKind::Lyx => ("a", "221h", 1, 5, 1, vec![]),
        SugarKind::Psi => ("ha", "222h", 2, 5, 2, vec![]),
        SugarKind::Sor => ("ha", "121h", 2, 6, 2, vec![]),
        SugarKind::Tag => ("ha", "112h", 2, 6, 2, vec![]),
        SugarKind::GalF => ("a", "2112h", 1, 4, 1, vec![]),
        SugarKind::Bac => (
            "a",
            "xxxxm",
            1,
            5,
            1,
            vec![modification(2, "NCC/3=O"), modification(4, "NCC/3=O")],
        ),
        SugarKind::Kdo => ("A", "1122h", 2, 6, 2, vec![]),
        SugarKind::Neu5Ac => ("Aad", "21122h", 2, 6, 2, vec![modification(5, "NCC/3=O")]),
        SugarKind::Neu5Gc => ("Aad", "21122h", 2, 6, 2, vec![modification(5, "NCCO/3=O")]),
        SugarKind::Kdn => ("Aad", "21122h", 2, 6, 2, vec![]),
        SugarKind::Unknown => ("a", "xxxxh", 1, 5, 1, vec![]),
    };
    Monosaccharide::new(
        skeleton.chars().filter(char::is_ascii_digit).count() as u8,
        skeleton.into(),
        vec![],
        ring,
        Some(ring_start),
        Some(ring_end),
        position,
        anomer,
        prefix.into(),
        modifications,
    )
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

    extract_glycans_from_pdb(
        &pdb,
        &raw_pdb_residue_names(&contents),
        &raw_pdb_bonds(&contents),
    )
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

    extract_glycans_from_pdb(
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

fn extract_glycans_from_pdb(
    pdb: &pdbtbx::PDB,
    raw_names: &HashMap<(isize, String), String>,
    raw_bonds: &[(usize, usize)],
) -> PdbResult<Vec<ExtractedGlycan>> {
    #[derive(Debug, Clone)]
    struct ResidueMeta {
        name: String,
        sequence: isize,
        chain: String,
        order: usize,
    }

    #[derive(Debug, Clone)]
    struct AtomMeta {
        residue: usize,
        name: String,
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
            let sequence = hierarchy.residue().id().0;
            ResidueMeta {
                name: raw_names
                    .get(&(sequence, parsed_name.to_string()))
                    .cloned()
                    .unwrap_or_else(|| parsed_name.to_string()),
                sequence,
                chain: hierarchy.chain().id().to_string(),
                order: next_order,
            }
        });
        atoms.insert(
            atom_key,
            AtomMeta {
                residue: residue_key,
                name: hierarchy.atom().name().to_string(),
                position: [
                    hierarchy.atom().x(),
                    hierarchy.atom().y(),
                    hierarchy.atom().z(),
                ],
            },
        );
        atoms_by_serial.insert(hierarchy.atom().serial_number(), atom_key);
    }

    let sugar_keys = residues
        .iter()
        .filter_map(|(key, residue)| is_sugar_residue(&residue.name).then_some(*key))
        .collect::<HashSet<_>>();
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
        for key in &ordered {
            let mut residue = residue_to_monosaccharide(&residues[key].name);
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
        extracted.push(ExtractedGlycan {
            attachment_site,
            graph,
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
        format!(
            "HETATM{:5} {:<4}{:1}{:3} {}{:4}{:1}   {:8.3}{:8.3}{:8.3}{:6.2}{:6.2}          {:>2}{:2}",
            serial, atom_name, "", res_name, chain, seq, "", x, y, z, 1.0f64, 0.0f64, "C", ""
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
            assert!(neu5gc
                .modifications
                .iter()
                .any(|modification| modification.position.0 == 5
                    && modification.descriptor == "NCCO/3=O"));
        }

        let arabinofuranose = residue_to_monosaccharide("0aU");
        assert_eq!(arabinofuranose.ring, RingClosure::Furanose);
        assert_eq!(arabinofuranose.skeleton_code, "211h");

        let glucosamine_sulfate = residue_to_monosaccharide("UYS");
        assert!(glucosamine_sulfate
            .modifications
            .iter()
            .any(|modification| modification.descriptor == "NSO/3=O/3=O"));

        let iduronate = residue_to_monosaccharide("IDR");
        assert_eq!(iduronate.skeleton_code, "2121Ah");

        assert_eq!(residue_to_monosaccharide("3hA").skeleton_code, "2211m");
        assert_eq!(residue_to_monosaccharide("3HA").skeleton_code, "1122m");
        assert_eq!(residue_to_monosaccharide("FUC").skeleton_code, "1221m");
        assert_eq!(residue_to_monosaccharide("RHM").skeleton_code, "1122m");
        assert_eq!(residue_to_monosaccharide("XYS").skeleton_code, "212h");
        assert_eq!(
            residue_to_monosaccharide("XYS").anomeric_symbol,
            AnomericSymbol::Alpha
        );
        assert_eq!(residue_to_monosaccharide("8SA").anomeric_position, 2);
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
