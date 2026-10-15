//! What residue names say: the bundled wwPDB CCD component table, concrete
//! registry names and GLYCAM residue codes.
//!
//! These readings are hints. The coordinates decide (see `geometry`); a name
//! fills in only what the coordinates cannot show, and a name that
//! contradicts them is reported.

use crabwurcs_core::{AnomericSymbol, Monosaccharide, ResidueKind, RingClosure, residue_from_kind};

pub(crate) const PDB_COMPONENTS: &str = include_str!("../data/pdb_carbohydrate_components.tsv");

/// Where a name reading comes from.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum NameSource {
    /// The bundled wwPDB CCD table: the residue's kind, anomer and ring.
    /// The table does not list O-substituents (IDS reads as IdoA).
    Ccd,
    /// A concrete registry name such as `Neu5Ac`.
    Registry,
    /// A GLYCAM residue code. GLYCAM keeps O-substituents (`SO3`, `MEX`,
    /// `ACX`) in residues of their own, so the code describes every
    /// substituent carried by the residue's own atoms.
    Glycam,
}

#[derive(Debug, Clone)]
pub(crate) struct NameReading {
    pub residue: Monosaccharide,
    pub source: NameSource,
}

impl NameReading {
    /// Whether the name lists every substituent made of the residue's own
    /// atoms, so that an extra one in the coordinates contradicts it.
    pub fn lists_own_substituents(&self) -> bool {
        self.source == NameSource::Glycam
    }

    /// Whether an amino group of the name may carry a group the name does
    /// not spell out: the CCD table and registry names give the parent kind
    /// (SGN is GlcN in the bundled table), and GLYCAM writes the N-acetyl
    /// groups of bacillosamine as modifications of its `BC` residue
    /// (`DBacp[2Ac,4Ac]`).
    pub fn amines_may_carry_groups(&self) -> bool {
        self.source != NameSource::Glycam || self.residue.residue_kind == Some(ResidueKind::Bac)
    }
}

pub(crate) fn is_sugar_name(name: &str, hetero: bool) -> bool {
    pdb_component_residue(name).is_some()
        // The GLYCAM code lookup is a loose heuristic (any 3-letter name can
        // look like a code), so in ATOM records it is applied only to codes
        // that begin with a linkage digit (`0aU`, `4YS`): amino acid names
        // such as VAL collide with GLYCAM codes and must never be treated as
        // glycans. GlycoShape's files write every residue as ATOM.
        || ((hetero || name.starts_with(|character: char| character.is_ascii_digit()))
            && glycam_residue_kind(name).is_some())
        || registry_named_residue(name).is_some()
}

/// Whether the name is a CCD component or a registry name: names that may
/// stand for a sugar even when its ring atoms are missing. A GLYCAM-looking
/// code alone is not enough; GLYCAM files are complete.
pub(crate) fn is_table_name(name: &str) -> bool {
    pdb_component_residue(name).is_some() || registry_named_residue(name).is_some()
}

pub(crate) fn registry_named_residue(name: &str) -> Option<Monosaccharide> {
    let kind = ResidueKind::from_name(name)?;
    (!kind.is_generic())
        .then(|| residue_from_kind(kind).ok())
        .flatten()
}

pub(crate) fn apply_declared_form(
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

pub(crate) fn pdb_component_residue(name: &str) -> Option<Monosaccharide> {
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
    residue.skeleton_code = crabwurcs_core::mirror_skeleton_code(&residue.skeleton_code);
    // The registry hint names the template's configuration, not its mirror
    // image (L-Fuc's template mirrored is D-Fuc): classify by chemistry.
    residue.residue_kind = None;
}

/// Decode a GLYCAM residue code: linkage character, sugar letter, form.
///
/// The case of the sugar letter is the configuration: upper case D, lower
/// case L (`0GA` α-D-Glcp, `0uA` α-L-IdopA, `0fA` α-L-Fucp). The registry
/// template is mirrored when the code's configuration is not the template's
/// own (the templates of Ara, Fuc, Rha, IdoA, Alt and Sor are L).
pub(crate) fn glycam_residue_kind(name: &str) -> Option<Monosaccharide> {
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
        // The Kdo, Kdn and Bac codes carry no anomer letter.
        _ if matches!(kind, ResidueKind::Bac | ResidueKind::Kdo | ResidueKind::Kdn) => {
            AnomericSymbol::Unknown
        }
        _ => return None,
    };
    let mut residue = residue_from_kind(kind).ok()?;
    // Sorbose has only the lower-case code; sialic acid letters are not a
    // configuration.
    let configuration = if matches!(
        kind,
        ResidueKind::Neu5Ac | ResidueKind::Neu5Gc | ResidueKind::Sor
    ) {
        kind.default_configuration()
    } else if middle.is_ascii_uppercase() {
        Some('D')
    } else {
        Some('L')
    };
    if configuration.is_some() && configuration != kind.default_configuration() {
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

/// The name reading of a residue, if its name is known.
pub(crate) fn name_reading(name: &str) -> Option<NameReading> {
    if let Some(residue) = pdb_component_residue(name) {
        return Some(NameReading {
            residue,
            source: NameSource::Ccd,
        });
    }
    if let Some(residue) = registry_named_residue(name) {
        return Some(NameReading {
            residue,
            source: NameSource::Registry,
        });
    }
    glycam_residue_kind(name).map(|residue| NameReading {
        residue,
        source: NameSource::Glycam,
    })
}

#[cfg(test)]
pub(crate) fn residue_to_monosaccharide(name: &str) -> Monosaccharide {
    name_reading(name)
        .map(|reading| reading.residue)
        .unwrap_or_else(|| {
            residue_from_kind(ResidueKind::Hex).expect("the generic hexose registry entry is valid")
        })
}
