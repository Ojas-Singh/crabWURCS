//! Glycans from PDB and mmCIF coordinate files.
//!
//! Every monosaccharide is read from its 3D coordinates: the ring and the
//! carbon numbering come from the bond graph, the configuration of every
//! stereocentre (and so D/L, the parent sugar and the anomer) from signed
//! volumes, and the substituents and glycosidic linkages from what the ring
//! atoms are bonded to. Residue names are hints: they fill in what the
//! coordinates cannot show (missing atoms, a centre too distorted to read),
//! and every disagreement between a name and the coordinates is reported as
//! a [`PdbWarning`]. [`NamePolicy`] chooses what happens on a disagreement.
//!
//! A residue whose name is a known sugar but whose ring atoms are missing is
//! read from its name, as before.

mod geometry;
mod names;
mod reading;

use crabwurcs_core::{
    CarbonPosition, Linkage, Modification, Monosaccharide, ResidueGraph, residue_from_kind,
};
use geometry::{CarbonReading, Group, RingUnit, Structure, StructureAtom};
use names::{NameReading, is_sugar_name, is_table_name, name_reading};
#[cfg(test)]
use names::{PDB_COMPONENTS, pdb_component_residue, residue_to_monosaccharide};
use pdbtbx::{
    ContainsAtomConformer, ContainsAtomConformerResidue, ContainsAtomConformerResidueChain,
};
use reading::{CoordinateReading, CoordinateSubstituent, residue_label, unique_residue};
use std::collections::{BTreeMap, HashMap, HashSet};
use thiserror::Error;

#[derive(Debug, Error)]
pub enum PdbError {
    #[error("failed to parse structure file: {0}")]
    ParseError(String),

    #[error("no carbohydrate residues found in structure")]
    NoGlycansFound,

    #[error(transparent)]
    Core(#[from] crabwurcs_core::CoreError),

    /// [`NamePolicy::Strict`]: residue names contradict the coordinates.
    #[error("residue names contradict the coordinates:\n{}", format_warnings(.0))]
    NamesContradictCoordinates(Vec<PdbWarning>),
}

fn format_warnings(warnings: &[PdbWarning]) -> String {
    warnings
        .iter()
        .map(|warning| format!("  {warning}"))
        .collect::<Vec<_>>()
        .join("\n")
}

pub type PdbResult<T> = Result<T, PdbError>;

/// What decides when a residue name and the coordinates disagree.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Hash)]
pub enum NamePolicy {
    /// The coordinates decide; a name fills in only what the coordinates
    /// cannot show. Every disagreement is reported.
    #[default]
    Coordinates,
    /// A name decides wherever it says something; the coordinates fill in
    /// the rest (unknown names, substituents the name does not describe).
    /// Every disagreement is still reported.
    Names,
    /// As [`NamePolicy::Coordinates`], but a disagreement is an error
    /// ([`PdbError::NamesContradictCoordinates`]).
    Strict,
}

impl std::str::FromStr for NamePolicy {
    type Err = String;

    fn from_str(value: &str) -> Result<Self, Self::Err> {
        match value.trim().to_ascii_lowercase().as_str() {
            "coordinates" => Ok(Self::Coordinates),
            "names" => Ok(Self::Names),
            "strict" => Ok(Self::Strict),
            other => Err(format!(
                "unknown name policy `{other}` (expected coordinates, names or strict)"
            )),
        }
    }
}

/// Options for [`read_glycans_from_file`] and [`read_glycans_from_str`].
#[derive(Debug, Clone, Default, PartialEq, Eq)]
#[non_exhaustive]
pub struct PdbReadOptions {
    pub policy: NamePolicy,
}

impl PdbReadOptions {
    pub fn with_policy(policy: NamePolicy) -> Self {
        Self { policy }
    }
}

/// The kind of a [`PdbWarning`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum PdbWarningKind {
    /// The residue name describes another monosaccharide than the
    /// coordinates: another configuration at a carbon, another anomer, ring
    /// size, backbone or substituent.
    NameContradictsCoordinates,
    /// Atoms needed to read part of the residue are missing; that part was
    /// taken from the residue name.
    MissingAtoms,
    /// A centre is too close to planar to read its configuration; it was
    /// taken from the residue name, or left undetermined.
    DistortedGeometry,
    /// A carbon has a valence pattern the reader does not interpret.
    UnreadableGeometry,
    /// No monosaccharide ring was found in the coordinates of a residue
    /// named as a sugar; it was read from its name.
    NoRing,
    /// A group on a ring heteroatom is not a substituent the reader knows; it
    /// was left out.
    UnrecognizedGroup,
}

impl PdbWarningKind {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::NameContradictsCoordinates => "name-contradicts-coordinates",
            Self::MissingAtoms => "missing-atoms",
            Self::DistortedGeometry => "distorted-geometry",
            Self::UnreadableGeometry => "unreadable-geometry",
            Self::NoRing => "no-ring",
            Self::UnrecognizedGroup => "unrecognized-group",
        }
    }
}

/// Something about one residue that its name and its coordinates do not
/// agree on, or that the coordinates cannot decide.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PdbWarning {
    /// Index of the glycan in [`PdbReading::glycans`].
    pub glycan_index: usize,
    /// Graph node of the residue in that glycan.
    pub node_index: usize,
    pub chain: String,
    pub sequence_number: isize,
    pub insertion_code: Option<String>,
    pub residue_name: String,
    pub kind: PdbWarningKind,
    /// What the residue name says (`α-L-IdopA`), when the name is known.
    pub named: Option<String>,
    /// What the coordinates say (`α-L-GulpA2S`).
    pub coordinates: Option<String>,
    /// Carbon positions concerned; empty for the whole residue.
    pub positions: Vec<u8>,
    /// What differs: `C2`, `anomer`, `ring`, `backbone`, `C2 substituent`.
    pub aspects: Vec<String>,
}

impl PdbWarning {
    fn location(&self) -> String {
        format!(
            "residue {} {}:{}{}",
            self.residue_name,
            self.chain,
            self.sequence_number,
            self.insertion_code.as_deref().unwrap_or("")
        )
    }
}

impl std::fmt::Display for PdbWarning {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let location = self.location();
        let aspects = self.aspects.join(", ");
        let named = self.named.as_deref().unwrap_or("?");
        let coordinates = self.coordinates.as_deref().unwrap_or("?");
        match self.kind {
            PdbWarningKind::NameContradictsCoordinates => write!(
                formatter,
                "{location}: named {named}, coordinates say {coordinates} ({aspects})"
            ),
            PdbWarningKind::MissingAtoms if self.named.is_some() => write!(
                formatter,
                "{location}: atoms missing ({aspects}); taken from the name {named}"
            ),
            PdbWarningKind::MissingAtoms => write!(
                formatter,
                "{location}: atoms missing ({aspects}); left undetermined"
            ),
            PdbWarningKind::DistortedGeometry if self.named.is_some() => write!(
                formatter,
                "{location}: too distorted to read ({aspects}); taken from the name {named}"
            ),
            PdbWarningKind::DistortedGeometry => write!(
                formatter,
                "{location}: too distorted to read ({aspects}); left undetermined"
            ),
            PdbWarningKind::UnreadableGeometry if self.named.is_some() => write!(
                formatter,
                "{location}: not readable from the coordinates ({aspects}); taken from the name {named}"
            ),
            PdbWarningKind::UnreadableGeometry => write!(
                formatter,
                "{location}: not readable from the coordinates ({aspects}); left undetermined"
            ),
            PdbWarningKind::NoRing => write!(
                formatter,
                "{location}: no sugar ring in the coordinates; read from its name as {named}"
            ),
            PdbWarningKind::UnrecognizedGroup => write!(
                formatter,
                "{location}: unrecognised group ({aspects}); left out"
            ),
        }
    }
}

/// One configuration read from a signed volume.
#[derive(Debug, Clone, PartialEq)]
pub struct StereoCentre {
    /// Carbon position.
    pub position: u8,
    /// The atom's name in the file (reported, never used to decide).
    pub atom: String,
    /// Fischer side: `'1'` (left) or `'2'` (right). For the anomeric carbon
    /// this is the side of its exocyclic heteroatom.
    pub code: char,
    /// Signed volume normalised by the three bond lengths (about ±0.77 for a
    /// tetrahedral centre, 0 for a planar one).
    pub volume: f64,
}

/// What the coordinates and the name say about one glycan residue.
#[derive(Debug, Clone, PartialEq)]
pub struct ResidueReading {
    pub glycan_index: usize,
    pub node_index: usize,
    pub chain: String,
    pub sequence_number: isize,
    pub insertion_code: Option<String>,
    pub residue_name: String,
    /// WURCS UniqueRES read from the coordinates alone (`x` where they cannot
    /// decide); `None` when no ring was found.
    pub from_coordinates: Option<String>,
    /// WURCS UniqueRES the residue name describes, when the name is known.
    pub from_name: Option<String>,
    /// Ring atoms: the ring oxygen, then the anomeric carbon and around.
    pub ring_atoms: Vec<String>,
    /// Backbone carbons C1..Cn.
    pub backbone_atoms: Vec<String>,
    pub centres: Vec<StereoCentre>,
}

/// The result of [`read_glycans_from_file`] / [`read_glycans_from_str`].
#[derive(Debug, Clone, Default)]
pub struct PdbReading {
    pub glycans: Vec<ExtractedGlycanWithProvenance>,
    /// Every disagreement between residue names and coordinates, and every
    /// part the coordinates could not decide.
    pub warnings: Vec<PdbWarning>,
    /// One record per graph node of every glycan.
    pub residues: Vec<ResidueReading>,
}

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

fn modification(position: u8, descriptor: &str) -> Modification {
    Modification {
        position: CarbonPosition(position),
        descriptor: descriptor.to_string(),
        probability: None,
    }
}

fn without_provenance(glycans: Vec<ExtractedGlycanWithProvenance>) -> Vec<ExtractedGlycan> {
    glycans
        .into_iter()
        .map(|glycan| ExtractedGlycan {
            attachment_site: glycan.attachment_site,
            graph: glycan.graph,
        })
        .collect()
}

/// Extract glycans (coordinates decide; warnings are dropped). See
/// [`read_glycans_from_file`] for the warnings and the policy.
pub fn extract_glycans_from_file(path: &std::path::Path) -> PdbResult<Vec<ExtractedGlycan>> {
    extract_glycans_with_provenance_from_file(path).map(without_provenance)
}

/// Extract glycans (coordinates decide; warnings are dropped). See
/// [`read_glycans_from_str`] for the warnings and the policy.
pub fn extract_glycans_from_str(contents: &str, is_mmcif: bool) -> PdbResult<Vec<ExtractedGlycan>> {
    extract_glycans_with_provenance_from_str(contents, is_mmcif).map(without_provenance)
}

/// The PDB records the extractor reads: the coordinate records (`ATOM`,
/// `HETATM`, `TER`) of the first model and the `CONECT` and `LINK` records,
/// each once, followed by `END`.
///
/// Everything else is skipped. Free-text records such as the banner
/// `REMARK` lines of GlycoShape's multi-model files are not in the column
/// layout the PDB format prescribes and would otherwise stop the parser,
/// and later models repeat the atoms of the first.
pub fn pdb_records_used(contents: &str) -> String {
    let mut kept = String::with_capacity(contents.len());
    let mut seen_bonds = HashSet::new();
    let mut first_model_done = false;
    for line in contents.lines() {
        let record = line.get(..6).unwrap_or(line).trim_end();
        let keep = match record {
            "ATOM" | "HETATM" | "TER" => !first_model_done,
            "MODEL" => !first_model_done,
            "ENDMDL" => {
                let first = !first_model_done;
                first_model_done = true;
                first
            }
            "CONECT" | "LINK" => seen_bonds.insert(line.trim_end()),
            _ => false,
        };
        if keep {
            kept.push_str(line.trim_end());
            kept.push('\n');
        }
    }
    kept.push_str("END\n");
    kept
}

fn is_mmcif_text(contents: &str) -> bool {
    contents.trim_start().starts_with("data_")
}

/// Extract glycans from a PDB/mmCIF file while retaining source residue
/// identities for every graph node (coordinates decide; warnings are
/// dropped). PDB files are read through [`pdb_records_used`].
pub fn extract_glycans_with_provenance_from_file(
    path: &std::path::Path,
) -> PdbResult<Vec<ExtractedGlycanWithProvenance>> {
    read_glycans_from_file(path, &PdbReadOptions::default()).map(|reading| reading.glycans)
}

/// Extract glycans from a PDB or mmCIF string while retaining source residue
/// identities for every graph node (coordinates decide; warnings are
/// dropped). PDB text is read through [`pdb_records_used`].
pub fn extract_glycans_with_provenance_from_str(
    contents: &str,
    is_mmcif: bool,
) -> PdbResult<Vec<ExtractedGlycanWithProvenance>> {
    read_glycans_from_str(contents, is_mmcif, &PdbReadOptions::default())
        .map(|reading| reading.glycans)
}

/// Read the glycans of a PDB/mmCIF file with their warnings and per-residue
/// readings. Only the first model is read.
pub fn read_glycans_from_file(
    path: &std::path::Path,
    options: &PdbReadOptions,
) -> PdbResult<PdbReading> {
    let contents = std::fs::read_to_string(path)
        .map_err(|error| PdbError::ParseError(format!("cannot read file: {error}")))?;
    let mmcif_extension = path
        .extension()
        .and_then(|extension| extension.to_str())
        .is_some_and(|extension| {
            extension.eq_ignore_ascii_case("cif") || extension.eq_ignore_ascii_case("mmcif")
        });
    if mmcif_extension {
        use pdbtbx::ReadOptions;

        let path_str = path.to_string_lossy();
        let (pdb, _errors) = ReadOptions::default()
            .set_level(pdbtbx::StrictnessLevel::Loose)
            .set_only_first_model(true)
            .read(path_str.as_ref())
            .map_err(|e| PdbError::ParseError(format!("cannot parse file: {:?}", e)))?;
        return read_structure(
            &pdb,
            &raw_pdb_residue_names(&contents),
            &raw_pdb_bonds(&contents),
            options,
        );
    }
    read_glycans_from_str(&contents, is_mmcif_text(&contents), options)
}

/// Read the glycans of PDB or mmCIF text with their warnings and per-residue
/// readings. PDB text is read through [`pdb_records_used`]; of mmCIF, the
/// first model is read.
pub fn read_glycans_from_str(
    contents: &str,
    is_mmcif: bool,
    options: &PdbReadOptions,
) -> PdbResult<PdbReading> {
    use pdbtbx::{Format, ReadOptions, StrictnessLevel};
    use std::io::BufReader;

    let filtered;
    let contents = if is_mmcif {
        contents
    } else {
        filtered = pdb_records_used(contents);
        filtered.as_str()
    };
    let format = if is_mmcif { Format::Mmcif } else { Format::Pdb };
    let reader = BufReader::new(contents.as_bytes());
    let (pdb, _errors) = ReadOptions::default()
        .set_level(StrictnessLevel::Loose)
        .set_format(format)
        .set_only_first_model(true)
        .read_raw(reader)
        .map_err(|e| PdbError::ParseError(format!("cannot parse: {:?}", e)))?;
    read_structure(
        &pdb,
        &raw_pdb_residue_names(contents),
        &raw_pdb_bonds(contents),
        options,
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
    let mut bonds = bonds.into_iter().collect::<Vec<_>>();
    bonds.sort_unstable();
    bonds
}

fn raw_pdb_residue_names(contents: &str) -> HashMap<(isize, String), String> {
    let mut names = HashMap::new();
    for line in contents.lines() {
        if !(line.starts_with("ATOM  ") || line.starts_with("HETATM")) || line.len() < 26 {
            continue;
        }
        let (Some(raw), Some(sequence)) = (line.get(17..20), line.get(22..26)) else {
            continue;
        };
        let raw = raw.trim();
        let Ok(sequence) = sequence.trim().parse::<isize>() else {
            continue;
        };
        names
            .entry((sequence, raw.to_ascii_uppercase()))
            .or_insert_with(|| raw.to_string());
    }
    names
}

fn element_from_name(name: &str) -> String {
    name.trim_start_matches(|character: char| character.is_ascii_digit())
        .chars()
        .next()
        .map(|character| character.to_ascii_uppercase().to_string())
        .unwrap_or_default()
}

#[derive(Debug, Clone)]
struct ResidueMeta {
    name: String,
    sequence: isize,
    insertion_code: Option<String>,
    chain: String,
    hetero: bool,
}

impl ResidueMeta {
    fn site(&self) -> String {
        format!("{}/{}/{}", self.chain, self.name, self.sequence)
    }
}

/// One monosaccharide of the structure: read from a ring in the coordinates,
/// or (ring atoms missing) from its residue name only.
struct Unit {
    residue: usize,
    ring: Option<RingUnit>,
    name: Option<NameReading>,
}

/// Index of a ring unit's backbone carbons: atom → (unit, position).
type BackboneIndex = HashMap<usize, (usize, u8)>;

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
struct Link {
    parent: usize,
    child: usize,
    parent_position: u8,
    child_position: u8,
}

fn read_structure(
    pdb: &pdbtbx::PDB,
    raw_names: &HashMap<(isize, String), String>,
    raw_bonds: &[(usize, usize)],
    options: &PdbReadOptions,
) -> PdbResult<PdbReading> {
    // Atoms of the first model; of alternative locations, the atoms without
    // one and those of the first.
    let mut residues: Vec<ResidueMeta> = Vec::new();
    let mut residue_index: HashMap<usize, usize> = HashMap::new();
    let mut first_alternative: HashMap<usize, Option<usize>> = HashMap::new();
    let mut atoms: Vec<StructureAtom> = Vec::new();
    let mut atom_index: HashMap<usize, usize> = HashMap::new();
    let mut atoms_by_serial: HashMap<usize, usize> = HashMap::new();
    let mut seen_names: HashSet<(usize, String)> = HashSet::new();
    for hierarchy in pdb.atoms_with_hierarchy() {
        let residue = hierarchy.residue();
        let residue_key = std::ptr::from_ref(residue) as usize;
        let conformer = hierarchy.conformer();
        let alternative = *first_alternative.entry(residue_key).or_insert_with(|| {
            residue
                .conformers()
                .find(|conformer| conformer.alternative_location().is_some())
                .map(|conformer| std::ptr::from_ref(conformer) as usize)
        });
        if conformer.alternative_location().is_some()
            && Some(std::ptr::from_ref(conformer) as usize) != alternative
        {
            continue;
        }
        let next = residues.len();
        let index = *residue_index.entry(residue_key).or_insert_with(|| {
            let parsed_name = residue.name().unwrap_or("UNK");
            let (sequence, insertion_code) = residue.id();
            residues.push(ResidueMeta {
                name: raw_names
                    .get(&(sequence, parsed_name.to_string()))
                    .cloned()
                    .unwrap_or_else(|| parsed_name.to_string()),
                sequence,
                insertion_code: insertion_code.map(str::to_string),
                chain: hierarchy.chain().id().to_string(),
                hetero: hierarchy.atom().hetero(),
            });
            next
        });
        let atom = hierarchy.atom();
        if !seen_names.insert((index, atom.name().to_string())) {
            continue;
        }
        let element = atom
            .element()
            .map(|element| element.symbol().to_ascii_uppercase())
            .unwrap_or_else(|| element_from_name(atom.name()));
        atom_index.insert(std::ptr::from_ref(atom) as usize, atoms.len());
        atoms_by_serial
            .entry(atom.serial_number())
            .or_insert(atoms.len());
        atoms.push(StructureAtom {
            residue: index,
            name: atom.name().to_string(),
            element,
            position: [atom.x(), atom.y(), atom.z()],
        });
    }

    let mut explicit = Vec::new();
    for (first, second, _) in pdb.bonds() {
        if let (Some(first), Some(second)) = (
            atom_index.get(&(std::ptr::from_ref(first) as usize)),
            atom_index.get(&(std::ptr::from_ref(second) as usize)),
        ) {
            explicit.push((*first, *second));
        }
    }
    for &(first, second) in raw_bonds {
        if let (Some(first), Some(second)) =
            (atoms_by_serial.get(&first), atoms_by_serial.get(&second))
        {
            explicit.push((*first, *second));
        }
    }
    explicit.sort_unstable();
    explicit.dedup();
    let structure = Structure::new(atoms, residues.len(), &explicit);

    // Monosaccharides: every ring the coordinates show, whatever the residue
    // is called; named sugars whose ring atoms are missing, from the name.
    let mut units: Vec<Unit> = Vec::new();
    for (index, residue) in residues.iter().enumerate() {
        let hint = is_sugar_name(&residue.name, residue.hetero)
            .then(|| name_reading(&residue.name))
            .flatten();
        let position_hint = |atom: usize| atom_position(&structure.atoms[atom].name, 'C');
        let rings = geometry::residue_units(
            &structure,
            index,
            if hint.is_some() { 2 } else { 3 },
            &position_hint,
        );
        if rings.is_empty() {
            if hint.is_some() && is_table_name(&residue.name) {
                units.push(Unit {
                    residue: index,
                    ring: None,
                    name: hint,
                });
            }
            continue;
        }
        // A residue with a sugar ring is no amino acid: a GLYCAM code in an
        // ATOM record (`QYS`, `UYB`) is a hint too.
        let hint = hint.or_else(|| name_reading(&residue.name));
        // A name describes one monosaccharide, not a component of several.
        let single = rings.len() == 1;
        for ring in rings {
            units.push(Unit {
                residue: index,
                ring: Some(ring),
                name: if single { hint.clone() } else { None },
            });
        }
    }
    if units.is_empty() {
        return Ok(PdbReading::default());
    }

    let mut backbone_index: BackboneIndex = HashMap::new();
    for (unit_index, unit) in units.iter().enumerate() {
        if let Some(ring) = &unit.ring {
            for (offset, carbon) in ring.backbone.iter().enumerate() {
                backbone_index.insert(*carbon, (unit_index, offset as u8 + 1));
            }
        }
    }
    let unit_of_residue = {
        let mut map: HashMap<usize, usize> = HashMap::new();
        for (unit_index, unit) in units.iter().enumerate() {
            map.entry(unit.residue).or_insert(unit_index);
        }
        map
    };
    let name_only_residues = units
        .iter()
        .filter(|unit| unit.ring.is_none())
        .map(|unit| unit.residue)
        .collect::<HashSet<_>>();

    // Read every ring unit: carbons, substituents, linkage, aglycone.
    let mut coordinate_readings: HashMap<usize, CoordinateReading> = HashMap::new();
    let mut links: Vec<Link> = Vec::new();
    // The attachment site is the residue an explicit record (CONECT, LINK)
    // bonds the glycan to, as before; a contact perceived from distances
    // never names it.
    let explicit_bonds = explicit
        .iter()
        .map(|&(first, second)| (first.min(second), first.max(second)))
        .collect::<HashSet<_>>();
    let is_explicit = |first: usize, second: usize| {
        explicit_bonds.contains(&(first.min(second), first.max(second)))
    };
    let mut aglycone: HashMap<usize, usize> = HashMap::new();
    let mut substituent_residues: HashSet<usize> = HashSet::new();
    let mut seen_bridges: HashSet<usize> = HashSet::new();
    for (unit_index, unit) in units.iter().enumerate() {
        let Some(ring) = &unit.ring else {
            continue;
        };
        let carbons = geometry::carbon_readings(&structure, ring);
        let mut substituents = Vec::new();
        let mut unrecognized = Vec::new();
        let other_sugar = |atom: usize| {
            backbone_index
                .get(&atom)
                .is_some_and(|(owner, _)| *owner != unit_index)
        };
        for index in 0..ring.backbone.len() {
            let position = index as u8 + 1;
            let carbon = ring.backbone[index];
            for atom in ring.exocyclic(&structure, index) {
                if !matches!(
                    structure.atoms[atom].element.as_str(),
                    "O" | "N" | "S" | "F"
                ) {
                    continue;
                }
                let group =
                    geometry::classify_group(&structure, unit.residue, carbon, atom, &other_sugar);
                let atom_residue = structure.atoms[atom].residue;
                let element = structure.atoms[atom].element.clone();
                if index == ring.anomeric {
                    match group {
                        Group::Sugar { carbon: target } => {
                            let (owner, target_position) = backbone_index[&target];
                            let owner_anomeric = units[owner]
                                .ring
                                .as_ref()
                                .is_some_and(|other| other.anomeric_carbon() == target);
                            if owner_anomeric {
                                // Two anomeric carbons share the oxygen (a 1-1
                                // or 1-2 link): the oxygen's residue is the
                                // parent, as the PDB records it.
                                if !seen_bridges.insert(atom) {
                                    continue;
                                }
                                let parent_is_owner = atom_residue == units[owner].residue
                                    || (atom_residue != unit.residue && owner < unit_index);
                                links.push(if parent_is_owner {
                                    Link {
                                        parent: owner,
                                        child: unit_index,
                                        parent_position: target_position,
                                        child_position: position,
                                    }
                                } else {
                                    Link {
                                        parent: unit_index,
                                        child: owner,
                                        parent_position: position,
                                        child_position: target_position,
                                    }
                                });
                            } else {
                                links.push(Link {
                                    parent: owner,
                                    child: unit_index,
                                    parent_position: target_position,
                                    child_position: position,
                                });
                            }
                        }
                        Group::Substituent {
                            descriptor,
                            external: false,
                            ..
                        } if atom_residue == unit.residue => {
                            substituents.push(CoordinateSubstituent {
                                position,
                                descriptor: descriptor.to_string(),
                                external: false,
                                bare: false,
                            });
                        }
                        Group::Bare | Group::Substituent { .. } => {
                            if atom_residue != unit.residue && is_explicit(carbon, atom) {
                                aglycone.entry(unit_index).or_insert(atom_residue);
                            }
                        }
                        Group::Other { atom: next } => {
                            let (residue, bonded) = if atom_residue != unit.residue {
                                (atom_residue, is_explicit(carbon, atom))
                            } else {
                                (structure.atoms[next].residue, is_explicit(atom, next))
                            };
                            if bonded && !name_only_residues.contains(&residue) {
                                aglycone.entry(unit_index).or_insert(residue);
                            }
                        }
                    }
                    continue;
                }
                match group {
                    Group::Bare => {
                        if element != "O" {
                            substituents.push(CoordinateSubstituent {
                                position,
                                descriptor: element,
                                external: atom_residue != unit.residue,
                                bare: true,
                            });
                        } else {
                            substituents.push(CoordinateSubstituent {
                                position,
                                descriptor: "O".into(),
                                external: false,
                                bare: true,
                            });
                        }
                    }
                    Group::Substituent {
                        descriptor,
                        external,
                        residue,
                    } => {
                        let external = external || atom_residue != unit.residue;
                        if external {
                            substituent_residues.insert(residue);
                            substituent_residues.insert(atom_residue);
                        }
                        substituents.push(CoordinateSubstituent {
                            position,
                            descriptor: descriptor.to_string(),
                            external,
                            bare: false,
                        });
                    }
                    Group::Sugar { carbon: target } => {
                        // The glycosidic oxygen of a child: recorded from the
                        // child's side. Anything else is not a linkage.
                        let (owner, _) = backbone_index[&target];
                        let child_anomeric = units[owner]
                            .ring
                            .as_ref()
                            .is_some_and(|other| other.anomeric_carbon() == target);
                        if !child_anomeric {
                            unrecognized.push(position);
                        }
                    }
                    Group::Other { atom: next } => {
                        if !name_only_residues.contains(&structure.atoms[next].residue) {
                            unrecognized.push(position);
                        }
                    }
                }
            }
        }
        unrecognized.sort_unstable();
        unrecognized.dedup();
        coordinate_readings.insert(
            unit_index,
            CoordinateReading {
                complete: structure.residue_has_hydrogens(unit.residue),
                anomeric_position: ring.anomeric_position(),
                ring_end: ring.ring_end_position(),
                carbons,
                substituents,
                unrecognized,
            },
        );
    }

    // Residues named as sugars whose ring atoms are missing: atom names
    // decide their linkages and substituents, as the name decides the rest.
    let name_only = NameOnly {
        structure: &structure,
        residues: &residues,
        units: &units,
        unit_of_residue: &unit_of_residue,
        name_only_residues: &name_only_residues,
        explicit: &explicit,
    };
    links.extend(name_only.links());
    let name_substituents = name_only.substituents();
    // Residues bonded by an explicit record to a glycan residue, other than
    // glycan residues and its substituents.
    let mut explicit_attachments: HashMap<usize, usize> = HashMap::new();
    for &(first, second) in &explicit {
        let (first_residue, second_residue) = (
            structure.atoms[first].residue,
            structure.atoms[second].residue,
        );
        for (sugar, other) in [
            (first_residue, second_residue),
            (second_residue, first_residue),
        ] {
            if unit_of_residue.contains_key(&sugar)
                && !unit_of_residue.contains_key(&other)
                && !substituent_residues.contains(&other)
                && !matches!(residues[other].name.as_str(), "SO3" | "MEX" | "ACX" | "PCX")
            {
                explicit_attachments.entry(sugar).or_insert(other);
            }
        }
    }

    // A child has one parent.
    links.sort_unstable();
    links.dedup();
    let mut has_parent: HashSet<usize> = HashSet::new();
    links.retain(|link| link.parent != link.child && has_parent.insert(link.child));

    // Connected glycans, in file order.
    let mut adjacency: HashMap<usize, Vec<usize>> = HashMap::new();
    for link in &links {
        adjacency.entry(link.parent).or_default().push(link.child);
        adjacency.entry(link.child).or_default().push(link.parent);
    }
    let mut assigned = vec![false; units.len()];
    let mut components: Vec<Vec<usize>> = Vec::new();
    for start in 0..units.len() {
        if assigned[start] {
            continue;
        }
        let mut stack = vec![start];
        let mut component = Vec::new();
        while let Some(current) = stack.pop() {
            if assigned[current] {
                continue;
            }
            assigned[current] = true;
            component.push(current);
            stack.extend(adjacency.get(&current).into_iter().flatten().copied());
        }
        component.sort_unstable();
        components.push(component);
    }

    let mut reading = PdbReading::default();
    let mut contradictions = Vec::new();
    for (glycan_index, component) in components.iter().enumerate() {
        let members = component.iter().copied().collect::<HashSet<_>>();
        let mut graph = ResidueGraph::new();
        let mut nodes = HashMap::new();
        let mut provenance = Vec::with_capacity(component.len());
        for &unit_index in component {
            let unit = &units[unit_index];
            let meta = &residues[unit.residue];
            let (residue, notes, record) = residue_for_unit(
                &structure,
                unit,
                coordinate_readings.get(&unit_index),
                name_substituents.get(&unit_index),
                options.policy,
            );
            let node = graph.add_residue(residue);
            nodes.insert(unit_index, node);
            provenance.push(PdbResidueReference {
                node_index: node.index(),
                chain: meta.chain.clone(),
                sequence_number: meta.sequence,
                insertion_code: meta.insertion_code.clone(),
            });
            for note in notes {
                let warning = PdbWarning {
                    glycan_index,
                    node_index: node.index(),
                    chain: meta.chain.clone(),
                    sequence_number: meta.sequence,
                    insertion_code: meta.insertion_code.clone(),
                    residue_name: meta.name.clone(),
                    kind: note.kind,
                    named: record.named.clone(),
                    coordinates: record.coordinates.clone(),
                    positions: note.positions,
                    aspects: note.aspects,
                };
                if warning.kind == PdbWarningKind::NameContradictsCoordinates {
                    contradictions.push(warning.clone());
                }
                reading.warnings.push(warning);
            }
            reading.residues.push(ResidueReading {
                glycan_index,
                node_index: node.index(),
                chain: meta.chain.clone(),
                sequence_number: meta.sequence,
                insertion_code: meta.insertion_code.clone(),
                residue_name: meta.name.clone(),
                from_coordinates: record.from_coordinates,
                from_name: record.from_name,
                ring_atoms: record.ring_atoms,
                backbone_atoms: record.backbone_atoms,
                centres: record.centres,
            });
        }
        for link in links
            .iter()
            .filter(|link| members.contains(&link.parent) && members.contains(&link.child))
        {
            graph.add_linkage(
                nodes[&link.parent],
                nodes[&link.child],
                Linkage::new(
                    CarbonPosition(link.parent_position),
                    CarbonPosition(link.child_position),
                ),
            );
        }
        let root = component
            .iter()
            .copied()
            .find(|unit| !has_parent.contains(unit))
            .unwrap_or(component[0]);
        graph.set_root(nodes[&root]);
        let attachment_site = aglycone
            .get(&root)
            .map(|residue| residues[*residue].site())
            .or_else(|| {
                component
                    .iter()
                    .find_map(|unit| explicit_attachments.get(&units[*unit].residue))
                    .map(|residue| residues[*residue].site())
            })
            .or_else(|| Some(residues[units[root].residue].site()));
        reading.glycans.push(ExtractedGlycanWithProvenance {
            attachment_site,
            graph,
            residues: provenance,
        });
    }
    if options.policy == NamePolicy::Strict && !contradictions.is_empty() {
        return Err(PdbError::NamesContradictCoordinates(contradictions));
    }
    Ok(reading)
}

struct UnitNote {
    kind: PdbWarningKind,
    positions: Vec<u8>,
    aspects: Vec<String>,
}

#[derive(Default)]
struct UnitRecord {
    named: Option<String>,
    coordinates: Option<String>,
    from_coordinates: Option<String>,
    from_name: Option<String>,
    ring_atoms: Vec<String>,
    backbone_atoms: Vec<String>,
    centres: Vec<StereoCentre>,
}

fn residue_for_unit(
    structure: &Structure,
    unit: &Unit,
    coordinates: Option<&CoordinateReading>,
    name_substituents: Option<&Vec<(u8, String)>>,
    policy: NamePolicy,
) -> (Monosaccharide, Vec<UnitNote>, UnitRecord) {
    let mut record = UnitRecord {
        named: unit.name.as_ref().map(|name| residue_label(&name.residue)),
        from_name: unit.name.as_ref().map(|name| unique_residue(&name.residue)),
        ..UnitRecord::default()
    };
    let (Some(ring), Some(coordinates)) = (&unit.ring, coordinates) else {
        // Ring atoms missing: the name decides.
        let mut residue = unit
            .name
            .as_ref()
            .map(|name| name.residue.clone())
            .unwrap_or_else(|| {
                residue_from_kind(crabwurcs_core::ResidueKind::Hex)
                    .expect("the generic hexose registry entry is valid")
            });
        if let Some(substituents) = name_substituents {
            residue.modifications.extend(
                substituents
                    .iter()
                    .map(|(position, descriptor)| modification(*position, descriptor)),
            );
            residue.modifications.sort_by(|first, second| {
                (first.position.0, first.descriptor.as_str())
                    .cmp(&(second.position.0, second.descriptor.as_str()))
            });
            residue.modifications.dedup();
        }
        let note = UnitNote {
            kind: PdbWarningKind::NoRing,
            positions: vec![],
            aspects: vec![],
        };
        return (residue, vec![note], record);
    };
    let name = |atom: usize| structure.atoms[atom].name.clone();
    record.ring_atoms = ring.ring.iter().map(|atom| name(*atom)).collect();
    record.backbone_atoms = ring.backbone.iter().map(|atom| name(*atom)).collect();
    record.centres = coordinates
        .carbons
        .iter()
        .enumerate()
        .filter_map(|(index, carbon)| {
            let (code, volume) = match carbon {
                CarbonReading::Stereo { code, volume } => (*code, *volume),
                CarbonReading::Anomeric {
                    code: Some((code, volume)),
                } => (*code, *volume),
                _ => return None,
            };
            Some(StereoCentre {
                position: index as u8 + 1,
                atom: name(ring.backbone[index]),
                code,
                volume,
            })
        })
        .collect();
    let combined = reading::combine(unit.name.as_ref(), coordinates);
    record.from_coordinates = Some(unique_residue(&combined.coordinates_only));
    record.coordinates = Some(residue_label(&combined.coordinates_first));
    // Notes of one kind are reported together.
    let mut grouped: BTreeMap<u8, UnitNote> = BTreeMap::new();
    let order = |kind: PdbWarningKind| -> u8 {
        match kind {
            PdbWarningKind::NameContradictsCoordinates => 0,
            PdbWarningKind::MissingAtoms => 1,
            PdbWarningKind::DistortedGeometry => 2,
            PdbWarningKind::UnreadableGeometry => 3,
            PdbWarningKind::NoRing => 4,
            PdbWarningKind::UnrecognizedGroup => 5,
        }
    };
    for note in combined.notes {
        let entry = grouped.entry(order(note.kind)).or_insert_with(|| UnitNote {
            kind: note.kind,
            positions: vec![],
            aspects: vec![],
        });
        entry.positions.extend(note.positions);
        entry.aspects.push(note.aspect);
    }
    let notes = grouped
        .into_values()
        .map(|mut note| {
            note.positions.sort_unstable();
            note.positions.dedup();
            note
        })
        .collect();
    let residue = match policy {
        NamePolicy::Names => combined.names_first,
        NamePolicy::Coordinates | NamePolicy::Strict => combined.coordinates_first,
    };
    (residue, notes, record)
}

/// The atom-name reading of residues named as sugars whose ring atoms are
/// missing (the extractor's behaviour before coordinates decided).
struct NameOnly<'a> {
    structure: &'a Structure,
    residues: &'a [ResidueMeta],
    units: &'a [Unit],
    unit_of_residue: &'a HashMap<usize, usize>,
    name_only_residues: &'a HashSet<usize>,
    explicit: &'a [(usize, usize)],
}

impl NameOnly<'_> {
    fn involves(&self, first: usize, second: usize) -> bool {
        self.name_only_residues.contains(&first) || self.name_only_residues.contains(&second)
    }

    /// Glycosidic bonds with a name-only residue at one end, oriented by
    /// atom names (`O4`–`C1`): the explicit records of the file, and bonds
    /// perceived from distances that end at the child's declared anomeric
    /// carbon.
    fn links(&self) -> Vec<Link> {
        let atoms = &self.structure.atoms;
        let mut pairs = self
            .explicit
            .iter()
            .map(|&(first, second)| (first, second, true))
            .collect::<Vec<_>>();
        for &residue in self.name_only_residues {
            for &atom in &self.structure.residue_atoms[residue] {
                for neighbor in self.structure.heavy_neighbors(atom) {
                    pairs.push((atom, neighbor, false));
                }
            }
        }
        let mut links = Vec::new();
        for (first, second, explicit) in pairs {
            let (first_residue, second_residue) = (atoms[first].residue, atoms[second].residue);
            if first_residue == second_residue || !self.involves(first_residue, second_residue) {
                continue;
            }
            let (Some(&first_unit), Some(&second_unit)) = (
                self.unit_of_residue.get(&first_residue),
                self.unit_of_residue.get(&second_residue),
            ) else {
                continue;
            };
            let Some((parent, child, parent_position, child_position)) = orient_glycosidic_bond(
                first_unit,
                &atoms[first].name,
                second_unit,
                &atoms[second].name,
            ) else {
                continue;
            };
            let declared = self.units[child]
                .name
                .as_ref()
                .map(|name| name.residue.anomeric_position)
                .or_else(|| {
                    self.units[child]
                        .ring
                        .as_ref()
                        .map(RingUnit::anomeric_position)
                });
            if !explicit && declared.is_some_and(|position| position != child_position) {
                continue;
            }
            links.push(Link {
                parent,
                child,
                parent_position,
                child_position,
            });
        }
        links
    }

    /// Substituents of name-only residues: `SO3` sulfur close to a ring
    /// oxygen, and `MEX`/`ACX`/`PCX` residues bonded to one.
    fn substituents(&self) -> HashMap<usize, Vec<(u8, String)>> {
        let atoms = &self.structure.atoms;
        let mut found: HashMap<usize, Vec<(u8, String)>> = HashMap::new();
        if self.name_only_residues.is_empty() {
            return found;
        }
        let oxygen = |atom: &StructureAtom| {
            self.name_only_residues
                .contains(&atom.residue)
                .then(|| atom_position(&atom.name, 'O'))
                .flatten()
        };
        let sulfurs = atoms
            .iter()
            .filter(|atom| self.residues[atom.residue].name == "SO3" && atom.name.starts_with('S'))
            .collect::<Vec<_>>();
        for atom in atoms {
            let Some(position) = oxygen(atom) else {
                continue;
            };
            for sulfur in &sulfurs {
                let squared = (0..3)
                    .map(|axis| (atom.position[axis] - sulfur.position[axis]).powi(2))
                    .sum::<f64>();
                if (0.8f64.powi(2)..=1.9f64.powi(2)).contains(&squared) {
                    found
                        .entry(self.unit_of_residue[&atom.residue])
                        .or_default()
                        .push((position, "OSO/3=O/3=O".into()));
                }
            }
        }
        for &(first, second) in self.explicit {
            for (sugar, other) in [(first, second), (second, first)] {
                let Some(position) = oxygen(&atoms[sugar]) else {
                    continue;
                };
                let descriptor = match self.residues[atoms[other].residue].name.as_str() {
                    "MEX" => "OC",
                    "ACX" => "OCC/3=O",
                    "PCX" => "OP^XOCCNC/7C/7C/3O/3=O",
                    _ => continue,
                };
                found
                    .entry(self.unit_of_residue[&atoms[sugar].residue])
                    .or_default()
                    .push((position, descriptor.into()));
            }
        }
        found
    }
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
mod tests;
