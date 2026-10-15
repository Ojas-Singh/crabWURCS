//! Name-independent reading of monosaccharides from 3D coordinates.
//!
//! Nothing here reads an atom name or a residue name to decide chemistry:
//! bonds come from elements and distances (plus explicit `CONECT`/`LINK`
//! records), rings and the carbon numbering from the bond graph, and every
//! configuration from a signed volume. Atom names are used in exactly one
//! place, as a tie-break hint when the graph cannot tell the anomeric carbon
//! from the ring-closing carbon (a pentopyranose with its O1 missing).
//!
//! # Configuration
//!
//! The backbone is numbered from the ring graph: the ring carbon bonded to
//! the ring oxygen that carries an exocyclic heteroatom is the anomeric
//! carbon; C1 precedes it when it carries an exocyclic carbon (a ketose); the
//! numbering then walks the ring away from the ring oxygen and continues
//! along the exocyclic carbon chain of the ring-closing carbon.
//!
//! At backbone carbon `Ck` with one exocyclic heteroatom `X` (the ring oxygen
//! counts as the heteroatom of the ring-closing carbon), the signed volume
//!
//! ```text
//! V = (C(k-1) - Ck) . ((C(k+1) - Ck) x (X - Ck))
//! ```
//!
//! is negative when `X` lies on the right of the Fischer projection (WURCS
//! `2`) and positive on the left (WURCS `1`). It is a property of the
//! configuration, not of the conformation: ring flips keep its sign, an
//! inversion changes it. The anomeric carbon uses the ring oxygen in place of
//! `C(k-1)`: in the formal Fischer projection of the cyclic form the ring
//! bond takes the upper position of the anomeric carbon, for aldoses (C1)
//! and ketoses (C2, whose C1 then lies sideways) alike. The anomer is α when
//! the anomeric heteroatom and the heteroatom of the anomeric reference atom
//! are on the same side (IUPAC 2-Carb-6.2). wwPDB CCD references fix both
//! conventions (see the crate's tests).
//!
//! The volume is normalised by the three bond lengths, so it is about ±0.77
//! for an ideal tetrahedral centre and 0 for a planar one.

use std::collections::HashMap;

/// Below this normalised volume a centre is too close to planar to read.
pub(crate) const DISTORTED_VOLUME: f64 = 0.15;

/// Bonds are perceived when `d < BOND_FACTOR * (r_i + r_j)`.
const BOND_FACTOR: f64 = 1.2;
/// Explicit `CONECT`/`LINK` bonds longer than this are ignored.
const EXPLICIT_BOND_LIMIT: f64 = 3.0;
const MAXIMUM_BACKBONE: usize = 10;
/// A C–O bond shorter than this is a C=O.
const DOUBLE_CO: f64 = 1.30;
/// An S–O or P–O bond shorter than this is an S=O or P=O.
const DOUBLE_SO: f64 = 1.50;

#[derive(Debug, Clone)]
pub(crate) struct StructureAtom {
    pub residue: usize,
    pub name: String,
    /// Upper-case element symbol.
    pub element: String,
    pub position: [f64; 3],
}

impl StructureAtom {
    pub fn is_hydrogen(&self) -> bool {
        matches!(self.element.as_str(), "H" | "D")
    }

    pub fn is(&self, element: &str) -> bool {
        self.element == element
    }
}

/// The atoms of one model with their perceived bonds.
pub(crate) struct Structure {
    pub atoms: Vec<StructureAtom>,
    pub neighbors: Vec<Vec<usize>>,
    pub residue_atoms: Vec<Vec<usize>>,
}

fn covalent_radius(element: &str) -> Option<f64> {
    Some(match element {
        "H" | "D" => 0.31,
        "C" => 0.76,
        "N" => 0.71,
        "O" => 0.66,
        "S" => 1.05,
        "P" => 1.07,
        "F" => 0.57,
        "CL" => 1.02,
        "BR" => 1.20,
        "I" => 1.39,
        "SE" => 1.20,
        "B" => 0.84,
        "SI" => 1.11,
        _ => return None,
    })
}

fn distance_squared(left: [f64; 3], right: [f64; 3]) -> f64 {
    (0..3).map(|axis| (left[axis] - right[axis]).powi(2)).sum()
}

fn subtract(left: [f64; 3], right: [f64; 3]) -> [f64; 3] {
    [left[0] - right[0], left[1] - right[1], left[2] - right[2]]
}

fn cross(left: [f64; 3], right: [f64; 3]) -> [f64; 3] {
    [
        left[1] * right[2] - left[2] * right[1],
        left[2] * right[0] - left[0] * right[2],
        left[0] * right[1] - left[1] * right[0],
    ]
}

fn dot(left: [f64; 3], right: [f64; 3]) -> f64 {
    left[0] * right[0] + left[1] * right[1] + left[2] * right[2]
}

fn length(vector: [f64; 3]) -> f64 {
    dot(vector, vector).sqrt()
}

impl Structure {
    /// Perceive covalent bonds from distances, organic elements only (a metal
    /// ion near a sugar is never bonded to it).
    ///
    /// Each hydrogen is bonded to its nearest heavy atom. Heavy atoms are
    /// bonded in order of confidence, within their valence: bonds inside a
    /// residue first, then the explicit `CONECT`/`LINK` bonds of the file,
    /// then contacts between residues, shortest (relative to the covalent
    /// radii) first. A clash between residues (a sulfate packed onto the next
    /// ring in a crystal model) therefore cannot take the place of a ring
    /// bond.
    pub fn new(atoms: Vec<StructureAtom>, residues: usize, explicit: &[(usize, usize)]) -> Self {
        let mut neighbors = vec![Vec::new(); atoms.len()];
        let mut residue_atoms = vec![Vec::new(); residues];
        for (index, atom) in atoms.iter().enumerate() {
            residue_atoms[atom.residue].push(index);
        }
        let radii = atoms
            .iter()
            .map(|atom| covalent_radius(&atom.element))
            .collect::<Vec<_>>();
        let maximum_radius = radii.iter().flatten().copied().fold(0.0, f64::max);
        let cell = (BOND_FACTOR * 2.0 * maximum_radius).max(1.0);
        let key = |position: [f64; 3]| {
            (
                (position[0] / cell).floor() as i64,
                (position[1] / cell).floor() as i64,
                (position[2] / cell).floor() as i64,
            )
        };
        let mut grid: HashMap<(i64, i64, i64), Vec<usize>> = HashMap::new();
        for (index, atom) in atoms.iter().enumerate() {
            if radii[index].is_some() {
                grid.entry(key(atom.position)).or_default().push(index);
            }
        }
        let mut hydrogen_partners: Vec<Vec<(f64, usize)>> = vec![Vec::new(); atoms.len()];
        let mut inside = Vec::new();
        let mut between = Vec::new();
        for (index, atom) in atoms.iter().enumerate() {
            let Some(radius) = radii[index] else {
                continue;
            };
            let (x, y, z) = key(atom.position);
            for dx in -1..=1 {
                for dy in -1..=1 {
                    for dz in -1..=1 {
                        let Some(cell_atoms) = grid.get(&(x + dx, y + dy, z + dz)) else {
                            continue;
                        };
                        for &other in cell_atoms {
                            if other <= index {
                                continue;
                            }
                            let other_atom = &atoms[other];
                            if atom.is_hydrogen() && other_atom.is_hydrogen() {
                                continue;
                            }
                            let other_radius = radii[other].expect("gridded atoms have a radius");
                            let squared = distance_squared(atom.position, other_atom.position);
                            let limit = BOND_FACTOR * (radius + other_radius);
                            if squared < 0.4f64.powi(2) || squared > limit * limit {
                                continue;
                            }
                            if atom.is_hydrogen() {
                                hydrogen_partners[index].push((squared, other));
                            } else if other_atom.is_hydrogen() {
                                hydrogen_partners[other].push((squared, index));
                            } else {
                                let ratio = squared.sqrt() / (radius + other_radius);
                                if atom.residue == other_atom.residue {
                                    inside.push((ratio, index, other));
                                } else {
                                    between.push((ratio, index, other));
                                }
                            }
                        }
                    }
                }
            }
        }
        for (hydrogen, mut partners) in hydrogen_partners.into_iter().enumerate() {
            partners.sort_by(|left, right| left.0.total_cmp(&right.0));
            if let Some(&(_, partner)) = partners.first() {
                neighbors[hydrogen].push(partner);
                neighbors[partner].push(hydrogen);
            }
        }
        // Heavy neighbours an atom can take; a carbon counts its hydrogens.
        // A carbon with a C=O (carboxyl, carbonyl) takes three, an oxygen of
        // an S=O or P=O one.
        let free = |neighbors: &[Vec<usize>], atom: usize| {
            let element = atoms[atom].element.as_str();
            let bonded_within = |elements: &[&str], limit: f64| {
                neighbors[atom].iter().any(|neighbor| {
                    elements.contains(&atoms[*neighbor].element.as_str())
                        && distance_squared(atoms[atom].position, atoms[*neighbor].position)
                            < limit * limit
                })
            };
            let capacity = match element {
                "C" if bonded_within(&["O"], DOUBLE_CO) => 3,
                "O" if bonded_within(&["S", "P"], DOUBLE_SO) => 1,
                "O" | "SE" => 2,
                "F" | "CL" | "BR" | "I" => 1,
                _ => 4,
            };
            let used = neighbors[atom]
                .iter()
                .filter(|neighbor| element == "C" || !atoms[**neighbor].is_hydrogen())
                .count();
            used < capacity
        };
        let by_ratio =
            |left: &(f64, usize, usize), right: &(f64, usize, usize)| left.0.total_cmp(&right.0);
        inside.sort_by(by_ratio);
        between.sort_by(by_ratio);
        for (_, first, second) in inside {
            if free(&neighbors, first) && free(&neighbors, second) {
                neighbors[first].push(second);
                neighbors[second].push(first);
            }
        }
        for &(first, second) in explicit {
            if first == second
                || first >= atoms.len()
                || second >= atoms.len()
                || radii[first].is_none()
                || radii[second].is_none()
                || atoms[first].is_hydrogen()
                || atoms[second].is_hydrogen()
                || neighbors[first].contains(&second)
                || distance_squared(atoms[first].position, atoms[second].position)
                    > EXPLICIT_BOND_LIMIT.powi(2)
            {
                continue;
            }
            neighbors[first].push(second);
            neighbors[second].push(first);
        }
        for (_, first, second) in between {
            if !neighbors[first].contains(&second)
                && free(&neighbors, first)
                && free(&neighbors, second)
            {
                neighbors[first].push(second);
                neighbors[second].push(first);
            }
        }
        for list in &mut neighbors {
            list.sort_unstable();
            list.dedup();
        }
        Self {
            atoms,
            neighbors,
            residue_atoms,
        }
    }

    pub fn heavy_neighbors(&self, atom: usize) -> impl Iterator<Item = usize> + '_ {
        self.neighbors[atom]
            .iter()
            .copied()
            .filter(|neighbor| !self.atoms[*neighbor].is_hydrogen())
    }

    pub fn hydrogen_count(&self, atom: usize) -> usize {
        self.neighbors[atom]
            .iter()
            .filter(|neighbor| self.atoms[**neighbor].is_hydrogen())
            .count()
    }

    pub fn residue_has_hydrogens(&self, residue: usize) -> bool {
        self.residue_atoms[residue]
            .iter()
            .any(|atom| self.atoms[*atom].is_hydrogen())
    }

    fn is_heteroatom(&self, atom: usize) -> bool {
        matches!(self.atoms[atom].element.as_str(), "O" | "N" | "S" | "F")
    }

    /// Normalised signed volume at `centre` for the ordered substituents
    /// `(top, bottom, side)`.
    pub fn normalised_volume(&self, centre: usize, top: usize, bottom: usize, side: usize) -> f64 {
        let origin = self.atoms[centre].position;
        let a = subtract(self.atoms[top].position, origin);
        let b = subtract(self.atoms[bottom].position, origin);
        let c = subtract(self.atoms[side].position, origin);
        let scale = length(a) * length(b) * length(c);
        if scale < 1e-9 {
            return 0.0;
        }
        dot(a, cross(b, c)) / scale
    }
}

/// Fischer side of a stereocentre from its normalised volume: `'2'` (right)
/// for a negative volume, `'1'` (left) for a positive one.
pub(crate) fn fischer_code(volume: f64) -> char {
    if volume < 0.0 { '2' } else { '1' }
}

/// One monosaccharide found in the coordinates: its ring and its numbered
/// carbon backbone.
#[derive(Debug, Clone)]
pub(crate) struct RingUnit {
    /// Ring oxygen first, then the anomeric carbon, then around the ring.
    pub ring: Vec<usize>,
    /// C1..Cn.
    pub backbone: Vec<usize>,
    /// Index of the anomeric carbon in `backbone` (0 for an aldose).
    pub anomeric: usize,
    /// Index of the ring-closing carbon in `backbone`.
    pub ring_end: usize,
}

impl RingUnit {
    pub fn ring_oxygen(&self) -> usize {
        self.ring[0]
    }

    pub fn anomeric_carbon(&self) -> usize {
        self.backbone[self.anomeric]
    }

    pub fn anomeric_position(&self) -> u8 {
        self.anomeric as u8 + 1
    }

    pub fn ring_end_position(&self) -> u8 {
        self.ring_end as u8 + 1
    }

    /// Exocyclic heavy atoms of a backbone carbon: neighbours that are neither
    /// backbone carbons of this unit nor its ring oxygen.
    pub fn exocyclic(&self, structure: &Structure, index: usize) -> Vec<usize> {
        let carbon = self.backbone[index];
        structure
            .heavy_neighbors(carbon)
            .filter(|neighbor| !self.backbone.contains(neighbor) && *neighbor != self.ring_oxygen())
            .collect()
    }
}

/// 5- and 6-membered rings of one residue made of one oxygen and carbons,
/// each as the ring oxygen followed by the shortest carbon path between its
/// two carbon neighbours.
fn residue_rings(structure: &Structure, residue: usize) -> Vec<Vec<usize>> {
    let in_residue = |atom: usize| structure.atoms[atom].residue == residue;
    let mut rings: Vec<Vec<usize>> = Vec::new();
    for &oxygen in &structure.residue_atoms[residue] {
        if !structure.atoms[oxygen].is("O") {
            continue;
        }
        let heavy = structure.heavy_neighbors(oxygen).collect::<Vec<_>>();
        if heavy.len() != 2
            || !heavy
                .iter()
                .all(|atom| in_residue(*atom) && structure.atoms[*atom].is("C"))
        {
            continue;
        }
        let (start, end) = (heavy[0], heavy[1]);
        let mut previous: HashMap<usize, Option<usize>> = HashMap::from([(start, None)]);
        let mut frontier = vec![start];
        let mut depth = 1;
        while !frontier.is_empty() && !previous.contains_key(&end) && depth < 5 {
            let mut next = Vec::new();
            for atom in frontier {
                for neighbor in structure.heavy_neighbors(atom) {
                    if previous.contains_key(&neighbor)
                        || neighbor == oxygen
                        || !in_residue(neighbor)
                        || !structure.atoms[neighbor].is("C")
                    {
                        continue;
                    }
                    previous.insert(neighbor, Some(atom));
                    next.push(neighbor);
                }
            }
            frontier = next;
            depth += 1;
        }
        if !previous.contains_key(&end) {
            continue;
        }
        let mut path = vec![end];
        while let Some(Some(before)) = previous.get(path.last().expect("non-empty path")) {
            path.push(*before);
        }
        path.reverse();
        if !(4..=5).contains(&path.len()) {
            continue;
        }
        let mut ring = vec![oxygen];
        ring.extend(path);
        let mut sorted = ring.clone();
        sorted.sort_unstable();
        if rings.iter().any(|known| {
            let mut other = known.clone();
            other.sort_unstable();
            other == sorted
        }) {
            continue;
        }
        rings.push(ring);
    }
    rings
}

/// The monosaccharide rings of `residue`, numbered.
///
/// `minimum_decorated` is the number of ring carbons that must carry an
/// exocyclic heteroatom (ordinary cyclic ethers have none). `position_hint`
/// maps an atom to the carbon number its name declares (`C1` → 1); it is
/// consulted only when the graph cannot identify the anomeric carbon.
pub(crate) fn residue_units(
    structure: &Structure,
    residue: usize,
    minimum_decorated: usize,
    position_hint: &dyn Fn(usize) -> Option<u8>,
) -> Vec<RingUnit> {
    let mut units: Vec<RingUnit> = Vec::new();
    for ring in residue_rings(structure, residue) {
        if units
            .iter()
            .any(|unit| unit.ring.iter().any(|atom| ring.contains(atom)))
        {
            continue;
        }
        if let Some(unit) = number_ring(structure, residue, &ring, minimum_decorated, position_hint)
        {
            units.push(unit);
        }
    }
    units
}

fn number_ring(
    structure: &Structure,
    residue: usize,
    ring: &[usize],
    minimum_decorated: usize,
    position_hint: &dyn Fn(usize) -> Option<u8>,
) -> Option<RingUnit> {
    let oxygen = ring[0];
    let carbons = &ring[1..];
    let exocyclic = |carbon: usize| {
        structure
            .heavy_neighbors(carbon)
            .filter(|neighbor| !ring.contains(neighbor))
            .collect::<Vec<_>>()
    };
    let has_heteroatom = |carbon: usize| {
        exocyclic(carbon)
            .iter()
            .any(|atom| structure.is_heteroatom(*atom))
    };
    let has_residue_carbon = |carbon: usize| {
        exocyclic(carbon)
            .iter()
            .any(|atom| structure.atoms[*atom].is("C") && structure.atoms[*atom].residue == residue)
    };
    let decorated = carbons
        .iter()
        .filter(|carbon| has_heteroatom(**carbon))
        .count();
    if decorated < minimum_decorated {
        return None;
    }
    let (first, last) = (carbons[0], carbons[carbons.len() - 1]);
    let anomeric_first = match (has_heteroatom(first), has_heteroatom(last)) {
        (true, false) => true,
        (false, true) => false,
        (false, false) if has_residue_carbon(first) != has_residue_carbon(last) => {
            !has_residue_carbon(first)
        }
        _ => match (position_hint(first), position_hint(last)) {
            (Some(left), Some(right)) if left != right => left < right,
            _ => return None,
        },
    };
    let mut ordered = vec![oxygen];
    if anomeric_first {
        ordered.extend(carbons.iter().copied());
    } else {
        ordered.extend(carbons.iter().rev().copied());
    }
    let anomeric_carbon = ordered[1];

    // A nucleoside: the anomeric carbon carries a ring nitrogen (a base).
    if exocyclic(anomeric_carbon).iter().any(|atom| {
        structure.atoms[*atom].is("N")
            && structure
                .heavy_neighbors(*atom)
                .filter(|neighbor| *neighbor != anomeric_carbon)
                .count()
                >= 2
    }) {
        return None;
    }

    let ketose_c1 = exocyclic(anomeric_carbon)
        .into_iter()
        .filter(|atom| structure.atoms[*atom].is("C") && structure.atoms[*atom].residue == residue)
        .collect::<Vec<_>>();
    if ketose_c1.len() > 1 {
        return None;
    }
    let mut backbone = Vec::with_capacity(MAXIMUM_BACKBONE);
    backbone.extend(ketose_c1.iter().copied());
    let anomeric = backbone.len();
    backbone.extend(ordered[1..].iter().copied());
    let ring_end = backbone.len() - 1;
    while backbone.len() < MAXIMUM_BACKBONE {
        let current = *backbone.last().expect("non-empty backbone");
        let next = structure
            .heavy_neighbors(current)
            .filter(|atom| {
                structure.atoms[*atom].is("C")
                    && structure.atoms[*atom].residue == residue
                    && !backbone.contains(atom)
                    && !ring.contains(atom)
            })
            .collect::<Vec<_>>();
        if next.len() != 1 {
            break;
        }
        backbone.push(next[0]);
    }
    Some(RingUnit {
        ring: ordered,
        backbone,
        anomeric,
        ring_end,
    })
}

/// What the coordinates say about one backbone carbon.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) enum CarbonReading {
    /// The anomeric carbon; `code` is its Fischer side when its exocyclic
    /// heteroatom is present.
    Anomeric { code: Option<(char, f64)> },
    /// A stereocentre with its Fischer side and normalised volume.
    Stereo { code: char, volume: f64 },
    /// A chain carbon without an exocyclic heteroatom: deoxy, or atoms
    /// missing. `hydrogens` counts the bonded hydrogens of the file.
    NoHeteroatom { hydrogens: usize },
    /// A terminal carbon: `h` (CH2-O), `m` (CH3) or `A` (carboxyl).
    Terminal { code: char, hydrogens: usize },
    /// A valence pattern this reader does not interpret (a C-branch).
    Unreadable,
}

impl CarbonReading {
    pub fn code(self) -> char {
        match self {
            Self::Anomeric { .. } => 'a',
            Self::Stereo { code, .. } => code,
            Self::NoHeteroatom { .. } => 'd',
            Self::Terminal { code, .. } => code,
            Self::Unreadable => 'x',
        }
    }
}

pub(crate) fn carbon_readings(structure: &Structure, unit: &RingUnit) -> Vec<CarbonReading> {
    let count = unit.backbone.len();
    (0..count)
        .map(|index| {
            let carbon = unit.backbone[index];
            let exocyclic = unit.exocyclic(structure, index);
            let heteroatoms = exocyclic
                .iter()
                .copied()
                .filter(|atom| structure.is_heteroatom(*atom))
                .collect::<Vec<_>>();
            if index == unit.anomeric {
                // The ring oxygen takes the upper position of the anomeric
                // carbon in the formal Fischer projection, for aldoses and
                // ketoses alike.
                let (top, bottom) = (unit.ring_oxygen(), unit.backbone[index + 1]);
                let code = (heteroatoms.len() == 1).then(|| {
                    let volume = structure.normalised_volume(carbon, top, bottom, heteroatoms[0]);
                    (fischer_code(volume), volume)
                });
                return CarbonReading::Anomeric { code };
            }
            let terminal = index == 0 || index + 1 == count;
            if terminal {
                // The ring oxygen is a substituent of a terminal ring carbon
                // (C5 of a pentopyranose).
                let mut oxygens = exocyclic
                    .iter()
                    .filter(|atom| structure.atoms[**atom].is("O"))
                    .count();
                if index == unit.ring_end {
                    oxygens += 1;
                }
                if exocyclic.iter().any(|atom| !structure.atoms[*atom].is("O")) {
                    return CarbonReading::Unreadable;
                }
                let code = match oxygens {
                    0 => 'm',
                    1 => 'h',
                    2 => 'A',
                    _ => return CarbonReading::Unreadable,
                };
                return CarbonReading::Terminal {
                    code,
                    hydrogens: structure.hydrogen_count(carbon),
                };
            }
            let (top, bottom) = (unit.backbone[index - 1], unit.backbone[index + 1]);
            let mut side = exocyclic.clone();
            if index == unit.ring_end {
                side.push(unit.ring_oxygen());
            }
            match side.as_slice() {
                [] => CarbonReading::NoHeteroatom {
                    hydrogens: structure.hydrogen_count(carbon),
                },
                [atom] if structure.is_heteroatom(*atom) => {
                    let volume = structure.normalised_volume(carbon, top, bottom, *atom);
                    CarbonReading::Stereo {
                        code: fischer_code(volume),
                        volume,
                    }
                }
                _ => CarbonReading::Unreadable,
            }
        })
        .collect()
}

/// What an exocyclic heteroatom of a backbone carbon carries.
#[derive(Debug, Clone, PartialEq)]
pub(crate) enum Group {
    /// Nothing beyond the heteroatom (OH, NH2, a missing group).
    Bare,
    /// A recognised substituent. `external` when its atoms beyond the
    /// heteroatom belong to another residue (GLYCAM's `SO3`, `MEX`, `ACX`).
    Substituent {
        descriptor: &'static str,
        external: bool,
        /// The residue of its first atom beyond the heteroatom.
        residue: usize,
    },
    /// Bonded to `carbon`, a backbone carbon of another monosaccharide.
    Sugar { carbon: usize },
    /// Anything else, starting at `atom` (an aglycone, a protein residue, an
    /// unrecognised group).
    Other { atom: usize },
}

fn is_bare(structure: &Structure, atom: usize, from: usize) -> bool {
    structure
        .heavy_neighbors(atom)
        .all(|neighbor| neighbor == from)
}

fn others(structure: &Structure, atom: usize, from: usize) -> Vec<usize> {
    structure
        .heavy_neighbors(atom)
        .filter(|neighbor| *neighbor != from)
        .collect()
}

fn is_acid_centre(structure: &Structure, centre: usize, from: usize) -> bool {
    others(structure, centre, from)
        .iter()
        .all(|oxygen| structure.atoms[*oxygen].is("O") && is_bare(structure, *oxygen, centre))
}

fn is_phosphocholine(structure: &Structure, phosphorus: usize, from: usize) -> bool {
    let mut bridged = false;
    for oxygen in others(structure, phosphorus, from) {
        if !structure.atoms[oxygen].is("O") {
            return false;
        }
        let beyond = others(structure, oxygen, phosphorus);
        if beyond.is_empty() {
            continue;
        }
        let chain = (|| {
            let [first] = beyond.as_slice() else {
                return None;
            };
            let [second] = others(structure, *first, oxygen)[..] else {
                return None;
            };
            let [nitrogen] = others(structure, second, *first)[..] else {
                return None;
            };
            (structure.atoms[*first].is("C")
                && structure.atoms[second].is("C")
                && structure.atoms[nitrogen].is("N"))
            .then_some(nitrogen)
            .map(|nitrogen| (nitrogen, second))
        })();
        let Some((nitrogen, second)) = chain else {
            return false;
        };
        let methyls = others(structure, nitrogen, second);
        if bridged
            || methyls.len() != 3
            || !methyls.iter().all(|carbon| {
                structure.atoms[*carbon].is("C") && is_bare(structure, *carbon, nitrogen)
            })
        {
            return false;
        }
        bridged = true;
    }
    bridged
}

/// Classify what heteroatom `atom`, bonded to backbone carbon `carbon` of a
/// unit in `residue`, carries. `sugar_carbon` says whether an atom is a
/// backbone carbon of another monosaccharide.
pub(crate) fn classify_group(
    structure: &Structure,
    residue: usize,
    carbon: usize,
    atom: usize,
    sugar_carbon: &dyn Fn(usize) -> bool,
) -> Group {
    let beyond = others(structure, atom, carbon);
    let [next] = beyond[..] else {
        return if beyond.is_empty() {
            match structure.atoms[atom].element.as_str() {
                "F" => Group::Substituent {
                    descriptor: "F",
                    external: structure.atoms[atom].residue != residue,
                    residue: structure.atoms[atom].residue,
                },
                _ => Group::Bare,
            }
        } else {
            Group::Other { atom: beyond[0] }
        };
    };
    if sugar_carbon(next) {
        return Group::Sugar { carbon: next };
    }
    let external = structure.atoms[next].residue != residue;
    let substituent = |descriptor| Group::Substituent {
        descriptor,
        external,
        residue: structure.atoms[next].residue,
    };
    let element = structure.atoms[next].element.as_str();
    match structure.atoms[atom].element.as_str() {
        "O" => match element {
            "S" if is_acid_centre(structure, next, atom) => substituent("OSO/3=O/3=O"),
            "P" if is_acid_centre(structure, next, atom) => substituent("OPO/3O/3=O"),
            "P" if is_phosphocholine(structure, next, atom) => {
                substituent("OP^XOCCNC/7C/7C/3O/3=O")
            }
            "C" => {
                let rest = others(structure, next, atom);
                if rest.is_empty() {
                    return substituent("OC");
                }
                if is_acetyl(structure, next, atom) {
                    return substituent("OCC/3=O");
                }
                Group::Other { atom: next }
            }
            _ => Group::Other { atom: next },
        },
        "N" => match element {
            "S" if is_acid_centre(structure, next, atom) => substituent("NSO/3=O/3=O"),
            "C" if is_acetyl(structure, next, atom) => substituent("NCC/3=O"),
            "C" if is_glycolyl(structure, next, atom) => substituent("NCCO/3=O"),
            _ => Group::Other { atom: next },
        },
        _ => Group::Other { atom: next },
    }
}

/// `carbonyl` is the C of an acetyl group `C(=O)CH3` reached from `from`.
fn is_acetyl(structure: &Structure, carbonyl: usize, from: usize) -> bool {
    let rest = others(structure, carbonyl, from);
    let oxygens = rest
        .iter()
        .filter(|atom| structure.atoms[**atom].is("O") && is_bare(structure, **atom, carbonyl))
        .count();
    let methyls = rest
        .iter()
        .filter(|atom| structure.atoms[**atom].is("C") && is_bare(structure, **atom, carbonyl))
        .count();
    rest.len() == 2 && oxygens == 1 && methyls == 1
}

/// `carbonyl` is the C of a glycolyl group `C(=O)CH2OH` reached from `from`.
fn is_glycolyl(structure: &Structure, carbonyl: usize, from: usize) -> bool {
    let rest = others(structure, carbonyl, from);
    let [first, second] = rest[..] else {
        return false;
    };
    let (oxygen, methylene) = if structure.atoms[first].is("O") {
        (first, second)
    } else {
        (second, first)
    };
    if !structure.atoms[oxygen].is("O")
        || !is_bare(structure, oxygen, carbonyl)
        || !structure.atoms[methylene].is("C")
    {
        return false;
    }
    let [hydroxyl] = others(structure, methylene, carbonyl)[..] else {
        return false;
    };
    structure.atoms[hydroxyl].is("O") && is_bare(structure, hydroxyl, methylene)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn atom(residue: usize, name: &str, element: &str, position: [f64; 3]) -> StructureAtom {
        StructureAtom {
            residue,
            name: name.into(),
            element: element.into(),
            position,
        }
    }

    #[test]
    fn metal_ions_are_never_bonded_and_hydrogens_keep_one_partner() {
        let structure = Structure::new(
            vec![
                atom(0, "O1", "O", [0.0, 0.0, 0.0]),
                atom(1, "CA", "CA", [0.0, 0.0, 2.3]),
                atom(0, "H1", "H", [0.0, 0.96, 0.0]),
                atom(0, "O2", "O", [0.0, 2.0, 0.0]),
            ],
            2,
            &[],
        );
        assert!(structure.neighbors[1].is_empty());
        assert_eq!(structure.neighbors[2], vec![0]);
        assert!(!structure.neighbors[0].contains(&3));
    }

    #[test]
    fn a_clash_between_residues_does_not_displace_a_ring_bond() {
        // A sulfate sulfur of residue 1 packed 1.6 Å from the ring oxygen of
        // residue 0, which already has its two carbons.
        let structure = Structure::new(
            vec![
                atom(0, "O5", "O", [0.0, 0.0, 0.0]),
                atom(0, "C1", "C", [1.43, 0.0, 0.0]),
                atom(0, "C5", "C", [-0.5, 1.34, 0.0]),
                atom(1, "S", "S", [0.0, -0.6, 1.48]),
            ],
            2,
            &[],
        );
        assert_eq!(structure.neighbors[0], vec![1, 2]);
        assert!(!structure.neighbors[3].contains(&0));
    }

    #[test]
    fn the_signed_volume_is_negative_for_a_heteroatom_on_the_right() {
        // Fischer projection at the origin: the chain runs up (top) and down
        // (bottom) away from the viewer, the right-hand group points right
        // towards the viewer.
        let structure = Structure::new(
            vec![
                atom(0, "C", "C", [0.0, 0.0, 0.0]),
                atom(0, "Ctop", "C", [0.0, 1.0, -0.7]),
                atom(0, "Cbottom", "C", [0.0, -1.0, -0.7]),
                atom(0, "Oright", "O", [1.0, 0.0, 0.7]),
            ],
            1,
            &[],
        );
        let volume = structure.normalised_volume(0, 1, 2, 3);
        assert!(volume < -DISTORTED_VOLUME);
        assert_eq!(fischer_code(volume), '2');
    }
}
