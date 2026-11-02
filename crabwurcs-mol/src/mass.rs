//! Formula, mass, composition and fragment ions of a glycan.
//!
//! The formula of each residue comes from the molecular constructor, so it is
//! the chemistry SMILES and MOL output is written from. Mass does not need one
//! concrete molecule, though: linkage positions and stereochemistry do not
//! change it, so compositions, `?` linkages and undefined fragments have a
//! mass even though they have no molecular graph. Only what really leaves the
//! mass open is an error: repeat ranges, probabilities and residues without
//! chemistry.

use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::fmt;
use std::str::FromStr;

use chematic::core::{Atom, BondOrder, Element, Molecule, MoleculeBuilder, bond_order_sum};
use crabwurcs_core::{
    Linkage, Monosaccharide, RepeatCount, ResidueGraph, ResidueKind,
    classify_residue_configuration, residue_from_kind,
};
use petgraph::Direction;
use petgraph::graph::NodeIndex;
use petgraph::visit::EdgeRef;
use thiserror::Error;

use crate::MolError;
use crate::construct::{self, add_modification};
use crate::map::MapGraph;

/// Mass of the electron in u (CODATA 2018).
const ELECTRON: f64 = 0.000_548_579_909;

/// Monoisotopic mass (AME2020) and conventional atomic weight (IUPAC 2021).
fn element_masses(symbol: &str) -> Option<(f64, f64)> {
    Some(match symbol {
        "H" => (1.007_825_032_23, 1.008),
        "Li" => (7.016_003_436_6, 6.94),
        "B" => (11.009_305_36, 10.81),
        "C" => (12.0, 12.011),
        "N" => (14.003_074_004_43, 14.007),
        "O" => (15.994_914_619_57, 15.999),
        "F" => (18.998_403_162_73, 18.998_403_162),
        "Na" => (22.989_769_282_0, 22.989_769_28),
        "Si" => (27.976_926_534_65, 28.085),
        "P" => (30.973_761_998_42, 30.973_761_998),
        "S" => (31.972_071_174_4, 32.06),
        "Cl" => (34.968_852_682, 35.45),
        "K" => (38.963_706_486_4, 39.098_3),
        "Se" => (79.916_521_8, 78.971),
        "Br" => (78.918_337_6, 79.904),
        "I" => (126.904_471_9, 126.904_47),
        _ => return None,
    })
}

#[derive(Debug, Error)]
pub enum MassError {
    #[error("the glycan has no residues")]
    Empty,
    #[error("residue {0} has no defined chemistry, so it has no mass")]
    UndefinedResidue(String),
    #[error("a probability makes this an ensemble of species, not one mass")]
    Ensemble,
    #[error("the repeat count is not a single number")]
    VariableRepeat,
    #[error("a glycan with more than one repeat unit is not supported")]
    MultipleRepeats,
    #[error("an exact repeat count must be positive")]
    EmptyRepeat,
    #[error("the residues are not joined into one glycan")]
    Disconnected,
    #[error("the glycan has no free reducing end to derivatize")]
    NoReducingEnd,
    #[error("permethylation is defined for free or reduced glycans without sulfate or phosphate")]
    UnsupportedPermethylation,
    #[error("no atomic mass is known for element {0}")]
    UnknownElement(String),
    #[error("fragment ions need one fully defined tree, but the glycan has {0}")]
    NotATree(&'static str),
    #[error("unknown {kind}: {value}")]
    UnknownName { kind: &'static str, value: String },
    #[error(transparent)]
    Molecule(#[from] MolError),
}

pub type MassResult<T> = Result<T, MassError>;

/// An elemental formula. Counts may be negative while it is a difference.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct Formula(BTreeMap<String, i64>);

impl Formula {
    /// A formula from `(element, count)` pairs.
    pub fn of(elements: &[(&str, i64)]) -> Self {
        let mut formula = Self::default();
        for (element, count) in elements {
            formula.add_element(element, *count);
        }
        formula
    }

    /// The number of atoms of `element`.
    pub fn count(&self, element: &str) -> i64 {
        self.0.get(element).copied().unwrap_or(0)
    }

    pub fn elements(&self) -> impl Iterator<Item = (&str, i64)> {
        self.0
            .iter()
            .map(|(element, count)| (element.as_str(), *count))
    }

    fn add_element(&mut self, element: &str, count: i64) {
        let total = self.0.entry(element.to_owned()).or_insert(0);
        *total += count;
        if *total == 0 {
            self.0.remove(element);
        }
    }

    fn add(&mut self, other: &Self, times: i64) {
        for (element, count) in &other.0 {
            self.add_element(element, count * times);
        }
    }

    fn mass(&self, pick: fn((f64, f64)) -> f64) -> MassResult<f64> {
        self.0.iter().try_fold(0.0, |total, (element, count)| {
            let masses = element_masses(element)
                .ok_or_else(|| MassError::UnknownElement(element.clone()))?;
            Ok(total + pick(masses) * *count as f64)
        })
    }

    pub fn monoisotopic_mass(&self) -> MassResult<f64> {
        self.mass(|(monoisotopic, _)| monoisotopic)
    }

    pub fn average_mass(&self) -> MassResult<f64> {
        self.mass(|(_, average)| average)
    }
}

/// Hill order: carbon, hydrogen, then the rest alphabetically.
impl fmt::Display for Formula {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        let write = |formatter: &mut fmt::Formatter<'_>, element: &str, count: i64| {
            if count == 1 {
                write!(formatter, "{element}")
            } else {
                write!(formatter, "{element}{count}")
            }
        };
        for element in ["C", "H"] {
            if self.count(element) != 0 {
                write(formatter, element, self.count(element))?;
            }
        }
        for (element, count) in self.elements() {
            if !matches!(element, "C" | "H") {
                write(formatter, element, count)?;
            }
        }
        Ok(())
    }
}

/// What the reducing end of a glycan carries.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum ReducingEnd {
    /// The free hemiacetal.
    #[default]
    Free,
    /// The alditol from borohydride reduction (+H2).
    Reduced,
    /// 2-aminobenzamide by reductive amination (+C7H8N2).
    TwoAb,
    /// 2-aminobenzoic acid by reductive amination (+C7H7NO).
    TwoAa,
    /// Procainamide by reductive amination (+C13H21N3).
    Procainamide,
}

impl ReducingEnd {
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Free => "free",
            Self::Reduced => "reduced",
            Self::TwoAb => "2-ab",
            Self::TwoAa => "2-aa",
            Self::Procainamide => "procainamide",
        }
    }

    /// What the derivatization adds to the free glycan.
    fn increment(self) -> Formula {
        match self {
            Self::Free => Formula::default(),
            Self::Reduced => Formula::of(&[("H", 2)]),
            Self::TwoAb => Formula::of(&[("C", 7), ("H", 8), ("N", 2)]),
            Self::TwoAa => Formula::of(&[("C", 7), ("H", 7), ("N", 1), ("O", 1)]),
            Self::Procainamide => Formula::of(&[("C", 13), ("H", 21), ("N", 3)]),
        }
    }
}

impl FromStr for ReducingEnd {
    type Err = MassError;

    fn from_str(value: &str) -> Result<Self, Self::Err> {
        match value.trim().to_ascii_lowercase().as_str() {
            "free" | "oh" | "-oh" => Ok(Self::Free),
            "reduced" | "alditol" => Ok(Self::Reduced),
            "2-ab" | "2ab" => Ok(Self::TwoAb),
            "2-aa" | "2aa" => Ok(Self::TwoAa),
            "procainamide" | "proa" => Ok(Self::Procainamide),
            other => Err(MassError::UnknownName {
                kind: "reducing end",
                value: other.to_owned(),
            }),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct MassOptions {
    pub reducing_end: ReducingEnd,
    /// Every free OH and NH (and carboxyl) methylated.
    pub permethylated: bool,
}

/// The ion a neutral glycan is observed as.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Adduct {
    Proton,
    Sodium,
    Potassium,
    Ammonium,
    Lithium,
    Deprotonated,
    Chloride,
    Formate,
    Acetate,
}

impl Adduct {
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Proton => "H",
            Self::Sodium => "Na",
            Self::Potassium => "K",
            Self::Ammonium => "NH4",
            Self::Lithium => "Li",
            Self::Deprotonated => "-H",
            Self::Chloride => "Cl",
            Self::Formate => "HCOO",
            Self::Acetate => "CH3COO",
        }
    }

    /// What one adduct adds, and its charge.
    fn change(self) -> (Formula, f64) {
        match self {
            Self::Proton => (Formula::of(&[("H", 1)]), 1.0),
            Self::Sodium => (Formula::of(&[("Na", 1)]), 1.0),
            Self::Potassium => (Formula::of(&[("K", 1)]), 1.0),
            Self::Ammonium => (Formula::of(&[("N", 1), ("H", 4)]), 1.0),
            Self::Lithium => (Formula::of(&[("Li", 1)]), 1.0),
            Self::Deprotonated => (Formula::of(&[("H", -1)]), -1.0),
            Self::Chloride => (Formula::of(&[("Cl", 1)]), -1.0),
            Self::Formate => (Formula::of(&[("C", 1), ("H", 1), ("O", 2)]), -1.0),
            Self::Acetate => (Formula::of(&[("C", 2), ("H", 3), ("O", 2)]), -1.0),
        }
    }

    /// The mass of one adduct: an ion carries one electron less per positive
    /// charge and one more per negative charge.
    fn mass(self, mass: fn(&Formula) -> MassResult<f64>) -> MassResult<f64> {
        let (formula, charge) = self.change();
        Ok(mass(&formula)? - charge * ELECTRON)
    }

    fn mz(
        self,
        neutral: f64,
        count: u32,
        mass: fn(&Formula) -> MassResult<f64>,
    ) -> MassResult<f64> {
        let count = f64::from(count.max(1));
        Ok((neutral + count * self.mass(mass)?) / count)
    }
}

impl FromStr for Adduct {
    type Err = MassError;

    fn from_str(value: &str) -> Result<Self, Self::Err> {
        let name = value.trim().to_ascii_lowercase();
        // `[M+Na]+`, `M+Na` and `Na+` all name sodium.
        let name = name.trim_start_matches('[').trim_start_matches('m');
        let name = name.trim_end_matches(['+', ']']);
        let name = name.trim_start_matches('+');
        match name {
            "h" | "proton" => Ok(Self::Proton),
            "na" | "sodium" => Ok(Self::Sodium),
            "k" | "potassium" => Ok(Self::Potassium),
            "nh4" | "ammonium" => Ok(Self::Ammonium),
            "li" | "lithium" => Ok(Self::Lithium),
            "-h" | "-h-" | "h-" | "-h]-" | "deprotonated" => Ok(Self::Deprotonated),
            "cl" | "cl-" | "cl]-" | "chloride" => Ok(Self::Chloride),
            "hcoo" | "hcoo-" | "hcoo]-" | "formate" => Ok(Self::Formate),
            "ch3coo" | "ch3coo-" | "ch3coo]-" | "acetate" => Ok(Self::Acetate),
            _ => Err(MassError::UnknownName {
                kind: "adduct",
                value: value.trim().to_owned(),
            }),
        }
    }
}

/// The neutral formula and masses of a glycan.
#[derive(Debug, Clone, PartialEq)]
pub struct GlycanMass {
    pub formula: Formula,
    pub monoisotopic: f64,
    pub average: f64,
}

impl GlycanMass {
    fn of(formula: Formula) -> MassResult<Self> {
        Ok(Self {
            monoisotopic: formula.monoisotopic_mass()?,
            average: formula.average_mass()?,
            formula,
        })
    }

    /// Monoisotopic m/z with `count` adducts, each carrying one charge.
    pub fn mz(&self, adduct: Adduct, count: u32) -> MassResult<f64> {
        adduct.mz(self.monoisotopic, count, Formula::monoisotopic_mass)
    }

    /// Average m/z with `count` adducts, each carrying one charge.
    pub fn average_mz(&self, adduct: Adduct, count: u32) -> MassResult<f64> {
        adduct.mz(self.average, count, Formula::average_mass)
    }
}

/// A formula with the OH and NH hydrogens permethylation would replace.
#[derive(Debug, Clone, Default)]
struct Tally {
    formula: Formula,
    methyl_sites: i64,
}

impl Tally {
    fn of_molecule(molecule: &Molecule) -> Self {
        let mut tally = Self::default();
        let mut quaternary = 0;
        for (index, atom) in molecule.atoms() {
            let mut hydrogens = i64::from(molecule.implicit_hydrogen_count(index));
            if atom.element == Element::N && bond_order_sum(molecule, index) >= 4 {
                hydrogens = 0;
                quaternary += 1;
            }
            tally.formula.add_element(atom.element.symbol(), 1);
            tally.formula.add_element("H", hydrogens);
            if matches!(atom.element, Element::O | Element::N) {
                tally.methyl_sites += hydrogens;
            }
        }
        // Phosphocholine is neutral as an inner salt: the ammonium's charge
        // is balanced by the phosphate, which has lost its proton.
        if tally.formula.count("P") != 0 {
            tally.formula.add_element("H", -quaternary);
            tally.methyl_sites -= quaternary;
        }
        tally
    }

    fn add(&mut self, other: &Self, times: i64) {
        self.formula.add(&other.formula, times);
        self.methyl_sites += other.methyl_sites * times;
    }

    /// Two hydroxyls condensing into a glycosidic bond.
    fn condensation() -> Self {
        Self {
            formula: Formula::of(&[("H", -2), ("O", -1)]),
            methyl_sites: -2,
        }
    }

    /// A hydroxyl, which a substituent or a bridge replaces.
    fn methanol() -> Self {
        Self {
            formula: Formula::of(&[("C", 1), ("H", 4), ("O", 1)]),
            methyl_sites: 1,
        }
    }
}

/// What a residue-level substituent changes on a backbone hydroxyl.
fn substituent_tally(descriptor: &str) -> MassResult<Tally> {
    let mut builder = MoleculeBuilder::new();
    let carbon = builder.add_atom(Atom::new(Element::C));
    add_modification(&mut builder, carbon, descriptor)?;
    let mut tally = Tally::of_molecule(&builder.build());
    tally.add(&Tally::methanol(), -1);
    Ok(tally)
}

/// What a MAP bridging two residues changes on their two hydroxyls.
fn bridge_tally(linkage: &Linkage, map_code: &str) -> MassResult<Tally> {
    let map = MapGraph::parse(map_code)?;
    let mut builder = MoleculeBuilder::new();
    let emitted = map.build(&mut builder)?;
    let parent_star = linkage.parent_modification_position.unwrap_or(1);
    let child_star =
        linkage
            .child_modification_position
            .unwrap_or(if emitted.attachments.contains_key(&2) {
                2
            } else {
                1
            });
    for star in [parent_star, child_star] {
        let (atom, order) =
            emitted
                .attachments
                .get(&star)
                .copied()
                .ok_or_else(|| MolError::InvalidMap {
                    offset: 0,
                    message: format!("MAP has no attachment star {star}"),
                })?;
        let carbon = builder.add_atom(Atom::new(Element::C));
        builder
            .add_bond(carbon, atom, order)
            .map_err(|cause| MolError::InvalidValence(cause.to_string()))?;
    }
    let mut tally = Tally::of_molecule(&builder.build());
    tally.add(&Tally::methanol(), -2);
    Ok(tally)
}

fn bond_tally(linkage: &Linkage) -> MassResult<Tally> {
    if linkage.parent_probability.is_some() || linkage.child_probability.is_some() {
        return Err(MassError::Ensemble);
    }
    match linkage.map_code.as_deref() {
        Some(map_code) => bridge_tally(linkage, map_code),
        None => Ok(Tally::condensation()),
    }
}

/// A backbone without an anomeric carbon (an alditol, glycerol), counted
/// carbon by carbon: the molecular constructor only builds sugars.
fn alditol_tally(residue: &Monosaccharide) -> MassResult<Tally> {
    let descriptors = residue.skeleton_code.chars().collect::<Vec<_>>();
    let mut tally = Tally::default();
    for (offset, descriptor) in descriptors.iter().copied().enumerate() {
        let terminal = offset == 0 || offset + 1 == descriptors.len();
        let oxygens = construct::default_oxygen_bonds(descriptor, terminal)?;
        let hydroxyls = oxygens
            .iter()
            .filter(|order| **order == BondOrder::Single)
            .count() as i64;
        let oxygen_bonds = oxygens
            .iter()
            .map(|order| i64::from(order.order_int()))
            .sum::<i64>();
        let backbone_bonds = if terminal { 1 } else { 2 };
        tally.formula.add_element("C", 1);
        tally.formula.add_element("O", oxygens.len() as i64);
        tally
            .formula
            .add_element("H", 4 - backbone_bonds - oxygen_bonds + hydroxyls);
        tally.methyl_sites += hydroxyls;
    }
    for modification in &residue.modifications {
        tally.add(&substituent_tally(&modification.descriptor)?, 1);
    }
    Ok(tally)
}

fn residue_tally(residue: &Monosaccharide) -> MassResult<Tally> {
    if let Some(name) = &residue.display_name {
        return Err(MassError::UndefinedResidue(name.clone()));
    }
    if let Some(kind) = residue.residue_kind
        && kind.unique_residue().is_none()
    {
        return Err(MassError::UndefinedResidue(
            kind.canonical_name().to_owned(),
        ));
    }
    if residue
        .modifications
        .iter()
        .any(|m| m.probability.is_some())
    {
        return Err(MassError::Ensemble);
    }
    if residue.anomeric_prefix.is_empty() {
        return alditol_tally(residue);
    }
    // A substituent at an unknown position has a mass all the same; only the
    // molecular constructor needs to know where it goes.
    let (placed, unplaced): (Vec<_>, Vec<_>) = residue
        .modifications
        .iter()
        .cloned()
        .partition(|modification| modification.position.0 != 0);
    let mut free = residue.clone();
    free.modifications = placed;
    let mut graph = ResidueGraph::new();
    graph.add_residue(free);
    let mut tally = Tally::of_molecule(&construct::construct_molecule(&graph)?);
    for modification in unplaced {
        tally.add(&substituent_tally(&modification.descriptor)?, 1);
    }
    Ok(tally)
}

/// Residue and bond tallies by node and edge index, and the repeat unit.
struct Parts {
    residues: Vec<Tally>,
    bonds: HashMap<usize, Tally>,
    /// The one exact repeat: its unit, the bond between units, and the count.
    repeat: Option<(BTreeSet<usize>, Tally, i64)>,
}

impl Parts {
    fn of(graph: &ResidueGraph) -> MassResult<Self> {
        let inner = graph.inner();
        let residues = inner
            .node_indices()
            .map(|node| residue_tally(&inner[node]))
            .collect::<MassResult<Vec<_>>>()?;
        let mut bonds = HashMap::new();
        let mut repeat = None;
        for edge in inner.edge_references() {
            let linkage = edge.weight();
            let tally = bond_tally(linkage)?;
            match &linkage.repeat {
                None => {
                    bonds.insert(edge.id().index(), tally);
                }
                Some(RepeatCount::Unknown | RepeatCount::Range { .. }) => {
                    return Err(MassError::VariableRepeat);
                }
                Some(RepeatCount::Exact(0)) => return Err(MassError::EmptyRepeat),
                Some(RepeatCount::Exact(count)) => {
                    if repeat.is_some() {
                        return Err(MassError::MultipleRepeats);
                    }
                    repeat = Some((graph.repeat_unit(edge.id()), tally, i64::from(*count)));
                }
            }
        }
        Ok(Self {
            residues,
            bonds,
            repeat,
        })
    }

    /// The residues of `nodes` and the bonds between them.
    fn subset(&self, graph: &ResidueGraph, nodes: &BTreeSet<usize>) -> Tally {
        let mut tally = Tally::default();
        for node in nodes {
            tally.add(&self.residues[*node], 1);
        }
        for edge in graph.inner().edge_references() {
            if nodes.contains(&edge.source().index())
                && nodes.contains(&edge.target().index())
                && let Some(bond) = self.bonds.get(&edge.id().index())
            {
                tally.add(bond, 1);
            }
        }
        tally
    }
}

/// Whether the reducing-end residue still has its anomeric carbon free.
fn has_free_reducing_end(graph: &ResidueGraph) -> MassResult<bool> {
    let Some(root) = graph.root() else {
        return Ok(false);
    };
    if graph
        .inner()
        .edge_references()
        .any(|edge| edge.weight().cyclic)
    {
        return Ok(false);
    }
    let residue = &graph.inner()[root];
    if residue.anomeric_prefix.is_empty() {
        return Ok(false);
    }
    let descriptors = construct::carbon_descriptors(residue)?;
    if !descriptors
        .iter()
        .any(|descriptor| matches!(descriptor, 'a' | 'u' | 'U' | 'o' | 'O' | 'k'))
    {
        return Ok(false);
    }
    let anomeric = construct::effective_anomeric_position(residue);
    Ok(!graph
        .inner()
        .edges_directed(root, Direction::Outgoing)
        .any(|edge| usize::from(edge.weight().parent_position.0) == anomeric))
}

/// Apply the reducing-end derivatization to a fragment that holds the root.
fn derivatize_reducing_end(
    tally: &mut Tally,
    graph: &ResidueGraph,
    options: &MassOptions,
) -> MassResult<()> {
    if options.reducing_end == ReducingEnd::Free {
        return Ok(());
    }
    if !has_free_reducing_end(graph)? {
        return Err(MassError::NoReducingEnd);
    }
    if options.permethylated && options.reducing_end != ReducingEnd::Reduced {
        return Err(MassError::UnsupportedPermethylation);
    }
    tally.formula.add(&options.reducing_end.increment(), 1);
    if options.reducing_end == ReducingEnd::Reduced {
        // The ring opens: the anomeric hydroxyl stays and the ring oxygen
        // becomes a hydroxyl of its own.
        tally.methyl_sites += 1;
    }
    Ok(())
}

fn permethylate(tally: &Tally, options: &MassOptions) -> MassResult<Formula> {
    let mut formula = tally.formula.clone();
    if options.permethylated {
        if formula.count("S") != 0 || formula.count("P") != 0 {
            return Err(MassError::UnsupportedPermethylation);
        }
        let sites = tally.methyl_sites;
        formula.add(&Formula::of(&[("C", 1), ("H", 2)]), sites);
    }
    Ok(formula)
}

/// The elemental formula of a neutral glycan.
pub fn glycan_formula(graph: &ResidueGraph, options: &MassOptions) -> MassResult<Formula> {
    if graph.node_count() == 0 {
        return Err(MassError::Empty);
    }
    let parts = Parts::of(graph)?;
    let all = graph
        .inner()
        .node_indices()
        .map(|node| node.index())
        .collect::<BTreeSet<_>>();
    let mut tally = parts.subset(graph, &all);

    if graph.is_composition() {
        // A composition is one glycan whose bonds are not stated.
        tally.add(&Tally::condensation(), graph.node_count() as i64 - 1);
    } else {
        let components = petgraph::algo::connected_components(graph.inner());
        if components > 1 + graph.undefined_linkages().len() {
            return Err(MassError::Disconnected);
        }
    }
    for _ in graph.undefined_linkages() {
        tally.add(&Tally::condensation(), 1);
    }
    for modification in graph.undefined_modifications() {
        let descriptor = modification.map_code.trim_start_matches('*');
        tally.add(&substituent_tally(descriptor)?, 1);
    }
    if let Some((unit, junction, count)) = &parts.repeat {
        tally.add(&parts.subset(graph, unit), count - 1);
        tally.add(junction, count - 1);
    }
    derivatize_reducing_end(&mut tally, graph, options)?;
    permethylate(&tally, options)
}

/// The formula and masses of a neutral glycan.
pub fn glycan_mass(graph: &ResidueGraph, options: &MassOptions) -> MassResult<GlycanMass> {
    GlycanMass::of(glycan_formula(graph, options)?)
}

/// Residue classes and substituents with their counts, in the order glycomics
/// writes them (`HexNAc4Hex5Fuc1NeuAc2`).
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct Composition(Vec<(String, u32)>);

impl Composition {
    pub fn classes(&self) -> &[(String, u32)] {
        &self.0
    }

    /// The count of one class.
    pub fn count(&self, class: &str) -> u32 {
        self.0
            .iter()
            .find(|(name, _)| name == class)
            .map_or(0, |(_, count)| *count)
    }
}

impl fmt::Display for Composition {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        for (class, count) in &self.0 {
            write!(formatter, "{class}{count}")?;
        }
        Ok(())
    }
}

/// The order of classes and substituents in a composition. Anything else
/// sorts alphabetically between the two groups.
const COMPOSITION_ORDER: &[&str] = &[
    "HexNAc", "Hex", "Fuc", "dHex", "NeuAc", "NeuGc", "Neu", "Kdn", "HexA", "HexN", "Pen",
];
const SUBSTITUENT_ORDER: &[&str] = &["S", "P", "Ac", "Me", "Sub"];

fn composition_class(residue: &Monosaccharide) -> String {
    let Some((kind, _)) = classify_residue_configuration(residue) else {
        return "Unknown".to_owned();
    };
    for (family, name) in [
        (ResidueKind::HexNAc, "HexNAc"),
        (ResidueKind::HexN, "HexN"),
        (ResidueKind::HexA, "HexA"),
        (ResidueKind::Hex, "Hex"),
        (ResidueKind::DHexNAc, "dHexNAc"),
        (ResidueKind::DDHex, "ddHex"),
        (ResidueKind::Pen, "Pen"),
    ] {
        if family.matches_family(kind) {
            return name.to_owned();
        }
    }
    match kind {
        ResidueKind::Fuc => "Fuc".to_owned(),
        _ if ResidueKind::DHex.matches_family(kind) => "dHex".to_owned(),
        ResidueKind::Neu5Ac => "NeuAc".to_owned(),
        ResidueKind::Neu5Gc => "NeuGc".to_owned(),
        ResidueKind::Assigned => residue
            .display_name
            .clone()
            .unwrap_or_else(|| "Unknown".to_owned()),
        other => other.canonical_name().to_owned(),
    }
}

fn substituent_class(descriptor: &str) -> &'static str {
    match descriptor {
        "OSO/3=O/3=O" | "NSO/3=O/3=O" => "S",
        "OPO/3O/3=O" => "P",
        "OCC/3=O" | "NCC/3=O" => "Ac",
        "OC" => "Me",
        _ => "Sub",
    }
}

/// The residue classes of a glycan and the substituents its residues carry
/// beyond their names (the sulfates of GlcNS6S, not the acetyl of GlcNAc).
pub fn glycan_composition(graph: &ResidueGraph) -> Composition {
    let mut counts = BTreeMap::<String, u32>::new();
    for residue in graph.inner().node_weights() {
        *counts.entry(composition_class(residue)).or_default() += 1;
        let named = classify_residue_configuration(residue)
            .and_then(|(kind, _)| residue_from_kind(kind).ok())
            .map(|template| template.modifications)
            .unwrap_or_default();
        for modification in &residue.modifications {
            let in_the_name = named.iter().any(|known| {
                known.position == modification.position
                    && known.descriptor == modification.descriptor
            });
            if !in_the_name {
                let class = substituent_class(&modification.descriptor);
                *counts.entry(class.to_owned()).or_default() += 1;
            }
        }
    }
    for modification in graph.undefined_modifications() {
        let class = substituent_class(modification.map_code.trim_start_matches('*'));
        *counts.entry(class.to_owned()).or_default() += 1;
    }

    let rank = |class: &str| {
        if let Some(index) = COMPOSITION_ORDER.iter().position(|known| *known == class) {
            (0, index)
        } else if let Some(index) = SUBSTITUENT_ORDER.iter().position(|known| *known == class) {
            (2, index)
        } else {
            (1, 0)
        }
    };
    let mut classes = counts.into_iter().collect::<Vec<_>>();
    classes.sort_by(|(left, _), (right, _)| {
        rank(left).cmp(&rank(right)).then_with(|| left.cmp(right))
    });
    Composition(classes)
}

/// Domon–Costello glycosidic cleavage ions.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum IonKind {
    /// Non-reducing fragment without the glycosidic oxygen.
    B,
    /// Non-reducing fragment with the glycosidic oxygen.
    C,
    /// Reducing fragment with the glycosidic oxygen.
    Y,
    /// Reducing fragment without the glycosidic oxygen.
    Z,
}

impl IonKind {
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::B => "B",
            Self::C => "C",
            Self::Y => "Y",
            Self::Z => "Z",
        }
    }

    /// Whether the ion keeps the reducing end.
    pub const fn keeps_reducing_end(self) -> bool {
        matches!(self, Self::Y | Self::Z)
    }
}

/// One fragment ion from cleaving one glycosidic bond.
#[derive(Debug, Clone, PartialEq)]
pub struct Fragment {
    pub kind: IonKind,
    /// `B2`, `Y3α`: residues counted along the chain from the end the ion
    /// keeps, with one Greek letter per branch point on the way to the bond
    /// (α for the heaviest branch).
    pub label: String,
    /// Edge index of the cleaved bond.
    pub edge: usize,
    /// Node indices of the residues the ion keeps.
    pub residues: BTreeSet<usize>,
    /// Neutral formula of the fragment.
    pub formula: Formula,
    pub monoisotopic: f64,
    /// Monoisotopic m/z with the adduct the fragments were asked for.
    pub mz: f64,
}

const BRANCH_LETTERS: [char; 8] = [
    '\u{03B1}', '\u{03B2}', '\u{03B3}', '\u{03B4}', '\u{03B5}', '\u{03B6}', '\u{03B7}', '\u{03B8}',
];

fn children(graph: &ResidueGraph, node: NodeIndex) -> Vec<NodeIndex> {
    graph
        .inner()
        .edges_directed(node, Direction::Outgoing)
        .map(|edge| edge.target())
        .collect()
}

/// Residues on the longest path from `node` down its subtree, `node` included.
fn chain_length(graph: &ResidueGraph, node: NodeIndex) -> usize {
    1 + children(graph, node)
        .into_iter()
        .map(|child| chain_length(graph, child))
        .max()
        .unwrap_or(0)
}

/// The depth and branch letters of every residue below `root`.
fn branch_paths(
    graph: &ResidueGraph,
    parts: &Parts,
    root: NodeIndex,
) -> MassResult<HashMap<usize, (usize, String)>> {
    let mut paths = HashMap::from([(root.index(), (0usize, String::new()))]);
    let mut stack = vec![root];
    while let Some(node) = stack.pop() {
        let (depth, letters) = paths[&node.index()].clone();
        let mut next = children(graph, node);
        if next.len() > 1 {
            // Heaviest branch first; the node index settles a tie.
            let mut weighed = next
                .iter()
                .map(|child| {
                    let subtree = parts.subset(graph, &graph.subtree(*child));
                    Ok((subtree.formula.monoisotopic_mass()?, *child))
                })
                .collect::<MassResult<Vec<_>>>()?;
            weighed.sort_by(|left, right| right.0.total_cmp(&left.0).then(left.1.cmp(&right.1)));
            next = weighed.into_iter().map(|(_, child)| child).collect();
        }
        let branching = next.len() > 1;
        for (rank, child) in next.iter().enumerate() {
            let mut path = letters.clone();
            if branching {
                path.push(BRANCH_LETTERS[rank.min(BRANCH_LETTERS.len() - 1)]);
            }
            paths.insert(child.index(), (depth + 1, path));
            stack.push(*child);
        }
    }
    Ok(paths)
}

/// Every B, C, Y and Z ion from one glycosidic cleavage, with `adduct`.
pub fn glycan_fragments(
    graph: &ResidueGraph,
    options: &MassOptions,
    adduct: Adduct,
) -> MassResult<Vec<Fragment>> {
    let Some(root) = graph.root() else {
        return Err(MassError::Empty);
    };
    if graph.is_composition() {
        return Err(MassError::NotATree("no linkages"));
    }
    if !graph.undefined_linkages().is_empty() || !graph.undefined_modifications().is_empty() {
        return Err(MassError::NotATree("undefined fragments"));
    }
    if graph
        .inner()
        .edge_weights()
        .any(|linkage| linkage.repeat.is_some() || linkage.cyclic)
    {
        return Err(MassError::NotATree("a repeat or a cycle"));
    }
    if graph.subtree(root).len() != graph.node_count() {
        return Err(MassError::Disconnected);
    }

    let parts = Parts::of(graph)?;
    let paths = branch_paths(graph, &parts, root)?;
    let all = graph.subtree(root);
    let water = Formula::of(&[("H", 2), ("O", 1)]);
    let methylene = Formula::of(&[("C", 1), ("H", 2)]);
    let methanol = Tally::methanol().formula;

    let mut fragments = Vec::new();
    for edge in graph.inner().edge_references() {
        if edge.weight().map_code.is_some() {
            continue;
        }
        let child = edge.target();
        let kept_non_reducing = graph.subtree(child);
        let kept_reducing = all
            .difference(&kept_non_reducing)
            .copied()
            .collect::<BTreeSet<_>>();

        let non_reducing = permethylate(&parts.subset(graph, &kept_non_reducing), options)?;
        let mut reducing_tally = parts.subset(graph, &kept_reducing);
        derivatize_reducing_end(&mut reducing_tally, graph, options)?;
        let reducing = permethylate(&reducing_tally, options)?;

        // A cleaved permethylated bond leaves a free hydroxyl where the bond
        // was, and a B ion loses methanol from the anomeric methoxy group.
        let mut c = non_reducing.clone();
        let mut y = reducing;
        let mut b = non_reducing;
        if options.permethylated {
            c.add(&methylene, -1);
            y.add(&methylene, -1);
            b.add(&methanol, -1);
        } else {
            b.add(&water, -1);
        }
        let mut z = y.clone();
        z.add(&water, -1);

        let (depth, letters) = &paths[&child.index()];
        let chain = chain_length(graph, child);
        for (kind, formula, number, residues) in [
            (IonKind::B, b, chain, &kept_non_reducing),
            (IonKind::C, c, chain, &kept_non_reducing),
            (IonKind::Y, y, *depth, &kept_reducing),
            (IonKind::Z, z, *depth, &kept_reducing),
        ] {
            let monoisotopic = formula.monoisotopic_mass()?;
            // The formula of a B ion is the neutral that takes the one adduct
            // carrying its charge, so B + H+ is the oxonium ion.
            let mz = adduct.mz(monoisotopic, 1, Formula::monoisotopic_mass)?;
            fragments.push(Fragment {
                kind,
                label: format!("{}{number}{letters}", kind.as_str()),
                edge: edge.id().index(),
                residues: residues.clone(),
                formula,
                monoisotopic,
                mz,
            });
        }
    }
    fragments.sort_by(|left, right| {
        left.kind
            .cmp(&right.kind)
            .then_with(|| left.residues.len().cmp(&right.residues.len()))
            .then_with(|| left.label.cmp(&right.label))
    });
    Ok(fragments)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crabwurcs_core::parse_wurcs;

    const GLCNAC: &str = "WURCS=2.0/1,1,0/[a2122h-1b_1-5_2*NCC/3=O]/1/";
    /// Man5GlcNAc2 with a free reducing end.
    const MAN5: &str = "WURCS=2.0/4,7,6/[u2122h_2*NCC/3=O][a2122h-1b_1-5_2*NCC/3=O][a1122h-1b_1-5][a1122h-1a_1-5]/1-2-3-4-4-4-4/a4-b1_b4-c1_c3-d1_c6-e1_e3-f1_e6-g1";

    fn close(found: f64, expected: f64, tolerance: f64) {
        assert!(
            (found - expected).abs() <= tolerance,
            "{found} is not within {tolerance} of {expected}"
        );
    }

    fn free() -> MassOptions {
        MassOptions::default()
    }

    #[test]
    fn formulas_print_in_hill_order() {
        assert_eq!(
            Formula::of(&[("O", 6), ("N", 1), ("H", 15), ("C", 8), ("S", 1)]).to_string(),
            "C8H15NO6S"
        );
    }

    #[test]
    fn monosaccharides_have_their_textbook_formulas() {
        for (wurcs, formula, mass) in [
            ("WURCS=2.0/1,1,0/[a2122h-1b_1-5]/1/", "C6H12O6", 180.06339),
            (GLCNAC, "C8H15NO6", 221.08994),
            ("WURCS=2.0/1,1,0/[a1221m-1a_1-5]/1/", "C6H12O5", 164.06847),
            ("WURCS=2.0/1,1,0/[a212h-1b_1-5]/1/", "C5H10O5", 150.05282),
            ("WURCS=2.0/1,1,0/[a2122A-1b_1-5]/1/", "C6H10O7", 194.04265),
            (
                "WURCS=2.0/1,1,0/[Aad21122h-2a_2-6_5*NCC/3=O]/1/",
                "C11H19NO9",
                309.10599,
            ),
            (
                "WURCS=2.0/1,1,0/[Aad21122h-2a_2-6_5*NCCO/3=O]/1/",
                "C11H19NO10",
                325.10090,
            ),
            ("WURCS=2.0/1,1,0/[AUd21122h]/1/", "C9H16O9", 268.07943),
        ] {
            let mass_of = glycan_mass(&parse_wurcs(wurcs).unwrap(), &free()).unwrap();
            assert_eq!(mass_of.formula.to_string(), formula, "{wurcs}");
            close(mass_of.monoisotopic, mass, 0.0001);
        }
    }

    #[test]
    fn man5_reproduces_the_published_masses() {
        let graph = parse_wurcs(MAN5).unwrap();
        let native = glycan_mass(&graph, &free()).unwrap();
        assert_eq!(native.formula.to_string(), "C46H78N2O36");
        close(native.monoisotopic, 1234.4333, 0.0002);
        close(native.average, 1235.08, 0.05);
        close(native.mz(Adduct::Sodium, 1).unwrap(), 1257.4226, 0.0003);
        close(native.mz(Adduct::Proton, 1).unwrap(), 1235.4406, 0.0003);
        close(
            native.mz(Adduct::Deprotonated, 1).unwrap(),
            1233.4261,
            0.0003,
        );
        close(native.mz(Adduct::Proton, 2).unwrap(), 618.2239, 0.0003);

        // Permethylated Man5 is the 1579.8 peak of every N-glycan MALDI.
        let methylated = glycan_mass(
            &graph,
            &MassOptions {
                permethylated: true,
                ..free()
            },
        )
        .unwrap();
        close(methylated.mz(Adduct::Sodium, 1).unwrap(), 1579.7835, 0.001);

        let reduced = glycan_mass(
            &graph,
            &MassOptions {
                reducing_end: ReducingEnd::Reduced,
                ..free()
            },
        )
        .unwrap();
        close(reduced.monoisotopic - native.monoisotopic, 2.01565, 0.0001);

        // Reduction adds a methylation site as well as two hydrogens.
        let both = glycan_mass(
            &graph,
            &MassOptions {
                reducing_end: ReducingEnd::Reduced,
                permethylated: true,
            },
        )
        .unwrap();
        close(
            both.monoisotopic - methylated.monoisotopic,
            2.01565 + 14.01565,
            0.0001,
        );
    }

    #[test]
    fn labels_add_their_reductive_amination_increment() {
        let graph = parse_wurcs(MAN5).unwrap();
        let native = glycan_mass(&graph, &free()).unwrap().monoisotopic;
        for (reducing_end, increment) in [
            (ReducingEnd::TwoAb, 120.06875),
            (ReducingEnd::TwoAa, 121.05276),
            (ReducingEnd::Procainamide, 219.17355),
        ] {
            let options = MassOptions {
                reducing_end,
                ..free()
            };
            let labelled = glycan_mass(&graph, &options).unwrap().monoisotopic;
            close(labelled - native, increment, 0.0001);
            assert!(matches!(
                glycan_mass(
                    &graph,
                    &MassOptions {
                        permethylated: true,
                        ..options
                    }
                ),
                Err(MassError::UnsupportedPermethylation)
            ));
        }
    }

    #[test]
    fn a_glycan_without_a_free_reducing_end_cannot_be_derivatized() {
        let reduced = MassOptions {
            reducing_end: ReducingEnd::Reduced,
            ..free()
        };
        // Trehalose: both anomeric carbons are in the bond.
        let trehalose = parse_wurcs("WURCS=2.0/1,2,1/[a2122h-1a_1-5]/1-1/a1-b1").unwrap();
        assert!(glycan_mass(&trehalose, &free()).is_ok());
        assert!(matches!(
            glycan_mass(&trehalose, &reduced),
            Err(MassError::NoReducingEnd)
        ));
        // An alditol has been reduced already.
        let alditol = parse_wurcs("WURCS=2.0/1,1,0/[h2122h]/1/").unwrap();
        assert!(matches!(
            glycan_mass(&alditol, &reduced),
            Err(MassError::NoReducingEnd)
        ));
    }

    #[test]
    fn sulfate_keeps_its_mass_but_stops_permethylation() {
        let graph =
            parse_wurcs("WURCS=2.0/1,1,0/[a2122h-1b_1-5_2*NCC/3=O_6*OSO/3=O/3=O]/1/").unwrap();
        let mass = glycan_mass(&graph, &free()).unwrap();
        assert_eq!(mass.formula.to_string(), "C8H15NO9S");
        assert!(matches!(
            glycan_mass(
                &graph,
                &MassOptions {
                    permethylated: true,
                    ..free()
                }
            ),
            Err(MassError::UnsupportedPermethylation)
        ));
    }

    #[test]
    fn uncertain_linkages_keep_their_mass() {
        let exact = parse_wurcs("WURCS=2.0/2,2,1/[u2122h][a2112h-1b_1-5]/1-2/a4-b1").unwrap();
        let unknown_position =
            parse_wurcs("WURCS=2.0/2,2,1/[u2122h][a2112h-1x_1-5]/1-2/a?-b1").unwrap();
        let alternatives =
            parse_wurcs("WURCS=2.0/2,2,1/[u2122h][a2112h-1b_1-5]/1-2/a3|a4-b1").unwrap();
        let expected = glycan_formula(&exact, &free()).unwrap();
        assert_eq!(
            glycan_formula(&unknown_position, &free()).unwrap(),
            expected
        );
        assert_eq!(glycan_formula(&alternatives, &free()).unwrap(), expected);
    }

    #[test]
    fn compositions_and_undefined_fragments_have_a_mass() {
        // Hex2HexNAc2: the residues of a glycan, without its linkages.
        let composition =
            parse_wurcs("WURCS=2.0/2,4,0+/[uxxxxh_2*NCC/3=O][uxxxxh]/1-1-2-2/").unwrap();
        let mass = glycan_mass(&composition, &free()).unwrap();
        close(
            mass.monoisotopic,
            2.0 * 203.07937 + 2.0 * 162.05282 + 18.01056,
            0.001,
        );

        let fragment =
            parse_wurcs("WURCS=2.0/2,3,2/[u2122h][a2122h-1b_1-5]/1-2-2/a4-b1_c1-a3|b3}").unwrap();
        assert_eq!(fragment.undefined_linkages().len(), 1);
        assert_eq!(
            glycan_formula(&fragment, &free()).unwrap().to_string(),
            "C18H32O16"
        );

        // Lactose plus one acetyl, C2H2O, on either residue.
        let substituent =
            parse_wurcs("WURCS=2.0/2,2,2/[u2122h][a2112h-1b_1-5]/1-2/a4-b1_a?|b?}*OCC/3=O")
                .unwrap();
        assert_eq!(substituent.undefined_modifications().len(), 1);
        assert_eq!(
            glycan_formula(&substituent, &free()).unwrap().to_string(),
            "C14H24O12"
        );
    }

    #[test]
    fn residues_that_are_not_bonded_are_not_one_glycan() {
        let graph = parse_wurcs("WURCS=2.0/2,2,1/[u2122h][u2112h]/1-2/a?|b?}*OCC/3=O").unwrap();
        assert!(matches!(
            glycan_formula(&graph, &free()),
            Err(MassError::Disconnected)
        ));
    }

    #[test]
    fn repeats_scale_an_exact_count_and_refuse_a_range() {
        let exact = parse_wurcs("WURCS=2.0/2,2,2/[a2122h-1a_1-5][a2122h-1b_1-5]/1-2/a4-b1_a1-b4~3")
            .unwrap();
        // Three disaccharide units joined into one chain.
        assert_eq!(
            glycan_formula(&exact, &free()).unwrap(),
            Tally::of_molecule(&construct::construct_molecule(&exact).unwrap()).formula
        );
        assert_eq!(
            glycan_formula(&exact, &free()).unwrap().to_string(),
            "C36H62O31"
        );

        let variable =
            parse_wurcs("WURCS=2.0/2,2,2/[a2122h-1a_1-5][a2122h-1b_1-5]/1-2/a4-b1_a1-b4~n")
                .unwrap();
        assert!(matches!(
            glycan_formula(&variable, &free()),
            Err(MassError::VariableRepeat)
        ));
    }

    #[test]
    fn a_probability_is_an_ensemble_not_a_mass() {
        let graph = parse_wurcs("WURCS=2.0/2,2,1/[u2122h][a2122h-1b_1-5]/1-2/a4%.5%-b1").unwrap();
        assert!(matches!(
            glycan_formula(&graph, &free()),
            Err(MassError::Ensemble)
        ));
    }

    #[test]
    fn residues_without_chemistry_have_no_mass() {
        let mut residue = crabwurcs_core::residue_from_kind(ResidueKind::Hex).unwrap();
        residue.residue_kind = None;
        residue.display_name = Some("Foo".into());
        let mut graph = ResidueGraph::new();
        graph.add_residue(residue);
        assert!(matches!(
            glycan_formula(&graph, &free()),
            Err(MassError::UndefinedResidue(name)) if name == "Foo"
        ));
        assert!(matches!(
            glycan_formula(&ResidueGraph::new(), &free()),
            Err(MassError::Empty)
        ));
    }

    /// Every record the molecular constructor can build has the formula of
    /// that molecule, and the corpus SMILES agree with both.
    #[test]
    fn the_formula_of_every_corpus_glycan_is_the_formula_of_its_molecule() {
        let mut checked = 0usize;
        for line in crate::SOURCE_CORPUS
            .lines()
            .chain(crate::DERIVED_CORPUS.lines())
        {
            let mut fields = line.split('\t');
            let wurcs = fields.next().unwrap();
            let smiles = fields.nth(3).filter(|smiles| !smiles.is_empty());
            let graph = parse_wurcs(wurcs).unwrap();
            let formula =
                glycan_formula(&graph, &free()).unwrap_or_else(|error| panic!("{wurcs}: {error}"));
            if let Ok(molecule) = construct::construct_molecule(&graph) {
                assert_eq!(formula, Tally::of_molecule(&molecule).formula, "{wurcs}");
                checked += 1;
            }
            if let Some(smiles) = smiles {
                let molecule = crate::parse_smiles(smiles).unwrap();
                assert_eq!(formula, Tally::of_molecule(&molecule).formula, "{smiles}");
            }
        }
        assert!(checked > 800, "only {checked} records were constructed");
    }

    #[test]
    fn compositions_read_like_glycomics() {
        let graph = parse_wurcs(MAN5).unwrap();
        let composition = glycan_composition(&graph);
        assert_eq!(composition.to_string(), "HexNAc2Hex5");
        assert_eq!(composition.count("Hex"), 5);
        assert_eq!(composition.count("Fuc"), 0);

        let sulfated = parse_wurcs(
            "WURCS=2.0/2,2,1/[u2122h_2*NSO/3=O/3=O_6*OSO/3=O/3=O][a2121A-1a_1-5_2*OSO/3=O/3=O]/1-2/a4-b1",
        )
        .unwrap();
        assert_eq!(glycan_composition(&sulfated).to_string(), "HexA1HexN1S3");

        let sialylated = parse_wurcs(
            "WURCS=2.0/4,5,4/[u2122h_2*NCC/3=O][a2112h-1b_1-5][a1221m-1a_1-5][Aad21122h-2a_2-6_5*NCC/3=O]/1-2-3-4-4/a4-b1_a3-c1_b3-d2_b6-e2",
        )
        .unwrap();
        assert_eq!(
            glycan_composition(&sialylated).to_string(),
            "HexNAc1Hex1Fuc1NeuAc2"
        );
    }

    #[test]
    fn adducts_parse_their_usual_spellings() {
        for (text, adduct) in [
            ("H", Adduct::Proton),
            ("[M+H]+", Adduct::Proton),
            ("Na", Adduct::Sodium),
            ("[M+Na]+", Adduct::Sodium),
            ("M+K", Adduct::Potassium),
            ("nh4", Adduct::Ammonium),
            ("Li+", Adduct::Lithium),
            ("[M-H]-", Adduct::Deprotonated),
            ("-H", Adduct::Deprotonated),
            ("[M+Cl]-", Adduct::Chloride),
            ("formate", Adduct::Formate),
            ("acetate", Adduct::Acetate),
        ] {
            assert_eq!(text.parse::<Adduct>().unwrap(), adduct, "{text}");
        }
        assert!("xenon".parse::<Adduct>().is_err());
        for adduct in [Adduct::Deprotonated, Adduct::Formate, Adduct::Acetate] {
            assert_eq!(adduct.as_str().parse::<Adduct>().unwrap(), adduct);
        }
    }

    #[test]
    fn reducing_ends_parse_their_names() {
        for reducing_end in [
            ReducingEnd::Free,
            ReducingEnd::Reduced,
            ReducingEnd::TwoAb,
            ReducingEnd::TwoAa,
            ReducingEnd::Procainamide,
        ] {
            assert_eq!(
                reducing_end.as_str().parse::<ReducingEnd>().unwrap(),
                reducing_end
            );
        }
        assert_eq!("2AB".parse::<ReducingEnd>().unwrap(), ReducingEnd::TwoAb);
        assert!("3-ab".parse::<ReducingEnd>().is_err());
    }

    /// Gal(b1-4)GlcNAc(b1-2)Man, the shape of an antenna.
    const LINEAR: &str =
        "WURCS=2.0/3,3,2/[u1122h][a2122h-1b_1-5_2*NCC/3=O][a2112h-1b_1-5]/1-2-3/a2-b1_b4-c1";

    fn ion<'a>(fragments: &'a [Fragment], label: &str) -> &'a Fragment {
        fragments
            .iter()
            .find(|fragment| fragment.label == label)
            .unwrap_or_else(|| panic!("no {label} in {fragments:?}"))
    }

    #[test]
    fn fragments_of_a_linear_glycan_number_from_their_own_end() {
        let graph = parse_wurcs(LINEAR).unwrap();
        let ions = glycan_fragments(&graph, &free(), Adduct::Proton).unwrap();
        assert_eq!(ions.len(), 8);

        // B1 is the Gal oxonium ion and B2 the LacNAc one (m/z 366, the
        // classic N-glycan marker).
        close(ion(&ions, "B1").mz, 163.0601, 0.0002);
        close(ion(&ions, "B2").mz, 366.1395, 0.0002);
        close(ion(&ions, "C1").mz, 181.0707, 0.0002);
        close(ion(&ions, "C2").mz, 384.1500, 0.0002);
        // Y1 keeps the reducing Man, Y2 the Man and the GlcNAc.
        close(ion(&ions, "Y1").mz, 181.0707, 0.0002);
        close(ion(&ions, "Y2").mz, 384.1500, 0.0002);
        close(ion(&ions, "Z1").mz, 163.0601, 0.0002);
        assert_eq!(ion(&ions, "B1").residues, BTreeSet::from([2]));
        assert_eq!(ion(&ions, "Y1").residues, BTreeSet::from([0]));
        assert!(ion(&ions, "Y1").kind.keeps_reducing_end());
        assert!(!ion(&ions, "B1").kind.keeps_reducing_end());
    }

    #[test]
    fn the_two_sides_of_a_cleavage_add_up_to_the_glycan() {
        let graph = parse_wurcs(MAN5).unwrap();
        let ions = glycan_fragments(&graph, &free(), Adduct::Proton).unwrap();
        let whole = glycan_mass(&graph, &free()).unwrap().monoisotopic;
        let water = Formula::of(&[("H", 2), ("O", 1)])
            .monoisotopic_mass()
            .unwrap();
        for edge in 0..graph.edge_count() {
            let of_kind = |kind| {
                ions.iter()
                    .find(|fragment| fragment.edge == edge && fragment.kind == kind)
                    .unwrap()
            };
            let (b, y) = (of_kind(IonKind::B), of_kind(IonKind::Y));
            close(b.monoisotopic + y.monoisotopic, whole, 0.0001);
            let c = of_kind(IonKind::C);
            close(c.monoisotopic - b.monoisotopic, water, 0.0001);
            assert_eq!(b.residues.len() + y.residues.len(), graph.node_count());
        }
    }

    #[test]
    fn branches_are_lettered_heaviest_first() {
        // The core Man's arms: α6 carries two Man, α3 is one Man.
        let graph = parse_wurcs(MAN5).unwrap();
        let ions = glycan_fragments(&graph, &free(), Adduct::Proton).unwrap();
        let heavy = ion(&ions, "B2\u{03B1}");
        let light = ion(&ions, "B1\u{03B2}");
        assert_eq!(heavy.residues.len(), 3);
        assert_eq!(light.residues.len(), 1);
        // Y counts from the reducing end whichever branch it cleaves.
        assert_eq!(ion(&ions, "Y3\u{03B2}").residues.len(), 6);
        // Two branch points on the way to a terminal Man of the α6 arm.
        assert!(
            ions.iter()
                .any(|fragment| fragment.label == "B1\u{03B1}\u{03B1}")
        );
    }

    #[test]
    fn permethylated_fragments_keep_a_free_hydroxyl_at_the_cleavage() {
        let graph = parse_wurcs(MAN5).unwrap();
        let options = MassOptions {
            permethylated: true,
            ..free()
        };
        let ions = glycan_fragments(&graph, &options, Adduct::Sodium).unwrap();
        // A terminal permethylated hexose is the B ion at 219 protonated.
        let protonated = glycan_fragments(&graph, &options, Adduct::Proton).unwrap();
        close(ion(&protonated, "B1\u{03B2}").mz, 219.1227, 0.0002);
        // B loses methanol and Y a methyl, so the two still add up to the
        // permethylated glycan.
        let whole = glycan_mass(&graph, &options).unwrap().monoisotopic;
        let b = ion(&ions, "B1\u{03B2}");
        let y = ions
            .iter()
            .find(|fragment| fragment.kind == IonKind::Y && fragment.edge == b.edge)
            .unwrap();
        close(b.monoisotopic + y.monoisotopic, whole, 0.0001);
    }

    #[test]
    fn fragments_need_one_defined_tree() {
        let composition =
            parse_wurcs("WURCS=2.0/2,4,0+/[uxxxxh_2*NCC/3=O][uxxxxh]/1-1-2-2/").unwrap();
        assert!(matches!(
            glycan_fragments(&composition, &free(), Adduct::Proton),
            Err(MassError::NotATree(_))
        ));
        let repeat =
            parse_wurcs("WURCS=2.0/2,2,2/[a2122h-1a_1-5][a2122h-1b_1-5]/1-2/a4-b1_a1-b4~3")
                .unwrap();
        assert!(matches!(
            glycan_fragments(&repeat, &free(), Adduct::Proton),
            Err(MassError::NotATree(_))
        ));
    }

    #[test]
    fn map_bridges_replace_two_hydroxyls() {
        // Two hexoses joined through a phosphodiester:
        // 2 C6H12O6 + H3PO4 - 2 H2O.
        let graph = parse_wurcs(
            "WURCS=2.0/2,2,1/[a2122h-1b_1-5][a2122h-1b_1-5]/1-2/a4n1-b1n2*1OPO*2/3O/3=O",
        )
        .unwrap();
        assert!(
            graph
                .inner()
                .edge_weights()
                .all(|edge| edge.map_code.is_some())
        );
        assert_eq!(
            glycan_formula(&graph, &free()).unwrap().to_string(),
            "C12H23O14P"
        );
    }

    #[test]
    fn alditols_have_no_anomeric_carbon_to_count() {
        for (wurcs, formula) in [
            ("WURCS=2.0/1,1,0/[h2122h]/1/", "C6H14O6"),
            ("WURCS=2.0/1,1,0/[hxh]/1/", "C3H8O3"),
        ] {
            assert_eq!(
                glycan_formula(&parse_wurcs(wurcs).unwrap(), &free())
                    .unwrap()
                    .to_string(),
                formula
            );
        }

        // Written as an alditol or reduced by the option, it is one glycan.
        let alditol = parse_wurcs("WURCS=2.0/2,2,1/[h2122h][a2122h-1b_1-5]/1-2/a4-b1").unwrap();
        let free_form = parse_wurcs("WURCS=2.0/2,2,1/[u2122h][a2122h-1b_1-5]/1-2/a4-b1").unwrap();
        let reduced = MassOptions {
            reducing_end: ReducingEnd::Reduced,
            ..free()
        };
        assert_eq!(
            glycan_formula(&alditol, &free()).unwrap(),
            glycan_formula(&free_form, &reduced).unwrap()
        );
        // The alditol's extra hydroxyl is a methylation site as well.
        let methylated = MassOptions {
            permethylated: true,
            ..free()
        };
        assert_eq!(
            glycan_formula(&alditol, &methylated).unwrap(),
            glycan_formula(
                &free_form,
                &MassOptions {
                    permethylated: true,
                    ..reduced
                }
            )
            .unwrap()
        );
    }

    #[test]
    fn phosphocholine_is_counted_as_its_inner_salt() {
        let plain = parse_wurcs("WURCS=2.0/1,1,0/[a2122h-1b_1-5]/1/").unwrap();
        let modified =
            parse_wurcs("WURCS=2.0/1,1,0/[a2122h-1b_1-5_6*OP^XOCCNC/7C/7C/3O/3=O]/1/").unwrap();
        let increment = glycan_mass(&modified, &free()).unwrap().monoisotopic
            - glycan_mass(&plain, &free()).unwrap().monoisotopic;
        // The published phosphocholine increment, C5H12NO3P.
        close(increment, 165.05548, 0.0001);
    }
}
