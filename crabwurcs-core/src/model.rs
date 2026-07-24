use petgraph::graph::{Graph, NodeIndex};

#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct CarbonPosition(pub u8);

#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Stereo {
    D,
    L,
    Unspecified,
}

#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RingClosure {
    Pyranose,
    Furanose,
    Open,
    /// The notation does not declare a ring closure.
    Unknown,
}

#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Modification {
    pub position: CarbonPosition,
    pub descriptor: String,
    pub probability: Option<Probability>,
}

#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ProbabilityValue {
    Unknown,
    /// Fraction in ten-thousandths (`5000` = 0.5 = 50%).
    Known(u16),
}

#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Probability {
    pub lower: ProbabilityValue,
    pub upper: ProbabilityValue,
}

impl Probability {
    pub fn parse_wurcs(value: &str) -> Option<Self> {
        let (lower, upper) = value
            .split_once('-')
            .map_or((value, value), |(lower, upper)| (lower, upper));
        Some(Self {
            lower: parse_probability_value(lower)?,
            upper: parse_probability_value(upper)?,
        })
    }

    pub fn to_wurcs(self) -> String {
        let lower = probability_value_to_wurcs(self.lower);
        if self.lower == self.upper {
            lower
        } else {
            format!("{}-{}", lower, probability_value_to_wurcs(self.upper))
        }
    }

    pub fn to_iupac_percent(self) -> String {
        let display = |value| match value {
            ProbabilityValue::Unknown => "?".to_string(),
            ProbabilityValue::Known(value) => {
                let percent = value as f64 / 100.0;
                if percent.fract() == 0.0 {
                    format!("{}", percent as u16)
                } else {
                    format!("{percent:.2}").trim_end_matches('0').to_string()
                }
            }
        };
        if self.lower == self.upper {
            format!("{}%", display(self.lower))
        } else {
            format!("{},{}%", display(self.lower), display(self.upper))
        }
    }

    pub fn parse_iupac_percent(value: &str) -> Option<Self> {
        let value = value.strip_suffix('%')?;
        let (lower, upper) = value
            .split_once(',')
            .map_or((value, value), |(lower, upper)| (lower, upper));
        Some(Self {
            lower: parse_percent_value(lower)?,
            upper: parse_percent_value(upper)?,
        })
    }
}

fn parse_percent_value(value: &str) -> Option<ProbabilityValue> {
    if value == "?" {
        return Some(ProbabilityValue::Unknown);
    }
    let (whole, fraction) = value.split_once('.').unwrap_or((value, ""));
    let whole: u16 = whole.parse().ok()?;
    if whole > 100 {
        return None;
    }
    let mut fractional = fraction.chars().take(2).collect::<String>();
    while fractional.len() < 2 {
        fractional.push('0');
    }
    let fractional: u16 = if fractional.is_empty() {
        0
    } else {
        fractional.parse().ok()?
    };
    Some(ProbabilityValue::Known(whole * 100 + fractional))
}

fn parse_probability_value(value: &str) -> Option<ProbabilityValue> {
    if value == "?" {
        return Some(ProbabilityValue::Unknown);
    }
    let value = value.strip_prefix('0').unwrap_or(value);
    if value == "1" {
        return Some(ProbabilityValue::Known(10_000));
    }
    let fraction = value.strip_prefix('.').unwrap_or(value);
    if fraction.is_empty() || !fraction.chars().all(|character| character.is_ascii_digit()) {
        return None;
    }
    let mut padded = fraction.chars().take(4).collect::<String>();
    while padded.len() < 4 {
        padded.push('0');
    }
    padded.parse().ok().map(ProbabilityValue::Known)
}

fn probability_value_to_wurcs(value: ProbabilityValue) -> String {
    match value {
        ProbabilityValue::Unknown => "?".to_string(),
        ProbabilityValue::Known(10_000) => "1".to_string(),
        ProbabilityValue::Known(value) => {
            let mut fraction = format!("{value:04}");
            while fraction.ends_with('0') {
                fraction.pop();
            }
            format!(".{fraction}")
        }
    }
}

#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AnomericSymbol {
    Alpha,
    Beta,
    Unknown,
    OpenChain,
}

impl AnomericSymbol {
    pub fn from_char(c: char) -> Self {
        match c {
            'a' => AnomericSymbol::Alpha,
            'b' => AnomericSymbol::Beta,
            'o' => AnomericSymbol::OpenChain,
            _ => AnomericSymbol::Unknown,
        }
    }

    pub fn to_char(self) -> char {
        match self {
            AnomericSymbol::Alpha => 'a',
            AnomericSymbol::Beta => 'b',
            AnomericSymbol::OpenChain => 'o',
            AnomericSymbol::Unknown => 'x',
        }
    }
}

#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Monosaccharide {
    pub backbone_length: u8,
    pub skeleton_code: String,
    pub stereo: Vec<Stereo>,
    pub ring: RingClosure,
    pub ring_start: Option<u8>,
    pub ring_end: Option<u8>,
    pub anomeric_position: u8,
    pub anomeric_symbol: AnomericSymbol,
    pub anomeric_prefix: String,
    pub modifications: Vec<Modification>,
    /// Original notation name for an SNFG-assigned residue whose chemistry is
    /// not present in the official registry. Standard WURCS cannot serialize
    /// this value without losing its identity.
    #[cfg_attr(
        feature = "serde",
        serde(default, skip_serializing_if = "Option::is_none")
    )]
    pub display_name: Option<String>,
    /// Semantic registry identity retained while a graph remains in memory.
    /// This disambiguates notation-level classes that share a WURCS skeleton,
    /// such as `Sia` and generic `NulO`.
    #[cfg_attr(
        feature = "serde",
        serde(default, skip_serializing_if = "Option::is_none")
    )]
    pub residue_kind: Option<crate::ResidueKind>,
}

impl Monosaccharide {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        backbone_length: u8,
        skeleton_code: String,
        stereo: Vec<Stereo>,
        ring: RingClosure,
        ring_start: Option<u8>,
        ring_end: Option<u8>,
        anomeric_position: u8,
        anomeric_symbol: AnomericSymbol,
        anomeric_prefix: String,
        modifications: Vec<Modification>,
    ) -> Self {
        Self {
            backbone_length,
            skeleton_code,
            stereo,
            ring,
            ring_start,
            ring_end,
            anomeric_position,
            anomeric_symbol,
            anomeric_prefix,
            modifications,
            display_name: None,
            residue_kind: None,
        }
    }

    pub fn with_display_name(mut self, display_name: impl Into<String>) -> Self {
        self.display_name = Some(display_name.into());
        self
    }
}

#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Linkage {
    pub parent_position: CarbonPosition,
    pub child_position: CarbonPosition,
    pub parent_position_alternatives: Vec<CarbonPosition>,
    pub child_position_alternatives: Vec<CarbonPosition>,
    /// Repeat count for a WURCS repeat-closing edge (`~n`, `~3`, `~2-5`).
    pub repeat: Option<RepeatCount>,
    /// True when this edge closes a cycle rather than belonging to the
    /// rooted spanning tree.
    pub cyclic: bool,
    pub parent_probability: Option<Probability>,
    pub child_probability: Option<Probability>,
    /// WURCS MAP code for a substituent bridging the two backbones.
    pub map_code: Option<String>,
    /// WURCS direction descriptor on the parent-side endpoint (`n`, `u`,
    /// `d`, `e`, `z`, `x`, ...).
    pub parent_direction: Option<String>,
    /// MAP star/attachment index on the parent-side endpoint.
    pub parent_modification_position: Option<u8>,
    /// WURCS direction descriptor on the child-side endpoint.
    pub child_direction: Option<String>,
    /// MAP star/attachment index on the child-side endpoint.
    pub child_modification_position: Option<u8>,
}

#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum RepeatCount {
    Unknown,
    Exact(u32),
    Range { min: Option<u32>, max: Option<u32> },
}

impl RepeatCount {
    pub fn parse(value: &str) -> Option<Self> {
        let value = value.trim();
        if value.eq_ignore_ascii_case("n") || value == "?" {
            return Some(Self::Unknown);
        }
        if let Ok(count) = value.parse() {
            return Some(Self::Exact(count));
        }
        let (min, max) = value.split_once('-')?;
        let parse_bound = |bound: &str| match bound {
            "" | "?" | "n" | "N" => Some(None),
            number => number.parse::<u32>().ok().map(Some),
        };
        Some(Self::Range {
            min: parse_bound(min)?,
            max: parse_bound(max)?,
        })
    }

    pub fn to_wurcs(&self) -> String {
        match self {
            Self::Unknown => "n".to_string(),
            Self::Exact(count) => count.to_string(),
            Self::Range { min, max } => format!(
                "{}-{}",
                min.map(|value| value.to_string())
                    .unwrap_or_else(|| "n".to_string()),
                max.map(|value| value.to_string())
                    .unwrap_or_else(|| "n".to_string())
            ),
        }
    }
}

impl Linkage {
    pub fn new(parent_position: CarbonPosition, child_position: CarbonPosition) -> Self {
        Self {
            parent_position,
            child_position,
            parent_position_alternatives: Vec::new(),
            child_position_alternatives: Vec::new(),
            repeat: None,
            cyclic: false,
            parent_probability: None,
            child_probability: None,
            map_code: None,
            parent_direction: None,
            parent_modification_position: None,
            child_direction: None,
            child_modification_position: None,
        }
    }

    pub fn with_alternatives(
        parent_positions: Vec<CarbonPosition>,
        child_positions: Vec<CarbonPosition>,
    ) -> Self {
        let parent_position = parent_positions
            .first()
            .copied()
            .unwrap_or(CarbonPosition(0));
        let child_position = child_positions
            .first()
            .copied()
            .unwrap_or(CarbonPosition(0));
        Self {
            parent_position,
            child_position,
            parent_position_alternatives: parent_positions.into_iter().skip(1).collect(),
            child_position_alternatives: child_positions.into_iter().skip(1).collect(),
            repeat: None,
            cyclic: false,
            parent_probability: None,
            child_probability: None,
            map_code: None,
            parent_direction: None,
            parent_modification_position: None,
            child_direction: None,
            child_modification_position: None,
        }
    }

    pub fn parent_positions(&self) -> impl Iterator<Item = CarbonPosition> + '_ {
        std::iter::once(self.parent_position)
            .chain(self.parent_position_alternatives.iter().copied())
    }

    pub fn child_positions(&self) -> impl Iterator<Item = CarbonPosition> + '_ {
        std::iter::once(self.child_position).chain(self.child_position_alternatives.iter().copied())
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct UndefinedParent {
    pub residue: NodeIndex,
    pub positions: Vec<CarbonPosition>,
}

/// A fragment/antenna whose parent is one of several candidate residues.
/// This is deliberately not represented as several ordinary graph edges:
/// doing so would falsely assert that every candidate bond exists.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct UndefinedLinkage {
    pub child: NodeIndex,
    pub child_positions: Vec<CarbonPosition>,
    pub parents: Vec<UndefinedParent>,
}

/// A substituent whose attachment is one of several candidate residues,
/// represented by WURCS as `a?|b?}*MAP` without a child backbone.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct UndefinedModification {
    pub parents: Vec<UndefinedParent>,
    pub map_code: String,
}

#[derive(Debug, Clone)]
pub struct ResidueGraph {
    graph: Graph<Monosaccharide, Linkage>,
    root: Option<NodeIndex>,
    // Preserve a parsed WURCS record verbatim while the graph is unchanged,
    // including uncommon constructs not yet covered by the editable model.
    source_wurcs: Option<String>,
    source_iupac: Option<String>,
    source_iupac_extended: Option<String>,
    source_glycam: Option<String>,
    composition: bool,
    undefined_linkages: Vec<UndefinedLinkage>,
    undefined_modifications: Vec<UndefinedModification>,
}

impl ResidueGraph {
    pub fn new() -> Self {
        Self {
            graph: Graph::new(),
            root: None,
            source_wurcs: None,
            source_iupac: None,
            source_iupac_extended: None,
            source_glycam: None,
            composition: false,
            undefined_linkages: Vec::new(),
            undefined_modifications: Vec::new(),
        }
    }

    fn clear_sources(&mut self) {
        self.source_wurcs = None;
        self.source_iupac = None;
        self.source_iupac_extended = None;
        self.source_glycam = None;
    }

    pub fn add_residue(&mut self, residue: Monosaccharide) -> NodeIndex {
        self.clear_sources();
        let idx = self.graph.add_node(residue);
        if self.root.is_none() {
            self.root = Some(idx);
        }
        idx
    }

    pub fn add_linkage(&mut self, parent: NodeIndex, child: NodeIndex, linkage: Linkage) {
        self.clear_sources();
        self.graph.add_edge(parent, child, linkage);
    }

    pub fn add_undefined_linkage(&mut self, linkage: UndefinedLinkage) {
        self.clear_sources();
        self.undefined_linkages.push(linkage);
    }

    pub fn undefined_linkages(&self) -> &[UndefinedLinkage] {
        &self.undefined_linkages
    }

    pub fn undefined_linkages_mut(&mut self) -> &mut Vec<UndefinedLinkage> {
        self.clear_sources();
        &mut self.undefined_linkages
    }

    pub fn add_undefined_modification(&mut self, modification: UndefinedModification) {
        self.clear_sources();
        self.undefined_modifications.push(modification);
    }

    pub fn undefined_modifications(&self) -> &[UndefinedModification] {
        &self.undefined_modifications
    }

    pub fn undefined_modifications_mut(&mut self) -> &mut Vec<UndefinedModification> {
        self.clear_sources();
        &mut self.undefined_modifications
    }

    pub fn root(&self) -> Option<NodeIndex> {
        self.root
    }

    pub fn set_root(&mut self, root: NodeIndex) {
        self.clear_sources();
        self.root = Some(root);
    }

    pub fn residue(&self, idx: NodeIndex) -> Option<&Monosaccharide> {
        self.graph.node_weight(idx)
    }

    pub fn residue_mut(&mut self, idx: NodeIndex) -> Option<&mut Monosaccharide> {
        self.clear_sources();
        self.graph.node_weight_mut(idx)
    }

    pub fn inner(&self) -> &Graph<Monosaccharide, Linkage> {
        &self.graph
    }

    pub fn inner_mut(&mut self) -> &mut Graph<Monosaccharide, Linkage> {
        self.clear_sources();
        &mut self.graph
    }

    pub(crate) fn set_source_wurcs(&mut self, source: String) {
        self.source_wurcs = Some(source);
    }

    pub fn source_wurcs(&self) -> Option<&str> {
        self.source_wurcs.as_deref()
    }

    pub fn set_source_iupac(&mut self, source: String) {
        self.source_iupac = Some(source);
    }

    pub fn source_iupac(&self) -> Option<&str> {
        self.source_iupac.as_deref()
    }

    pub fn set_source_iupac_extended(&mut self, source: String) {
        self.source_iupac = None;
        self.source_iupac_extended = Some(source);
    }

    pub fn source_iupac_extended(&self) -> Option<&str> {
        self.source_iupac_extended.as_deref()
    }

    pub fn set_source_glycam(&mut self, source: String) {
        self.source_iupac = None;
        self.source_glycam = Some(source);
    }

    pub fn source_glycam(&self) -> Option<&str> {
        self.source_glycam.as_deref()
    }

    pub fn is_composition(&self) -> bool {
        self.composition
    }

    pub fn set_composition(&mut self, composition: bool) {
        self.clear_sources();
        self.composition = composition;
    }

    pub fn node_count(&self) -> usize {
        self.graph.node_count()
    }

    pub fn edge_count(&self) -> usize {
        self.graph.edge_count()
    }
}

impl Default for ResidueGraph {
    fn default() -> Self {
        Self::new()
    }
}
