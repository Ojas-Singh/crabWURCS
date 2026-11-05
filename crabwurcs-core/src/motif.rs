use crate::{
    AnomericSymbol, CarbonPosition, Linkage, Modification, Monosaccharide, ResidueGraph,
    ResidueKind, RingClosure, classify_residue,
};
use petgraph::Direction;
use petgraph::graph::{EdgeIndex, NodeIndex};
use petgraph::visit::EdgeRef;
use std::collections::{BTreeSet, HashMap, HashSet};
use std::str::FromStr;
use thiserror::Error;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MotifMatch {
    pub node_indices: BTreeSet<usize>,
    pub edge_indices: BTreeSet<usize>,
}

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum MotifError {
    #[error("motif is empty")]
    Empty,
    #[error("motif compositions are not supported")]
    Composition,
    #[error("motif must be a connected graph")]
    Disconnected,
    #[error("motif must be a directed tree")]
    NotTree,
    #[error("motif cycles and repeat closures are not supported")]
    CycleOrRepeat,
    #[error("motif contains an undefined linkage")]
    UndefinedLinkage,
    #[error("motif contains an undefined modification")]
    UndefinedModification,
    #[error("motif has {expected} residues but {found} positions were given")]
    PositionCount { expected: usize, found: usize },
}

/// Where a motif residue has to sit in the glycan it matches, as GlycoDraw's
/// `highlight_termini_list` (`terminal`, `internal`, `flexible`), plus the
/// reducing end.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub enum MotifPosition {
    /// Anywhere.
    #[default]
    Flexible,
    /// Nothing is bonded to it on the non-reducing side. A candidate parent
    /// of an undefined fragment may carry one, so it is not terminal.
    Terminal,
    /// At least one residue is bonded to it on the non-reducing side.
    Internal,
    /// The reducing-end residue of the glycan.
    Reducing,
}

impl MotifPosition {
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Flexible => "flexible",
            Self::Terminal => "terminal",
            Self::Internal => "internal",
            Self::Reducing => "reducing",
        }
    }

    fn holds(self, target: &ResidueGraph, node: NodeIndex) -> bool {
        match self {
            Self::Flexible => true,
            Self::Reducing => target.root() == Some(node),
            Self::Terminal | Self::Internal => {
                let has_child = target
                    .inner()
                    .edges_directed(node, Direction::Outgoing)
                    .any(|edge| edge.weight().repeat.is_none() && !edge.weight().cyclic);
                if self == Self::Internal {
                    return has_child;
                }
                !has_child
                    && !target.undefined_linkages().iter().any(|undefined| {
                        undefined
                            .parents
                            .iter()
                            .any(|parent| parent.residue == node)
                    })
            }
        }
    }
}

impl FromStr for MotifPosition {
    type Err = String;

    fn from_str(value: &str) -> Result<Self, Self::Err> {
        match value.trim().to_ascii_lowercase().as_str() {
            "flexible" | "any" => Ok(Self::Flexible),
            "terminal" => Ok(Self::Terminal),
            "internal" => Ok(Self::Internal),
            "reducing" | "root" => Ok(Self::Reducing),
            other => Err(format!(
                "unknown motif position: {other} (expected flexible, terminal, internal or reducing)"
            )),
        }
    }
}

/// A motif together with the positions its residues must take.
#[derive(Debug, Clone)]
pub struct MotifQuery {
    motif: ResidueGraph,
    /// By motif node index; empty when every residue is flexible.
    positions: Vec<MotifPosition>,
}

impl MotifQuery {
    /// A motif whose residues may sit anywhere.
    pub fn new(motif: ResidueGraph) -> Self {
        Self {
            motif,
            positions: Vec::new(),
        }
    }

    /// A motif with one position per residue, listed in the order the
    /// residues are read in its IUPAC or GLYCAM text, from the non-reducing
    /// end (the order of GlycoDraw's `highlight_termini_list`): for
    /// `Gal(b1-4)GlcNAc`, `[Terminal, Flexible]` asks for a terminal Gal.
    pub fn with_positions(
        motif: ResidueGraph,
        positions: &[MotifPosition],
    ) -> Result<Self, MotifError> {
        if positions.len() != motif.node_count() {
            return Err(MotifError::PositionCount {
                expected: motif.node_count(),
                found: positions.len(),
            });
        }
        Ok(Self {
            motif,
            // Notation readers number residues from the reducing end.
            positions: positions.iter().rev().copied().collect(),
        })
    }

    pub fn motif(&self) -> &ResidueGraph {
        &self.motif
    }

    /// Every occurrence in `target`; see [`find_motif_matches`].
    pub fn find(&self, target: &ResidueGraph) -> Result<Vec<MotifMatch>, MotifError> {
        find_matches(target, &self.motif, &self.positions)
    }

    /// Whether `target` holds the motif at least once.
    pub fn contained_in(&self, target: &ResidueGraph) -> Result<bool, MotifError> {
        Ok(!self.find(target)?.is_empty())
    }

    /// The number of occurrences in `target`. Occurrences that use the same
    /// residues and bonds count once.
    pub fn count_in(&self, target: &ResidueGraph) -> Result<usize, MotifError> {
        Ok(self.find(target)?.len())
    }
}

/// Find every injective, directed, non-induced occurrence of `motif` in
/// `target`.
///
/// Extra branches and extra residue modifications in the target are allowed.
/// Unknown motif anomers and carbon positions are wildcards. The returned
/// edge indices contain only motif edges, not boundary edges between a match
/// and the rest of the target graph.
pub fn find_motif_matches(
    target: &ResidueGraph,
    motif: &ResidueGraph,
) -> Result<Vec<MotifMatch>, MotifError> {
    find_matches(target, motif, &[])
}

fn find_matches(
    target: &ResidueGraph,
    motif: &ResidueGraph,
    positions: &[MotifPosition],
) -> Result<Vec<MotifMatch>, MotifError> {
    validate_motif(motif)?;

    let motif_graph = motif.inner();
    let target_graph = target.inner();
    if target_graph.node_count() < motif_graph.node_count() {
        return Ok(Vec::new());
    }

    let mut candidates = HashMap::<NodeIndex, Vec<NodeIndex>>::new();
    for motif_node in motif_graph.node_indices() {
        let motif_residue = &motif_graph[motif_node];
        let position = positions
            .get(motif_node.index())
            .copied()
            .unwrap_or_default();
        let compatible = target_graph
            .node_indices()
            .filter(|target_node| {
                residue_matches(motif_residue, &target_graph[*target_node])
                    && position.holds(target, *target_node)
            })
            .collect::<Vec<_>>();
        if compatible.is_empty() {
            return Ok(Vec::new());
        }
        candidates.insert(motif_node, compatible);
    }

    let mut order = motif_graph.node_indices().collect::<Vec<_>>();
    order.sort_by_key(|node| {
        (
            candidates[node].len(),
            std::cmp::Reverse(motif_graph.neighbors_undirected(*node).count()),
            node.index(),
        )
    });

    let mut matches = Vec::new();
    search_node_mappings(
        target,
        motif,
        &order,
        &candidates,
        0,
        &mut HashMap::new(),
        &mut HashSet::new(),
        &mut matches,
    );
    Ok(matches)
}

fn validate_motif(motif: &ResidueGraph) -> Result<(), MotifError> {
    let graph = motif.inner();
    if graph.node_count() == 0 {
        return Err(MotifError::Empty);
    }
    if motif.is_composition() {
        return Err(MotifError::Composition);
    }
    if !motif.undefined_linkages().is_empty() {
        return Err(MotifError::UndefinedLinkage);
    }
    if !motif.undefined_modifications().is_empty() {
        return Err(MotifError::UndefinedModification);
    }
    if petgraph::algo::is_cyclic_directed(graph)
        || graph
            .edge_weights()
            .any(|linkage| linkage.cyclic || linkage.repeat.is_some())
    {
        return Err(MotifError::CycleOrRepeat);
    }
    if graph.node_count() == 1 {
        return Ok(());
    }
    let start = graph.node_indices().next().ok_or(MotifError::Empty)?;
    let mut visited = HashSet::new();
    let mut stack = vec![start];
    while let Some(node) = stack.pop() {
        if visited.insert(node) {
            stack.extend(graph.neighbors_undirected(node));
        }
    }
    if visited.len() != graph.node_count() {
        return Err(MotifError::Disconnected);
    }
    if graph.edge_count() != graph.node_count() - 1 {
        return Err(MotifError::NotTree);
    }

    let roots = graph
        .node_indices()
        .filter(|node| {
            graph
                .neighbors_directed(*node, Direction::Incoming)
                .next()
                .is_none()
        })
        .count();
    if roots != 1
        || graph
            .node_indices()
            .any(|node| graph.neighbors_directed(node, Direction::Incoming).count() > 1)
    {
        return Err(MotifError::NotTree);
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn search_node_mappings(
    target: &ResidueGraph,
    motif: &ResidueGraph,
    order: &[NodeIndex],
    candidates: &HashMap<NodeIndex, Vec<NodeIndex>>,
    depth: usize,
    mapping: &mut HashMap<NodeIndex, NodeIndex>,
    used_targets: &mut HashSet<NodeIndex>,
    matches: &mut Vec<MotifMatch>,
) {
    if depth == order.len() {
        add_edge_mappings(target, motif, mapping, matches);
        return;
    }

    let motif_node = order[depth];
    for &target_node in &candidates[&motif_node] {
        if used_targets.contains(&target_node) {
            continue;
        }
        mapping.insert(motif_node, target_node);
        if partial_edges_match(target, motif, mapping, motif_node) {
            used_targets.insert(target_node);
            search_node_mappings(
                target,
                motif,
                order,
                candidates,
                depth + 1,
                mapping,
                used_targets,
                matches,
            );
            used_targets.remove(&target_node);
        }
        mapping.remove(&motif_node);
    }
}

fn partial_edges_match(
    target: &ResidueGraph,
    motif: &ResidueGraph,
    mapping: &HashMap<NodeIndex, NodeIndex>,
    new_node: NodeIndex,
) -> bool {
    motif
        .inner()
        .edges_directed(new_node, Direction::Outgoing)
        .chain(motif.inner().edges_directed(new_node, Direction::Incoming))
        .all(|edge| {
            let Some(&source) = mapping.get(&edge.source()) else {
                return true;
            };
            let Some(&target_node) = mapping.get(&edge.target()) else {
                return true;
            };
            target
                .inner()
                .edges_connecting(source, target_node)
                .any(|candidate| linkage_matches(edge.weight(), candidate.weight()))
        })
}

fn add_edge_mappings(
    target: &ResidueGraph,
    motif: &ResidueGraph,
    mapping: &HashMap<NodeIndex, NodeIndex>,
    matches: &mut Vec<MotifMatch>,
) {
    let mut edge_candidates = Vec::<Vec<EdgeIndex>>::new();
    for motif_edge in motif.inner().edge_references() {
        let source = mapping[&motif_edge.source()];
        let target_node = mapping[&motif_edge.target()];
        let compatible = target
            .inner()
            .edges_connecting(source, target_node)
            .filter(|edge| linkage_matches(motif_edge.weight(), edge.weight()))
            .map(|edge| edge.id())
            .collect::<Vec<_>>();
        if compatible.is_empty() {
            return;
        }
        edge_candidates.push(compatible);
    }

    let nodes = mapping
        .values()
        .map(|node| node.index())
        .collect::<BTreeSet<_>>();
    choose_edges(&edge_candidates, 0, &mut BTreeSet::new(), &nodes, matches);
}

fn choose_edges(
    candidates: &[Vec<EdgeIndex>],
    depth: usize,
    chosen: &mut BTreeSet<usize>,
    nodes: &BTreeSet<usize>,
    matches: &mut Vec<MotifMatch>,
) {
    if depth == candidates.len() {
        let found = MotifMatch {
            node_indices: nodes.clone(),
            edge_indices: chosen.clone(),
        };
        if !matches.contains(&found) {
            matches.push(found);
        }
        return;
    }
    for edge in &candidates[depth] {
        if chosen.insert(edge.index()) {
            choose_edges(candidates, depth + 1, chosen, nodes, matches);
            chosen.remove(&edge.index());
        }
    }
}

fn residue_matches(pattern: &Monosaccharide, candidate: &Monosaccharide) -> bool {
    let Some(pattern_kind) = classify_residue(pattern) else {
        return false;
    };
    let Some(candidate_kind) = classify_residue(candidate) else {
        return false;
    };
    if !pattern_kind.matches_family(candidate_kind) {
        return false;
    }

    if pattern_kind == ResidueKind::Assigned
        && pattern.display_name.as_deref() != candidate.display_name.as_deref()
    {
        return false;
    }

    if !pattern_kind.is_generic()
        && (pattern.skeleton_code != candidate.skeleton_code
            || prefix_class(&pattern.anomeric_prefix) != prefix_class(&candidate.anomeric_prefix)
            || (pattern.ring != RingClosure::Unknown && pattern.ring != candidate.ring))
    {
        return false;
    }

    if pattern.anomeric_symbol != AnomericSymbol::Unknown
        && pattern.anomeric_symbol != candidate.anomeric_symbol
    {
        return false;
    }
    if pattern.anomeric_position != 0 && pattern.anomeric_position != candidate.anomeric_position {
        return false;
    }

    pattern.modifications.iter().all(|required| {
        candidate
            .modifications
            .iter()
            .any(|actual| modification_required(required, actual))
    })
}

fn prefix_class(value: &str) -> char {
    value
        .chars()
        .next()
        .filter(|character| matches!(character, 'A' | 'h'))
        .unwrap_or('a')
}

fn modification_required(required: &Modification, actual: &Modification) -> bool {
    (required.position == CarbonPosition(0) || required.position == actual.position)
        && required.descriptor == actual.descriptor
        && required
            .probability
            .is_none_or(|probability| actual.probability == Some(probability))
}

fn linkage_matches(pattern: &Linkage, candidate: &Linkage) -> bool {
    !candidate.cyclic
        && candidate.repeat.is_none()
        && position_matches(
            pattern.parent_position,
            &pattern.parent_position_alternatives,
            candidate.parent_position,
            &candidate.parent_position_alternatives,
        )
        && position_matches(
            pattern.child_position,
            &pattern.child_position_alternatives,
            candidate.child_position,
            &candidate.child_position_alternatives,
        )
        && pattern
            .map_code
            .as_ref()
            .is_none_or(|map| candidate.map_code.as_ref() == Some(map))
        && pattern
            .parent_probability
            .is_none_or(|probability| candidate.parent_probability == Some(probability))
        && pattern
            .child_probability
            .is_none_or(|probability| candidate.child_probability == Some(probability))
}

fn position_matches(
    pattern: CarbonPosition,
    pattern_alternatives: &[CarbonPosition],
    candidate: CarbonPosition,
    candidate_alternatives: &[CarbonPosition],
) -> bool {
    let pattern_positions = std::iter::once(pattern)
        .chain(pattern_alternatives.iter().copied())
        .collect::<Vec<_>>();
    if pattern_positions.contains(&CarbonPosition(0)) {
        return true;
    }
    std::iter::once(candidate)
        .chain(candidate_alternatives.iter().copied())
        .filter(|position| *position != CarbonPosition(0))
        .any(|candidate_position| pattern_positions.contains(&candidate_position))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{RepeatCount, UndefinedLinkage, UndefinedParent, residue_from_kind};

    fn one_residue(kind: ResidueKind) -> ResidueGraph {
        let mut graph = ResidueGraph::new();
        graph.add_residue(residue_from_kind(kind).unwrap());
        graph
    }

    #[test]
    fn generic_classes_expand_but_concrete_kinds_do_not() {
        let glcnac = one_residue(ResidueKind::GlcNAc);
        assert_eq!(
            find_motif_matches(&glcnac, &one_residue(ResidueKind::HexNAc))
                .unwrap()
                .len(),
            1
        );
        assert!(
            find_motif_matches(&glcnac, &one_residue(ResidueKind::HexN))
                .unwrap()
                .is_empty()
        );
        assert!(
            find_motif_matches(&glcnac, &one_residue(ResidueKind::GalNAc))
                .unwrap()
                .is_empty()
        );
    }

    #[test]
    fn motif_modifications_are_a_required_subset() {
        let mut target = one_residue(ResidueKind::GlcNAc);
        target
            .residue_mut(target.root().unwrap())
            .unwrap()
            .modifications
            .push(Modification {
                position: CarbonPosition(6),
                descriptor: "OSO/3=O/3=O".into(),
                probability: None,
            });
        assert_eq!(
            find_motif_matches(&target, &one_residue(ResidueKind::GlcNAc))
                .unwrap()
                .len(),
            1
        );

        let mut sulfated = one_residue(ResidueKind::GlcNAc);
        sulfated
            .residue_mut(sulfated.root().unwrap())
            .unwrap()
            .modifications
            .push(Modification {
                position: CarbonPosition(6),
                descriptor: "OSO/3=O/3=O".into(),
                probability: None,
            });
        assert_eq!(find_motif_matches(&target, &sulfated).unwrap().len(), 1);
        assert!(
            find_motif_matches(&one_residue(ResidueKind::GlcNAc), &sulfated)
                .unwrap()
                .is_empty()
        );
    }

    #[test]
    fn concrete_stereochemistry_must_agree() {
        let target = one_residue(ResidueKind::Glc);
        let mut inverted = one_residue(ResidueKind::Glc);
        inverted
            .residue_mut(inverted.root().unwrap())
            .unwrap()
            .skeleton_code = "a1211h".into();
        assert!(find_motif_matches(&target, &inverted).unwrap().is_empty());
    }

    /// `Gal(b1-4)GlcNAc`, with a sialic acid on the Gal when `capped`.
    fn lacnac(capped: bool) -> ResidueGraph {
        let mut graph = ResidueGraph::new();
        let glcnac = graph.add_residue(residue_from_kind(ResidueKind::GlcNAc).unwrap());
        let gal = graph.add_residue(residue_from_kind(ResidueKind::Gal).unwrap());
        graph.add_linkage(
            glcnac,
            gal,
            Linkage::new(CarbonPosition(4), CarbonPosition(1)),
        );
        if capped {
            let neu5ac = graph.add_residue(residue_from_kind(ResidueKind::Neu5Ac).unwrap());
            graph.add_linkage(
                gal,
                neu5ac,
                Linkage::new(CarbonPosition(3), CarbonPosition(2)),
            );
        }
        graph
    }

    fn lacnac_at(positions: &[MotifPosition]) -> MotifQuery {
        MotifQuery::with_positions(lacnac(false), positions).unwrap()
    }

    #[test]
    fn motif_positions_are_read_from_the_non_reducing_end() {
        use MotifPosition::*;
        let free = lacnac(false);
        let capped = lacnac(true);

        for (positions, on_free, on_capped) in [
            ([Flexible, Flexible], 1, 1),
            ([Terminal, Flexible], 1, 0),
            ([Internal, Flexible], 0, 1),
            ([Flexible, Reducing], 1, 1),
            ([Reducing, Flexible], 0, 0),
            ([Flexible, Terminal], 0, 0),
        ] {
            let query = lacnac_at(&positions);
            assert_eq!(query.count_in(&free).unwrap(), on_free, "{positions:?}");
            assert_eq!(query.count_in(&capped).unwrap(), on_capped, "{positions:?}");
            assert_eq!(query.contained_in(&free).unwrap(), on_free == 1);
        }
    }

    #[test]
    fn a_query_without_positions_finds_what_the_plain_matcher_finds() {
        let target = lacnac(true);
        let motif = lacnac(false);
        assert_eq!(
            MotifQuery::new(motif.clone()).find(&target).unwrap(),
            find_motif_matches(&target, &motif).unwrap()
        );
    }

    #[test]
    fn a_position_is_required_for_every_motif_residue() {
        assert_eq!(
            MotifQuery::with_positions(lacnac(false), &[MotifPosition::Terminal]).unwrap_err(),
            MotifError::PositionCount {
                expected: 2,
                found: 1
            }
        );
    }

    #[test]
    fn a_candidate_parent_of_an_undefined_fragment_is_not_terminal() {
        let mut target = lacnac(false);
        let gal = NodeIndex::new(1);
        let fuc = target.add_residue(residue_from_kind(ResidueKind::Fuc).unwrap());
        target.add_undefined_linkage(UndefinedLinkage {
            child: fuc,
            child_positions: vec![CarbonPosition(1)],
            parents: vec![UndefinedParent {
                residue: gal,
                positions: vec![CarbonPosition(0)],
            }],
        });
        assert_eq!(
            lacnac_at(&[MotifPosition::Terminal, MotifPosition::Flexible])
                .count_in(&target)
                .unwrap(),
            0
        );
    }

    #[test]
    fn motif_positions_parse_their_names() {
        for position in [
            MotifPosition::Flexible,
            MotifPosition::Terminal,
            MotifPosition::Internal,
            MotifPosition::Reducing,
        ] {
            assert_eq!(position.as_str().parse::<MotifPosition>(), Ok(position));
        }
        assert!("sideways".parse::<MotifPosition>().is_err());
    }

    #[test]
    fn unsupported_motif_topologies_return_precise_errors() {
        let mut composition = one_residue(ResidueKind::Hex);
        composition.set_composition(true);
        assert_eq!(
            find_motif_matches(&one_residue(ResidueKind::Hex), &composition),
            Err(MotifError::Composition)
        );

        let mut disconnected = ResidueGraph::new();
        disconnected.add_residue(residue_from_kind(ResidueKind::Glc).unwrap());
        disconnected.add_residue(residue_from_kind(ResidueKind::Gal).unwrap());
        assert_eq!(
            find_motif_matches(&disconnected, &disconnected),
            Err(MotifError::Disconnected)
        );

        let mut repeated = ResidueGraph::new();
        let parent = repeated.add_residue(residue_from_kind(ResidueKind::Glc).unwrap());
        let child = repeated.add_residue(residue_from_kind(ResidueKind::Gal).unwrap());
        let mut linkage = Linkage::new(CarbonPosition(4), CarbonPosition(1));
        linkage.repeat = Some(RepeatCount::Unknown);
        repeated.add_linkage(parent, child, linkage);
        assert_eq!(
            find_motif_matches(&repeated, &repeated),
            Err(MotifError::CycleOrRepeat)
        );

        let mut cycle = ResidueGraph::new();
        let first = cycle.add_residue(residue_from_kind(ResidueKind::Glc).unwrap());
        let second = cycle.add_residue(residue_from_kind(ResidueKind::Gal).unwrap());
        cycle.add_linkage(
            first,
            second,
            Linkage::new(CarbonPosition(4), CarbonPosition(1)),
        );
        cycle.add_linkage(
            second,
            first,
            Linkage::new(CarbonPosition(3), CarbonPosition(1)),
        );
        assert_eq!(
            find_motif_matches(&cycle, &cycle),
            Err(MotifError::CycleOrRepeat)
        );

        let mut undefined = one_residue(ResidueKind::Glc);
        undefined.add_undefined_linkage(UndefinedLinkage {
            child: undefined.root().unwrap(),
            child_positions: vec![CarbonPosition(1)],
            parents: Vec::new(),
        });
        assert_eq!(
            find_motif_matches(&undefined, &undefined),
            Err(MotifError::UndefinedLinkage)
        );
    }
}
