use crate::{
    AnomericSymbol, CarbonPosition, Linkage, Modification, Monosaccharide, ResidueGraph,
    ResidueKind, RingClosure, classify_residue,
};
use petgraph::Direction;
use petgraph::graph::{EdgeIndex, NodeIndex};
use petgraph::visit::EdgeRef;
use std::collections::{BTreeSet, HashMap, HashSet};
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
    validate_motif(motif)?;

    let motif_graph = motif.inner();
    let target_graph = target.inner();
    if target_graph.node_count() < motif_graph.node_count() {
        return Ok(Vec::new());
    }

    let mut candidates = HashMap::<NodeIndex, Vec<NodeIndex>>::new();
    for motif_node in motif_graph.node_indices() {
        let motif_residue = &motif_graph[motif_node];
        let compatible = target_graph
            .node_indices()
            .filter(|target_node| residue_matches(motif_residue, &target_graph[*target_node]))
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

    if !pattern_kind.is_generic() {
        if pattern.skeleton_code != candidate.skeleton_code
            || prefix_class(&pattern.anomeric_prefix) != prefix_class(&candidate.anomeric_prefix)
            || (pattern.ring != RingClosure::Unknown && pattern.ring != candidate.ring)
        {
            return false;
        }
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
    if pattern_positions
        .iter()
        .any(|position| *position == CarbonPosition(0))
    {
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
    use crate::{RepeatCount, UndefinedLinkage, residue_from_kind};

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
