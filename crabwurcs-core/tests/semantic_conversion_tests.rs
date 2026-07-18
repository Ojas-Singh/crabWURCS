use crabwurcs_core::{parse_wurcs, Monosaccharide, ResidueGraph};
use crabwurcs_iupac::parse_iupac_condensed;
use petgraph::graph::NodeIndex;
use petgraph::visit::EdgeRef;
use serde::Deserialize;
use std::collections::{HashMap, HashSet};

#[derive(Deserialize)]
struct Entry {
    archetype: Archetype,
}
#[derive(Deserialize)]
struct Archetype {
    iupac: Option<String>,
    wurcs: Option<String>,
}

fn residue_key(residue: &Monosaccharide, include_ring: bool) -> String {
    let skeleton = residue.skeleton_code.trim_end_matches(['h', 'm', 'x']);
    let mut mods: Vec<_> = residue
        .modifications
        .iter()
        .map(|m| format!("{}{:?}*{}", m.position.0, m.probability, m.descriptor))
        .collect();
    mods.sort();
    format!(
        "{}:{}:{}:{:?}:{}",
        skeleton,
        if include_ring {
            format!("{:?}", residue.ring)
        } else {
            "*".to_string()
        },
        residue.anomeric_position,
        residue.anomeric_symbol,
        mods.join("_")
    )
}

fn subtree(
    graph: &ResidueGraph,
    node: NodeIndex,
    visited: &mut HashSet<NodeIndex>,
    include_ring: bool,
) -> String {
    if !visited.insert(node) {
        return "cycle".into();
    }
    let mut children: Vec<_> = graph
        .inner()
        .edges(node)
        .map(|edge| {
            let child = edge
                .weight()
                .child_positions()
                .map(|p| p.0.to_string())
                .collect::<Vec<_>>()
                .join("/");
            let parent = edge
                .weight()
                .parent_positions()
                .map(|p| p.0.to_string())
                .collect::<Vec<_>>()
                .join("/");
            format!(
                "{}-{}:{:?}:{}:{:?}:{:?}:{:?}:{:?}:{:?}:{:?}:{:?}:{}",
                child,
                parent,
                edge.weight().repeat,
                edge.weight().cyclic,
                edge.weight().parent_probability,
                edge.weight().child_probability,
                edge.weight().map_code,
                edge.weight().parent_direction,
                edge.weight().parent_modification_position,
                edge.weight().child_direction,
                edge.weight().child_modification_position,
                subtree(graph, edge.target(), visited, include_ring)
            )
        })
        .collect();
    children.sort();
    format!(
        "{}[{}]",
        residue_key(graph.residue(node).unwrap(), include_ring),
        children.join(",")
    )
}

fn signature(graph: &ResidueGraph, include_ring: bool) -> String {
    let rooted = subtree(
        graph,
        graph.root().unwrap(),
        &mut HashSet::new(),
        include_ring,
    );
    let mut undefined = graph
        .undefined_linkages()
        .iter()
        .map(|linkage| {
            let child = graph
                .residue(linkage.child)
                .map(|residue| residue_key(residue, include_ring))
                .unwrap_or_default();
            let mut parents = linkage
                .parents
                .iter()
                .map(|parent| {
                    let residue = graph
                        .residue(parent.residue)
                        .map(|residue| residue_key(residue, include_ring))
                        .unwrap_or_default();
                    let positions = parent
                        .positions
                        .iter()
                        .map(|position| position.0.to_string())
                        .collect::<Vec<_>>()
                        .join("/");
                    format!("{residue}@{positions}")
                })
                .collect::<Vec<_>>();
            parents.sort();
            let child_positions = linkage
                .child_positions
                .iter()
                .map(|position| position.0.to_string())
                .collect::<Vec<_>>()
                .join("/");
            format!("{child}@{child_positions}->{}", parents.join("|"))
        })
        .collect::<Vec<_>>();
    undefined.sort();
    let mut undefined_modifications = graph
        .undefined_modifications()
        .iter()
        .map(|modification| {
            let mut parents = modification
                .parents
                .iter()
                .map(|parent| {
                    let residue = graph
                        .residue(parent.residue)
                        .map(|residue| residue_key(residue, include_ring))
                        .unwrap_or_default();
                    let positions = parent
                        .positions
                        .iter()
                        .map(|position| position.0.to_string())
                        .collect::<Vec<_>>()
                        .join("/");
                    format!("{residue}@{positions}")
                })
                .collect::<Vec<_>>();
            parents.sort();
            format!("{}->{}", parents.join("|"), modification.map_code)
        })
        .collect::<Vec<_>>();
    undefined_modifications.sort();
    format!(
        "composition={};{rooted};undefined={};undefined_modifications={}",
        graph.is_composition(),
        undefined.join(";"),
        undefined_modifications.join(";")
    )
}

#[test]
fn glycoshape_iupac_and_wurcs_are_semantically_equivalent() {
    let entries: HashMap<String, Entry> =
        serde_json::from_str(include_str!("../../GLYCOSHAPE.json")).unwrap();
    let mut mismatches = Vec::new();
    let mut information_limited = 0;
    let mut tested = 0;
    for (id, entry) in entries {
        let (Some(iupac), Some(wurcs)) = (entry.archetype.iupac, entry.archetype.wurcs) else {
            continue;
        };
        tested += 1;
        let from_iupac = parse_iupac_condensed(&iupac).unwrap();
        let from_wurcs = parse_wurcs(&wurcs).unwrap();
        if signature(&from_iupac, true) != signature(&from_wurcs, true) {
            let has_undeclared_ring = from_wurcs
                .inner()
                .node_weights()
                .any(|residue| residue.ring == crabwurcs_core::RingClosure::Unknown);
            if has_undeclared_ring && signature(&from_iupac, false) == signature(&from_wurcs, false)
            {
                information_limited += 1;
            } else {
                mismatches.push(id);
            }
        }
    }
    eprintln!(
        "semantic matches: {} exact, {} limited by undeclared WURCS ring, {} hard mismatches",
        tested - information_limited - mismatches.len(),
        information_limited,
        mismatches.len()
    );
    assert!(
        mismatches.is_empty(),
        "{} hard semantic mismatches: {:?}",
        mismatches.len(),
        &mismatches[..mismatches.len().min(30)]
    );
}
