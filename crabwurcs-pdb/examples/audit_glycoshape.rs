use crabwurcs_core::{Monosaccharide, ResidueGraph};
use crabwurcs_iupac::parse_iupac_condensed;
use crabwurcs_pdb::extract_glycans_from_file;
use petgraph::graph::NodeIndex;
use petgraph::visit::EdgeRef;
use serde::Deserialize;
use std::collections::{HashMap, HashSet};
use std::path::{Path, PathBuf};

#[derive(Deserialize)]
struct Entry {
    archetype: Archetype,
}

#[derive(Deserialize)]
struct Archetype {
    iupac: Option<String>,
}

fn residue_key(residue: &Monosaccharide, is_root: bool) -> String {
    let mut modifications = residue
        .modifications
        .iter()
        .map(|modification| {
            format!(
                "{}:{:?}:{}",
                modification.position.0, modification.probability, modification.descriptor
            )
        })
        .collect::<Vec<_>>();
    modifications.sort();
    format!(
        "{}:{:?}:{}:{}:{}",
        residue.skeleton_code,
        residue.ring,
        if is_root {
            "*".to_string()
        } else {
            residue.anomeric_position.to_string()
        },
        if is_root {
            '*'
        } else {
            residue.anomeric_symbol.to_char()
        },
        modifications.join("_")
    )
}

fn subtree(
    graph: &ResidueGraph,
    node: NodeIndex,
    root: NodeIndex,
    visited: &mut HashSet<NodeIndex>,
) -> String {
    if !visited.insert(node) {
        return "cycle".into();
    }
    let mut children = graph
        .inner()
        .edges(node)
        .map(|edge| {
            format!(
                "{}-{}:{}",
                edge.weight().child_position.0,
                edge.weight().parent_position.0,
                subtree(graph, edge.target(), root, visited)
            )
        })
        .collect::<Vec<_>>();
    children.sort();
    format!(
        "{}[{}]",
        residue_key(graph.residue(node).unwrap(), node == root),
        children.join(",")
    )
}

fn signature(graph: &ResidueGraph) -> String {
    let root = graph.root().expect("graph root");
    subtree(graph, root, root, &mut HashSet::new())
}

fn audit_file(id: &str, format: &str, path: &Path, expected: &ResidueGraph) -> Result<(), String> {
    let glycans = extract_glycans_from_file(path).map_err(|error| error.to_string())?;
    if glycans.len() != 1 {
        return Err(format!(
            "{id} {format}: expected one glycan, extracted {}",
            glycans.len()
        ));
    }
    let actual = &glycans[0].graph;
    let actual_signature = signature(actual);
    let expected_signature = signature(expected);
    if actual_signature == expected_signature {
        Ok(())
    } else {
        Err(format!(
            "{id} {format}: expected {} residues/{} links, extracted {}/{}\n  expected: {expected_signature}\n  actual:   {actual_signature}",
            expected.node_count(),
            expected.edge_count(),
            actual.node_count(),
            actual.edge_count()
        ))
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut arguments = std::env::args_os().skip(1);
    let json_path = PathBuf::from(arguments.next().unwrap_or_else(|| "GLYCOSHAPE.json".into()));
    let structure_directory = PathBuf::from(
        arguments
            .next()
            .unwrap_or_else(|| "glycoshape-structures".into()),
    );
    let entries: HashMap<String, Entry> =
        serde_json::from_str(&std::fs::read_to_string(json_path)?)?;

    let mut checked = 0usize;
    let mut mismatches = Vec::new();
    for (id, entry) in entries {
        let Some(iupac) = entry.archetype.iupac else {
            continue;
        };
        let expected = parse_iupac_condensed(&iupac)?;
        for format in ["PDB", "GLYCAM"] {
            let path = structure_directory.join(format!("{id}-{format}.pdb"));
            if !path.exists() {
                continue;
            }
            checked += 1;
            if let Err(error) = audit_file(&id, format, &path, &expected) {
                mismatches.push(error);
            }
        }
    }

    println!(
        "GlycoShape structure audit: {checked} checked, {} semantic mismatches",
        mismatches.len()
    );
    for mismatch in &mismatches {
        eprintln!("{mismatch}");
    }
    if mismatches.is_empty() {
        Ok(())
    } else {
        Err("structure audit failed".into())
    }
}
