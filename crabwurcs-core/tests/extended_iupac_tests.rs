use crabwurcs_core::parse_wurcs;
use crabwurcs_iupac::{parse_iupac_extended, write_iupac_extended};
use serde::Deserialize;
use std::collections::HashMap;

#[derive(Deserialize)]
struct Entry {
    archetype: Archetype,
}

#[derive(Deserialize)]
struct Archetype {
    iupac_extended: Option<String>,
    wurcs: Option<String>,
}

fn entries() -> HashMap<String, Entry> {
    serde_json::from_str(include_str!("../../GLYCOSHAPE.json")).unwrap()
}

#[test]
fn glycoshape_extended_iupac_is_lossless_and_structurally_complete() {
    let mut tested = 0;
    let mut failures = Vec::new();
    for (id, entry) in entries() {
        let (Some(extended), Some(wurcs)) = (entry.archetype.iupac_extended, entry.archetype.wurcs)
        else {
            continue;
        };
        tested += 1;
        match parse_iupac_extended(&extended) {
            Ok(graph) => {
                let expected = parse_wurcs(&wurcs).unwrap();
                if graph.node_count() != expected.node_count()
                    || graph.edge_count() != expected.edge_count()
                    || !matches!(write_iupac_extended(&graph), Ok(ref output) if output == &extended)
                {
                    failures.push(format!(
                        "{id}: graph {}/{} expected {}/{}",
                        graph.node_count(),
                        graph.edge_count(),
                        expected.node_count(),
                        expected.edge_count()
                    ));
                }
            }
            Err(error) => failures.push(format!("{id}: {error}")),
        }
    }
    assert_eq!(tested, 839);
    assert!(
        failures.is_empty(),
        "{} failures: {:?}",
        failures.len(),
        &failures[..failures.len().min(20)]
    );
}
