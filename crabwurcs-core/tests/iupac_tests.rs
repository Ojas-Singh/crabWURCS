#[cfg(test)]
mod iupac_tests {
    use crabwurcs_core::parse_wurcs;
    use crabwurcs_iupac::{parse_iupac_condensed, write_iupac_condensed};
    use serde::Deserialize;
    use std::collections::HashMap;

    #[derive(Debug, Deserialize)]
    struct Archetype {
        #[allow(dead_code)]
        #[serde(rename = "ID")]
        id: String,
        wurcs: Option<String>,
        iupac: Option<String>,
    }

    #[derive(Debug, Deserialize)]
    struct Entry {
        archetype: Archetype,
    }

    fn load_entries() -> Vec<(String, Archetype)> {
        let path = concat!(env!("CARGO_MANIFEST_DIR"), "/../GLYCOSHAPE.json");
        let data = std::fs::read_to_string(path).expect("Cannot read GLYCOSHAPE.json");
        let raw: HashMap<String, Entry> =
            serde_json::from_str(&data).expect("Cannot parse GLYCOSHAPE.json");
        raw.into_iter().map(|(k, v)| (k, v.archetype)).collect()
    }

    #[test]
    fn test_iupac_parse_on_all_entries() {
        let entries = load_entries();
        let mut with_iupac = 0;
        let mut parse_ok = 0;
        let mut parse_fail = 0;
        let mut failures: Vec<(String, String)> = Vec::new();

        for (key, archetype) in &entries {
            let iupac = match &archetype.iupac {
                Some(i) => i,
                None => continue,
            };
            with_iupac += 1;

            if iupac.is_empty() {
                continue;
            }

            match parse_iupac_condensed(iupac) {
                Ok(graph) => {
                    if graph.node_count() == 0 {
                        parse_fail += 1;
                        failures.push((key.clone(), format!("0 nodes: {}", iupac)));
                    } else {
                        parse_ok += 1;
                    }
                }
                Err(e) => {
                    parse_fail += 1;
                    failures.push((key.clone(), format!("{:?}: {}", e, iupac)));
                }
            }
        }

        let tested = parse_ok + parse_fail;
        let pass_rate = (parse_ok as f64 / tested as f64) * 100.0;
        println!("\n======= GLYCOSHAPE IUPAC Parse Results =======");
        println!(
            "Total with IUPAC: {}, Tested: {}, Parse OK: {}, Fail: {} ({:.1}%)",
            with_iupac, tested, parse_ok, parse_fail, pass_rate
        );

        if !failures.is_empty() {
            println!("\nParse failures (showing first 50):");
            for (key, desc) in &failures[..std::cmp::min(50, failures.len())] {
                let short = if desc.len() > 120 { &desc[..120] } else { desc };
                println!("  {} | {}", key, short);
            }
            if failures.len() > 50 {
                println!("  ... and {} more", failures.len() - 50);
            }

            // Categorize failures
            println!("\n--- Failure Categories ---");
            let unsupported = failures
                .iter()
                .filter(|(_, d)| d.contains("UnsupportedToken"))
                .count();
            let other = failures.len() - unsupported;
            println!("  UnsupportedToken: {}", unsupported);
            println!("  Other errors: {}", other);
        }

        assert!(
            pass_rate >= 50.0,
            "IUPAC parse pass rate {:.1}% is below 50% threshold. {} passed, {} failed.",
            pass_rate,
            parse_ok,
            parse_fail
        );
    }

    #[test]
    fn test_iupac_to_wurcs_via_graph() {
        let entries = load_entries();
        let mut tested = 0;
        let mut iupac_parse_fail = 0;
        let mut wurcs_parse_fail = 0;
        let mut graphs_match = 0;
        let mut graphs_differ = 0;
        let mut diffs: Vec<(String, String, String, usize, usize)> = Vec::new();

        for (key, archetype) in &entries {
            let iupac = match &archetype.iupac {
                Some(i) => i,
                None => continue,
            };
            let wurcs = match &archetype.wurcs {
                Some(w) => w,
                None => continue,
            };
            if iupac.is_empty() || wurcs.is_empty() {
                continue;
            }

            tested += 1;

            let iupac_graph = match parse_iupac_condensed(iupac) {
                Ok(g) => g,
                Err(_) => {
                    iupac_parse_fail += 1;
                    continue;
                }
            };

            let wurcs_graph = match parse_wurcs(wurcs) {
                Ok(g) => g,
                Err(_) => {
                    wurcs_parse_fail += 1;
                    continue;
                }
            };

            let iupac_nodes = iupac_graph.node_count();
            let wurcs_nodes = wurcs_graph.node_count();

            if iupac_nodes == wurcs_nodes {
                graphs_match += 1;
            } else {
                graphs_differ += 1;
                diffs.push((
                    key.clone(),
                    iupac.clone(),
                    wurcs.clone(),
                    iupac_nodes,
                    wurcs_nodes,
                ));
            }
        }

        let compared = graphs_match + graphs_differ;
        let match_rate = if compared > 0 {
            (graphs_match as f64 / compared as f64) * 100.0
        } else {
            0.0
        };

        println!("\n======= GLYCOSHAPE IUPAC→WURCS Graph Comparison =======");
        println!(
            "Entries with both IUPAC+WURCS: {}. IUPAC parse fail: {}, WURCS parse fail: {}",
            tested, iupac_parse_fail, wurcs_parse_fail
        );
        println!(
            "Node count matches: {} / {} ({:.1}%)",
            graphs_match, compared, match_rate
        );

        if !diffs.is_empty() {
            println!("\n--- Node count mismatches (first 30) ---");
            for (key, iupac, _wurcs, i_nodes, w_nodes) in &diffs[..std::cmp::min(30, diffs.len())] {
                let short_iupac = if iupac.len() > 80 {
                    &iupac[..80]
                } else {
                    iupac
                };
                println!(
                    "  {}: IUPAC→{} nodes, WURCS→{} nodes | {}",
                    key, i_nodes, w_nodes, short_iupac
                );
            }
            if diffs.len() > 30 {
                println!("  ... and {} more", diffs.len() - 30);
            }
        }
    }

    #[test]
    fn test_iupac_roundtrip() {
        let entries = load_entries();
        let mut tested = 0;
        let mut parse_ok = 0;
        let mut write_ok = 0;
        let mut roundtrip_match = 0;
        let mut roundtrip_diff = 0;
        let mut diffs: Vec<(String, String, String)> = Vec::new();

        for (key, archetype) in &entries {
            let iupac = match &archetype.iupac {
                Some(i) => i,
                None => continue,
            };
            if iupac.is_empty() {
                continue;
            }
            tested += 1;

            let graph = match parse_iupac_condensed(iupac) {
                Ok(g) => {
                    parse_ok += 1;
                    g
                }
                Err(_) => continue,
            };

            let output = match write_iupac_condensed(&graph) {
                Ok(o) => {
                    write_ok += 1;
                    o
                }
                Err(_) => continue,
            };

            if output.trim() == iupac.trim() {
                roundtrip_match += 1;
            } else {
                roundtrip_diff += 1;
                diffs.push((key.clone(), iupac.clone(), output));
            }
        }

        println!("\n======= GLYCOSHAPE IUPAC Roundtrip =======");
        println!(
            "Total IUPAC entries: {}, Parsed: {}, Written: {}",
            tested, parse_ok, write_ok
        );

        let with_output = roundtrip_match + roundtrip_diff;
        let match_rate = if with_output > 0 {
            (roundtrip_match as f64 / with_output as f64) * 100.0
        } else {
            0.0
        };

        println!(
            "Roundtrip match: {} / {} ({:.1}%)",
            roundtrip_match, with_output, match_rate
        );

        if !diffs.is_empty() {
            println!("\n--- Roundtrip differences (first 30) ---");
            for (key, orig, out) in &diffs[..std::cmp::min(30, diffs.len())] {
                let short_orig = if orig.len() > 80 { &orig[..80] } else { orig };
                let short_out = if out.len() > 80 { &out[..80] } else { out };
                println!("  {} | {}", key, short_orig);
                println!("      -> {}", short_out);
            }
            if diffs.len() > 30 {
                println!("  ... and {} more", diffs.len() - 30);
            }
        }
    }

    #[test]
    fn test_iupac_only_entries_parse() {
        let entries = load_entries();
        let mut iupac_only = 0;
        let mut parse_ok = 0;
        let mut parse_fail = 0;
        let mut failures: Vec<(String, String)> = Vec::new();

        for (key, archetype) in &entries {
            let iupac = match &archetype.iupac {
                Some(i) => i,
                None => continue,
            };

            if archetype.wurcs.is_some() {
                continue;
            }

            iupac_only += 1;

            if iupac.is_empty() {
                continue;
            }

            match parse_iupac_condensed(iupac) {
                Ok(graph) => {
                    if graph.node_count() == 0 {
                        parse_fail += 1;
                        failures.push((key.clone(), format!("0 nodes for {}", iupac)));
                    } else {
                        parse_ok += 1;
                    }
                }
                Err(e) => {
                    parse_fail += 1;
                    failures.push((key.clone(), format!("{:?}: {}", e, iupac)));
                }
            }
        }

        let pass_rate = (parse_ok as f64 / iupac_only as f64) * 100.0;
        println!("\n======= IUPAC-Only Entries ({}) =======", iupac_only);
        println!(
            "Parse OK: {}, Fail: {} ({:.1}%)",
            parse_ok, parse_fail, pass_rate
        );

        if !failures.is_empty() {
            println!("\n--- Failures (first 30) ---");
            for (key, desc) in &failures[..std::cmp::min(30, failures.len())] {
                let short = if desc.len() > 120 { &desc[..120] } else { desc };
                println!("  {} | {}", key, short);
            }
        }
    }
}
