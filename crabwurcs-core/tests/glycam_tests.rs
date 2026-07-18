#[cfg(test)]
mod glycam_tests {
    use crabwurcs_iupac::{
        parse_glycam, parse_iupac_condensed, write_glycam, write_iupac_condensed,
    };
    use serde::Deserialize;
    use std::collections::HashMap;

    #[derive(Debug, Deserialize)]
    struct Archetype {
        #[allow(dead_code)]
        #[serde(rename = "ID")]
        id: String,
        wurcs: Option<String>,
        iupac: Option<String>,
        #[allow(dead_code)]
        iupac_extended: Option<String>,
        glycam: Option<String>,
        smiles: Option<String>,
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
    fn test_glycam_parse_all() {
        let entries = load_entries();
        let mut tested = 0;
        let mut parse_ok = 0;
        let mut parse_fail = 0;
        let mut failures: Vec<(String, String)> = Vec::new();

        for (key, archetype) in &entries {
            let glycam = match &archetype.glycam {
                Some(g) => g,
                None => continue,
            };
            if glycam.is_empty() {
                continue;
            }
            tested += 1;

            match parse_glycam(glycam) {
                Ok(graph) => {
                    if graph.node_count() == 0 {
                        parse_fail += 1;
                        failures.push((key.clone(), format!("0 nodes: {}", glycam)));
                    } else {
                        parse_ok += 1;
                    }
                }
                Err(e) => {
                    parse_fail += 1;
                    failures.push((key.clone(), format!("{:?}: {}", e, glycam)));
                }
            }
        }

        let pass_rate = (parse_ok as f64 / tested as f64) * 100.0;
        println!("\n======= GLYCAM Parse Results =======");
        println!(
            "Total: {}, Parse OK: {}, Fail: {} ({:.1}%)",
            tested, parse_ok, parse_fail, pass_rate
        );

        if !failures.is_empty() {
            println!("\n--- Failures (first 30) ---");
            for (key, desc) in &failures[..std::cmp::min(30, failures.len())] {
                let short = if desc.len() > 100 { &desc[..100] } else { desc };
                println!("  {} | {}", key, short);
            }
            if failures.len() > 30 {
                println!("  ... and {} more", failures.len() - 30);
            }
        }

        assert!(
            pass_rate >= 50.0,
            "GLYCAM parse pass rate {:.1}% below 50%",
            pass_rate
        );
    }

    #[test]
    fn test_glycam_to_iupac_node_match() {
        let entries = load_entries();
        let mut tested = 0;
        let mut glycam_ok = 0;
        let mut iupac_ok = 0;
        let mut node_match = 0;
        let mut node_mismatch = 0;
        let mut mismatches: Vec<(String, usize, usize, String, String)> = Vec::new();

        for (key, archetype) in &entries {
            let glycam = match &archetype.glycam {
                Some(g) => g,
                None => continue,
            };
            let iupac = match &archetype.iupac {
                Some(i) => i,
                None => continue,
            };
            if glycam.is_empty() || iupac.is_empty() {
                continue;
            }
            tested += 1;

            let g_graph = match parse_glycam(glycam) {
                Ok(g) => {
                    glycam_ok += 1;
                    g
                }
                Err(_) => continue,
            };
            let i_graph = match parse_iupac_condensed(iupac) {
                Ok(g) => {
                    iupac_ok += 1;
                    g
                }
                Err(_) => continue,
            };

            if g_graph.node_count() == i_graph.node_count() {
                node_match += 1;
            } else {
                node_mismatch += 1;
                mismatches.push((
                    key.clone(),
                    g_graph.node_count(),
                    i_graph.node_count(),
                    glycam.clone(),
                    iupac.clone(),
                ));
            }
        }

        let compared = node_match + node_mismatch;
        let match_rate = if compared > 0 {
            (node_match as f64 / compared as f64) * 100.0
        } else {
            0.0
        };

        println!("\n======= GLYCAM→IUPAC Node Count Match =======");
        println!(
            "Tested: {}, GLYCAM parse ok: {}, IUPAC parse ok: {}",
            tested, glycam_ok, iupac_ok
        );
        println!(
            "Node count matches: {} / {} ({:.1}%)",
            node_match, compared, match_rate
        );

        if !mismatches.is_empty() {
            println!("\n--- Mismatches (first 30) ---");
            for (key, gn, im, glycam, iupac) in &mismatches[..std::cmp::min(30, mismatches.len())] {
                println!(
                    "  {} | GLYCAM→{} IUPAC→{} | G:{} | I:{}",
                    key, gn, im, glycam, iupac
                );
            }
        }
    }

    #[test]
    fn test_glycam_to_wurcs_node_match() {
        let entries = load_entries();
        use crabwurcs_core::parse_wurcs;
        let mut tested = 0;
        let mut glyph_ok = 0;
        let mut wurcs_ok = 0;
        let mut node_match = 0;
        let mut node_mismatch = 0;
        let mut mismatches: Vec<(String, usize, usize)> = Vec::new();

        for (key, archetype) in &entries {
            let glycam = match &archetype.glycam {
                Some(g) => g,
                None => continue,
            };
            let wurcs = match &archetype.wurcs {
                Some(w) => w,
                None => continue,
            };
            if glycam.is_empty() || wurcs.is_empty() {
                continue;
            }
            tested += 1;

            let g_graph = match parse_glycam(glycam) {
                Ok(g) => {
                    glyph_ok += 1;
                    g
                }
                Err(_) => continue,
            };
            let w_graph = match parse_wurcs(wurcs) {
                Ok(g) => {
                    wurcs_ok += 1;
                    g
                }
                Err(_) => continue,
            };

            if g_graph.node_count() == w_graph.node_count() {
                node_match += 1;
            } else {
                node_mismatch += 1;
                mismatches.push((key.clone(), g_graph.node_count(), w_graph.node_count()));
            }
        }

        let compared = node_match + node_mismatch;
        let match_rate = if compared > 0 {
            (node_match as f64 / compared as f64) * 100.0
        } else {
            0.0
        };

        println!("\n======= GLYCAM→WURCS Node Count Match =======");
        println!(
            "Tested: {}, GLYCAM ok: {}, WURCS ok: {}",
            tested, glyph_ok, wurcs_ok
        );
        println!(
            "Node count matches: {} / {} ({:.1}%)",
            node_match, compared, match_rate
        );

        if !mismatches.is_empty() {
            println!("\n--- Mismatches (first 20) ---");
            for (key, gn, wn) in &mismatches[..std::cmp::min(20, mismatches.len())] {
                println!("  {} | GLYCAM→{} WURCS→{}", key, gn, wn);
            }
        }
    }

    #[test]
    fn test_anomer_fidelity() {
        let entries = load_entries();
        let mut tested = 0;
        let mut anomer_matches = 0;
        let mut anomer_differs = 0;

        for (_key, archetype) in &entries {
            let iupac = match &archetype.iupac {
                Some(i) => i,
                None => continue,
            };
            if iupac.is_empty() {
                continue;
            }
            tested += 1;

            let alpha_count_orig = iupac.matches("(a").count();
            let beta_count_orig = iupac.matches("(b").count();

            let graph = match parse_iupac_condensed(iupac) {
                Ok(g) => g,
                Err(_) => continue,
            };
            let output = match write_iupac_condensed(&graph) {
                Ok(o) => o,
                Err(_) => continue,
            };

            let alpha_count_out = output.matches("(a").count();
            let beta_count_out = output.matches("(b").count();

            if alpha_count_orig == alpha_count_out && beta_count_orig == beta_count_out {
                anomer_matches += 1;
            } else {
                anomer_differs += 1;
            }
        }

        let match_rate = (anomer_matches as f64 / tested as f64) * 100.0;
        println!("\n======= Anomer Fidelity (alpha/beta counts) =======");
        println!(
            "Tested: {}, Match: {}, Differ: {} ({:.1}%)",
            tested, anomer_matches, anomer_differs, match_rate
        );
    }

    #[test]
    fn test_smiles_present() {
        let entries = load_entries();
        let mut with_smiles = 0;
        let mut without = 0;
        for (_, archetype) in &entries {
            if archetype.smiles.is_some() {
                with_smiles += 1;
            } else {
                without += 1;
            }
        }
        println!("\n======= SMILES Availability =======");
        println!("Entries with SMILES: {}, without: {}", with_smiles, without);
        assert!(with_smiles > 0, "No SMILES data found");
    }

    #[test]
    fn test_smiles_to_wurcs_skip_unless_rdkit() {
        use crabwurcs_core::parse_wurcs;
        let entries = load_entries();
        let mut tested = 0;
        let mut reached = false;

        for (_key, archetype) in &entries {
            let smiles = match &archetype.smiles {
                Some(s) => s,
                None => continue,
            };
            let wurcs = match &archetype.wurcs {
                Some(w) => w,
                None => continue,
            };
            if smiles.is_empty() || wurcs.is_empty() {
                continue;
            }

            // Without RDKit, we can't convert SMILES→WURCS
            // But we CAN test the inverse: WURCS→graph→write
            let graph = match parse_wurcs(wurcs) {
                Ok(g) => g,
                Err(_) => continue,
            };
            let _ = graph;
            reached = true;

            if tested > 100 {
                break;
            }
            tested += 1;
        }

        println!("\n======= SMILES→WURCS (requires RDKit) =======");
        println!(
            "Entries with SMILES+WURCS: tested {} samples (full test needs RDKit backend)",
            tested
        );
        assert!(reached, "No SMILES+WURCS entries found");
    }

    #[test]
    fn test_glycam_roundtrip() {
        let entries = load_entries();
        let mut tested = 0;
        let mut parse_ok = 0;
        let mut roundtrip_match = 0;
        let mut roundtrip_diff = 0;

        for (key, archetype) in &entries {
            let glycam = match &archetype.glycam {
                Some(g) => g,
                None => continue,
            };
            if glycam.is_empty() {
                continue;
            }
            tested += 1;

            let graph = match parse_glycam(glycam) {
                Ok(g) => {
                    parse_ok += 1;
                    g
                }
                Err(_) => continue,
            };
            let output = match write_glycam(&graph) {
                Ok(o) => o,
                Err(_) => continue,
            };

            if output.trim() == glycam.trim() {
                roundtrip_match += 1;
            } else {
                roundtrip_diff += 1;
                if roundtrip_diff <= 10 {
                    println!("GLYCAM RT diff {}: ORIG='{}' OUT='{}'", key, glycam, output);
                }
            }
        }

        let with_output = roundtrip_match + roundtrip_diff;
        let match_rate = if with_output > 0 {
            (roundtrip_match as f64 / with_output as f64) * 100.0
        } else {
            0.0
        };

        println!("\n======= GLYCAM Roundtrip =======");
        println!(
            "Total: {}, Parsed: {}, Match: {}, Diff: {} ({:.1}%)",
            tested, parse_ok, roundtrip_match, roundtrip_diff, match_rate
        );
    }
}
