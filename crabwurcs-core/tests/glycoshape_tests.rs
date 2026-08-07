#[cfg(test)]
mod glycoshape_tests {
    use crabwurcs_core::{parse_wurcs, write_wurcs, write_wurcs_canonical};
    use serde::Deserialize;
    use std::collections::HashMap;

    #[derive(Debug, Deserialize)]
    struct Archetype {
        #[allow(dead_code)]
        #[serde(rename = "ID")]
        id: String,
        #[allow(dead_code)]
        name: Option<String>,
        wurcs: Option<String>,
        #[allow(dead_code)]
        iupac: Option<String>,
        #[allow(dead_code)]
        iupac_extended: Option<String>,
        #[allow(dead_code)]
        glytoucan: Option<String>,
    }

    #[derive(Debug, Deserialize)]
    struct Entry {
        archetype: Archetype,
    }

    fn load_entries() -> Vec<(String, Archetype)> {
        let path = std::env::var("CRABWURCS_GLYCOSHAPE_JSON").unwrap_or_else(|_| {
            concat!(
                env!("CARGO_MANIFEST_DIR"),
                "/tests/fixtures/glycoshape.sample.json"
            )
            .to_owned()
        });
        let data = std::fs::read_to_string(path).expect("Cannot read GLYCOSHAPE.json");
        let raw: HashMap<String, Entry> =
            serde_json::from_str(&data).expect("Cannot parse GLYCOSHAPE.json");
        raw.into_iter().map(|(k, v)| (k, v.archetype)).collect()
    }

    #[test]
    fn test_glycoshape_wurcs_roundtrip() {
        let entries = load_entries();
        let total = entries.len();
        let mut passed = 0;
        let mut skipped = 0;
        let mut parse_failures = Vec::new();
        let mut write_failures = Vec::new();
        let mut roundtrip_differences = Vec::new();

        for (key, archetype) in &entries {
            let wurcs = match &archetype.wurcs {
                Some(w) => w,
                None => {
                    skipped += 1;
                    continue;
                }
            };

            let graph = match parse_wurcs(wurcs) {
                Ok(g) => g,
                Err(e) => {
                    parse_failures.push((key.clone(), wurcs.clone(), format!("{:?}", e)));
                    continue;
                }
            };

            let output = match write_wurcs(&graph) {
                Ok(o) => o,
                Err(e) => {
                    write_failures.push((key.clone(), wurcs.clone(), format!("{:?}", e)));
                    continue;
                }
            };

            if output.trim() != wurcs.trim() {
                roundtrip_differences.push((key.clone(), wurcs.clone(), output));
            } else {
                passed += 1;
            }
        }

        let parse_fail_count = parse_failures.len();
        let write_fail_count = write_failures.len();
        let diff_count = roundtrip_differences.len();
        let fail_count = parse_fail_count + write_fail_count + diff_count;
        let tested = passed + fail_count;

        if fail_count > 0 {
            let mut report = String::new();
            report.push_str("\n======= GLYCOSHAPE WURCS Roundtrip Results =======\n");
            report.push_str(&format!(
                "Total entries: {}, Skipped (no WURCS): {}, Tested: {}\n",
                total, skipped, tested
            ));
            report.push_str(&format!(
                "Passed: {}, Failed: {} ({:.1}% pass rate)\n",
                passed,
                fail_count,
                (passed as f64 / tested as f64) * 100.0
            ));
            report.push_str(&format!(
                "  Parse failures: {}, Write failures: {}, Roundtrip diffs: {}\n",
                parse_fail_count, write_fail_count, diff_count
            ));

            if parse_fail_count > 0 {
                report.push_str("\n--- Parse Failures ---\n");
                for (key, wurcs, err) in &parse_failures {
                    report.push_str(&format!("  {} | {} | {}\n", key, err, wurcs));
                }
            }

            if write_fail_count > 0 {
                report.push_str("\n--- Write Failures ---\n");
                for (key, wurcs, err) in &write_failures {
                    report.push_str(&format!("  {} | {} | {}\n", key, err, wurcs));
                }
            }

            if diff_count > 0 {
                report.push_str("\n--- Roundtrip Differences ---\n");
                for (key, original, output) in
                    &roundtrip_differences[..std::cmp::min(50, diff_count)]
                {
                    report.push_str(&format!(
                        "  {} | ORIG: {} | OUT:  {}\n",
                        key, original, output
                    ));
                }
                if diff_count > 50 {
                    report.push_str(&format!("  ... and {} more differences\n", diff_count - 50));
                }
            }

            println!("{}", report);

            let pass_rate = (passed as f64 / tested as f64) * 100.0;
            assert!(
                pass_rate >= 90.0,
                "GLYCOSHAPE WURCS roundtrip pass rate {:.1}% is below 90% threshold. {} passed, {} failed.",
                pass_rate,
                passed,
                fail_count
            );
        } else {
            println!(
                "GLYCOSHAPE WURCS Roundtrip: {}/{} entries passed (100%!)",
                passed, tested
            );
        }
    }

    #[test]
    fn test_glycoshape_wurcs_parse_only() {
        let entries = load_entries();
        let mut passed = 0;
        let mut skipped = 0;
        let mut failures: Vec<(String, String)> = Vec::new();

        for (key, archetype) in &entries {
            let wurcs = match &archetype.wurcs {
                Some(w) => w,
                None => {
                    skipped += 1;
                    continue;
                }
            };

            match parse_wurcs(wurcs) {
                Ok(graph) => {
                    if graph.node_count() == 0 {
                        failures.push((key.clone(), format!("parsed but got 0 nodes: {}", wurcs)));
                    } else {
                        passed += 1;
                    }
                }
                Err(_e) => {
                    failures.push((key.clone(), wurcs.clone()));
                }
            }
        }

        let tested = passed + failures.len();
        let pass_rate = (passed as f64 / tested as f64) * 100.0;
        println!(
            "GLYCOSHAPE WURCS Parse: {}/{} entries parsed ({:.1}%) [{} skipped no WURCS]",
            passed, tested, pass_rate, skipped
        );

        if !failures.is_empty() {
            println!("Parse failures ({}):", failures.len());
            for (key, wurcs) in &failures[..std::cmp::min(30, failures.len())] {
                println!("  {} | {}", key, wurcs);
            }
            if failures.len() > 30 {
                println!("  ... and {} more", failures.len() - 30);
            }
        }

        assert!(
            pass_rate >= 95.0,
            "GLYCOSHAPE WURCS parse pass rate {:.1}% is below 95% threshold",
            pass_rate
        );
    }

    #[test]
    fn test_glycoshape_wurcs_to_wurcs_stable() {
        let entries = load_entries();
        let total = entries.len();
        let mut parse_fail = 0;
        let mut stable = 0;
        let mut unstable = 0;

        for (_, archetype) in &entries {
            let wurcs = match &archetype.wurcs {
                Some(w) => w,
                None => continue,
            };

            let graph = match parse_wurcs(wurcs) {
                Ok(g) => g,
                Err(_) => {
                    parse_fail += 1;
                    continue;
                }
            };

            let round1 = match write_wurcs(&graph) {
                Ok(o) => o,
                Err(_) => {
                    parse_fail += 1;
                    continue;
                }
            };

            let graph2 = match parse_wurcs(&round1) {
                Ok(g) => g,
                Err(_) => {
                    unstable += 1;
                    continue;
                }
            };

            let round2 = match write_wurcs(&graph2) {
                Ok(o) => o,
                Err(_) => {
                    unstable += 1;
                    continue;
                }
            };

            if round1 == round2 {
                stable += 1;
            } else {
                unstable += 1;
            }
        }

        println!(
            "GLYCOSHAPE Stability: parse_fail={}, stable={}, unstable={}, total={}",
            parse_fail, stable, unstable, total
        );

        let covered = stable + unstable;
        if covered > 0 {
            let stability = (stable as f64 / covered as f64) * 100.0;
            println!(
                "Stability rate: {:.1}% (stable after one roundtrip)",
                stability
            );
            assert!(
                stability >= 95.0,
                "WURCS parse→write stability {:.1}% below 95%",
                stability
            );
        }
    }

    #[test]
    fn test_glycoshape_wurcs_samples() {
        let entries = load_entries();

        let known_ids = vec![
            "G97131OU", "G00027JG", "G00020MO", "G00003VQ", "G00019BE", "G00048ZA", "G00025YC",
            "G00173GD", "G00367NK", "G00393YJ", "G77095CK",
        ];

        let mut found = 0;
        for (_, archetype) in &entries {
            if let Some(ref gt) = archetype.glytoucan
                && known_ids.contains(&gt.as_str())
            {
                let wurcs = match &archetype.wurcs {
                    Some(w) => w,
                    None => continue,
                };
                let graph = parse_wurcs(wurcs)
                    .unwrap_or_else(|_| panic!("Failed to parse known glycan {}: {}", gt, wurcs));
                let output = write_wurcs(&graph)
                    .unwrap_or_else(|_| panic!("Failed to write known glycan {}", gt));
                assert!(!output.is_empty(), "Empty output for known glycan {}", gt);
                found += 1;
            }
        }

        println!(
            "Validated {}/{} known glycans in GLYCOSHAPE",
            found,
            known_ids.len()
        );
        assert!(
            found >= 1,
            "Only found {}/{} known glycans",
            found,
            known_ids.len()
        );
    }

    /// WURCS for the alpha/beta anomers must differ from the archetype at the
    /// reducing end and survive a parse→write round-trip. A reducing-end
    /// monosaccharide with a modification (e.g. `u2122h_2*NCC/3=O`) was being
    /// emitted with the anomeric+ring info spliced into the modification
    /// segment (`a2122h_2-1a_1-5*NCC/3=O`) instead of the header
    /// (`a2122h-1a_1-5_2*NCC/3=O`).
    #[derive(Debug, Deserialize)]
    struct Variant {
        #[allow(dead_code)]
        #[serde(rename = "ID")]
        id: Option<String>,
        wurcs: Option<String>,
    }

    #[derive(Debug, Deserialize)]
    struct TriEntry {
        archetype: Variant,
        alpha: Variant,
        beta: Variant,
    }

    #[test]
    #[ignore = "maintenance audit: set CRABWURCS_GLYCOSHAPE_JSON to the full corpus"]
    fn test_glycoshape_alpha_beta_roundtrip() {
        use crabwurcs_core::model::AnomericSymbol;
        let path = std::env::var("CRABWURCS_GLYCOSHAPE_JSON")
            .expect("set CRABWURCS_GLYCOSHAPE_JSON to run the full corpus audit");
        let data = std::fs::read_to_string(path).expect("Cannot read GLYCOSHAPE.json");
        let raw: std::collections::HashMap<String, TriEntry> =
            serde_json::from_str(&data).expect("Cannot parse GLYCOSHAPE.json");

        // Read the parsed reducing-end anomer (root residue) of a WURCS string.
        fn root_anomer(wurcs: &str) -> Option<AnomericSymbol> {
            let g = parse_wurcs(wurcs).ok()?;
            let root = g.root()?;
            Some(g.residue(root)?.anomeric_symbol)
        }

        fn roundtrip_root_anomer(wurcs: &str) -> Option<AnomericSymbol> {
            let graph = parse_wurcs(wurcs).ok()?;
            let canonical = write_wurcs_canonical(&graph).ok()?;
            root_anomer(&canonical)
        }

        let mut missing = Vec::new();
        // Cases where alpha/beta fail to encode the anomer at the reducing end:
        // crabWURCS parses them back as Unknown — so alpha, beta and the
        // archetype all collapse to the same (undefined-anomer) structure.
        let mut alpha_not_alpha = Vec::new();
        let mut beta_not_beta = Vec::new();
        let mut archetype_not_unknown = Vec::new();
        let mut alpha_eq_arch = Vec::new();
        let mut beta_eq_arch = Vec::new();
        let mut alpha_eq_beta = Vec::new();

        for (gid, e) in &raw {
            let aw = e.archetype.wurcs.as_deref();
            let al = e.alpha.wurcs.as_deref();
            let be = e.beta.wurcs.as_deref();
            if aw.is_none() || al.is_none() || be.is_none() {
                missing.push(gid.clone());
                continue;
            }
            let (aw, al, be) = (aw.unwrap(), al.unwrap(), be.unwrap());
            if al == aw {
                alpha_eq_arch.push(gid.clone());
            }
            if be == aw {
                beta_eq_arch.push(gid.clone());
            }
            if al == be {
                alpha_eq_beta.push(gid.clone());
            }
            if root_anomer(aw) != Some(AnomericSymbol::Unknown)
                || roundtrip_root_anomer(aw) != Some(AnomericSymbol::Unknown)
            {
                archetype_not_unknown.push(gid.clone());
            }
            if root_anomer(al) != Some(AnomericSymbol::Alpha)
                || roundtrip_root_anomer(al) != Some(AnomericSymbol::Alpha)
            {
                alpha_not_alpha.push(gid.clone());
            }
            if root_anomer(be) != Some(AnomericSymbol::Beta)
                || roundtrip_root_anomer(be) != Some(AnomericSymbol::Beta)
            {
                beta_not_beta.push(gid.clone());
            }
        }
        missing.sort();

        let mut report = String::new();
        report.push_str("\n======= alpha/beta reducing-end anomer audit =======\n");
        report.push_str(&format!(
            "entries: {}, missing-variant: {}\n",
            raw.len(),
            missing.len()
        ));
        report.push_str(&format!(
            "string equality -> alpha==archetype: {}, beta==archetype: {}, alpha==beta: {}\n",
            alpha_eq_arch.len(),
            beta_eq_arch.len(),
            alpha_eq_beta.len()
        ));
        report.push_str(&format!(
            "parsed reducing-end anomer wrong -> archetype not Unknown: {}, alpha not Alpha: {}, beta not Beta: {}\n",
            archetype_not_unknown.len(),
            alpha_not_alpha.len(),
            beta_not_beta.len()
        ));
        report.push_str(&format!(
            "ALPHA reducing-end anomer != Alpha ({}): {:?}\n",
            alpha_not_alpha.len(),
            alpha_not_alpha
        ));
        report.push_str(&format!(
            "BETA  reducing-end anomer != Beta ({}): {:?}\n",
            beta_not_beta.len(),
            beta_not_beta
        ));
        println!("{}", report);

        assert!(
            archetype_not_unknown.is_empty(),
            "archetype not Unknown: {:?}",
            archetype_not_unknown
        );
        assert_eq!(
            missing,
            ["GS00898"],
            "unexpected entries without alpha/beta WURCS"
        );
        assert!(alpha_eq_beta.is_empty(), "alpha==beta: {:?}", alpha_eq_beta);
        assert!(
            alpha_eq_arch.is_empty() && beta_eq_arch.is_empty(),
            "alpha==archetype: {:?}, beta==archetype: {:?}",
            alpha_eq_arch,
            beta_eq_arch
        );
        assert!(
            alpha_not_alpha.is_empty(),
            "{} alpha entries parse as non-Alpha reducing end: {:?}",
            alpha_not_alpha.len(),
            alpha_not_alpha
        );
        assert!(
            beta_not_beta.is_empty(),
            "{} beta entries parse as non-Beta reducing end: {:?}",
            beta_not_beta.len(),
            beta_not_beta
        );
    }
}
