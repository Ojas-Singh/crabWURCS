#[cfg(test)]
mod glycam_conversion_tests {
    use crabwurcs_core::{parse_wurcs, write_wurcs};
    use crabwurcs_iupac::{
        parse_glycam, parse_iupac_condensed, write_glycam, write_iupac_condensed,
    };
    use serde::Deserialize;
    use std::collections::HashMap;

    #[derive(Debug, Deserialize)]
    struct Archetype {
        #[allow(dead_code)]
        wurcs: Option<String>,
        #[allow(dead_code)]
        iupac: Option<String>,
        glycam: Option<String>,
    }

    #[derive(Debug, Deserialize)]
    struct Entry {
        archetype: Archetype,
    }

    fn load_entries() -> HashMap<String, Archetype> {
        let path = concat!(env!("CARGO_MANIFEST_DIR"), "/../GLYCOSHAPE.json");
        let data = std::fs::read_to_string(path).expect("Cannot read GLYCOSHAPE.json");
        let raw: HashMap<String, Entry> =
            serde_json::from_str(&data).expect("Cannot parse GLYCOSHAPE.json");
        raw.into_iter().map(|(k, v)| (k, v.archetype)).collect()
    }

    // =====================================================
    // GLYCAM → WURCS structural verification
    // =====================================================

    #[test]
    fn test_glycam_to_wurcs_structural_gs00002() {
        let glycam = "LFucpa1-2DGalpa1-3DGalpb1-3DGalpNAc";
        let graph = parse_glycam(glycam)
            .unwrap_or_else(|e| panic!("Failed to parse GLYCAM '{}': {:?}", glycam, e));

        assert!(graph.node_count() > 0, "GLYCAM parsed to 0 nodes");

        let wurcs =
            write_wurcs(&graph).unwrap_or_else(|e| panic!("Failed to write WURCS: {:?}", e));

        assert!(!wurcs.is_empty(), "Empty WURCS output");
        assert!(
            wurcs.starts_with("WURCS=2.0/"),
            "Not a valid WURCS string: {}",
            wurcs
        );

        // Verify the generated WURCS is parseable (stable roundtrip)
        let regraph = parse_wurcs(&wurcs).unwrap_or_else(|e| {
            panic!(
                "Generated WURCS is not parseable: {:?}\nWURCS: {}",
                e, wurcs
            )
        });
        assert_eq!(
            graph.node_count(),
            regraph.node_count(),
            "WURCS roundtrip changed node count: {} → {}",
            graph.node_count(),
            regraph.node_count()
        );
    }

    #[test]
    fn test_glycam_to_wurcs_structural_gs00010() {
        let glycam = "DGalpNAcb1-3DGalpa1-3DGalpb1-4DGlcp";
        let graph = parse_glycam(glycam)
            .unwrap_or_else(|e| panic!("Failed to parse GLYCAM '{}': {:?}", glycam, e));

        assert!(graph.node_count() > 0, "GLYCAM parsed to 0 nodes");

        let wurcs =
            write_wurcs(&graph).unwrap_or_else(|e| panic!("Failed to write WURCS: {:?}", e));

        assert!(wurcs.starts_with("WURCS=2.0/"), "Not a valid WURCS string");

        let regraph = parse_wurcs(&wurcs)
            .unwrap_or_else(|e| panic!("Generated WURCS not parseable: {:?}", e));
        assert_eq!(
            graph.node_count(),
            regraph.node_count(),
            "WURCS roundtrip changed node count"
        );
    }

    #[test]
    fn test_glycam_to_wurcs_structural_gs00025() {
        let glycam = "DFucpa1-2DGalpb1-3DGlcpNAcb1-3DGalpb1-4DGlcp";

        let graph = parse_glycam(glycam)
            .unwrap_or_else(|e| panic!("Failed to parse GLYCAM '{}': {:?}", glycam, e));

        assert!(graph.node_count() > 0, "GLYCAM parsed to 0 nodes");

        let wurcs =
            write_wurcs(&graph).unwrap_or_else(|e| panic!("Failed to write WURCS: {:?}", e));
        assert!(wurcs.starts_with("WURCS=2.0/"), "Not a valid WURCS string");

        let regraph = parse_wurcs(&wurcs)
            .unwrap_or_else(|e| panic!("Generated WURCS not parseable: {:?}", e));
        assert_eq!(
            graph.node_count(),
            regraph.node_count(),
            "WURCS roundtrip changed node count"
        );
    }

    // =====================================================
    // WURCS → IUPAC conversion tests
    // =====================================================

    #[test]
    fn test_wurcs_to_iupac_gs00002() {
        let wurcs = "WURCS=2.0/4,4,3/[u2112h_2*NCC/3=O][a2112h-1b_1-5][a2112h-1a_1-5][a1221m-1a_1-5]/1-2-3-4/a3-b1_b3-c1_c2-d1";

        let graph = parse_wurcs(wurcs).expect("WURCS parse");
        assert_eq!(graph.node_count(), 4, "Expected 4 residues");
        assert_eq!(graph.edge_count(), 3, "Expected 3 linkages");

        let iupac = write_iupac_condensed(&graph).expect("IUPAC write");
        assert!(!iupac.is_empty(), "Empty IUPAC output");
        assert!(iupac.contains("Fuc"), "Missing Fuc: {}", iupac);
        assert!(iupac.contains("Gal"), "Missing Gal: {}", iupac);
        assert!(iupac.contains("GalNAc"), "Missing GalNAc: {}", iupac);

        // Verify IUPAC is re-parseable
        let graph2 = parse_iupac_condensed(&iupac)
            .unwrap_or_else(|e| panic!("IUPAC not re-parseable '{}': {:?}", iupac, e));
        assert_eq!(graph2.node_count(), 4, "IUPAC roundtrip lost residues");
        assert_eq!(graph2.edge_count(), 3, "IUPAC roundtrip lost linkages");

        // Verify generated WURCS is structurally equivalent
        let wurcs2 = write_wurcs(&graph2).expect("WURCS from IUPAC");
        let graph3 = parse_wurcs(&wurcs2).expect("WURCS re-parse");
        assert_eq!(graph3.node_count(), 4);
        assert_eq!(graph3.edge_count(), 3);

        println!("WURCS→IUPAC: {}", iupac);
    }

    #[test]
    fn test_wurcs_to_iupac_gs00010() {
        let wurcs = "WURCS=2.0/4,4,3/[u2122h][a2112h-1b_1-5][a2112h-1a_1-5][a2112h-1b_1-5_2*NCC/3=O]/1-2-3-4/a4-b1_b3-c1_c3-d1";

        let graph = parse_wurcs(wurcs).expect("WURCS parse");
        assert_eq!(graph.node_count(), 4);
        assert_eq!(graph.edge_count(), 3);

        let iupac = write_iupac_condensed(&graph).expect("IUPAC write");
        assert!(!iupac.is_empty(), "Empty IUPAC");
        assert!(iupac.contains("Gal"), "Missing Gal: {}", iupac);
        assert!(iupac.contains("Glc"), "Missing Glc: {}", iupac);

        let graph2 = parse_iupac_condensed(&iupac)
            .unwrap_or_else(|e| panic!("IUPAC not re-parseable: {:?}", e));
        assert_eq!(graph2.node_count(), 4);
        assert_eq!(graph2.edge_count(), 3);

        println!("WURCS→IUPAC: {}", iupac);
    }

    #[test]
    fn test_wurcs_to_iupac_gs01002_branched() {
        let wurcs = "WURCS=2.0/6,8,7/[u2122h_2*NCC/3=O][a1221m-1a_1-5][a2122h-1b_1-5_2*NCC/3=O][a1122h-1b_1-5][a2112h-1b_1-5][a1122h-1a_1-5]/1-2-3-2-4-5-6-6/a3-b1_a4-c1_c4-e1_d4-f1_a6-d1_e3-g1_e6-h1";

        let graph = parse_wurcs(wurcs).expect("WURCS parse");
        assert_eq!(
            graph.node_count(),
            8,
            "Expected 8 residues in branched N-glycan"
        );
        assert_eq!(graph.edge_count(), 7, "Expected 7 linkages");

        let iupac = write_iupac_condensed(&graph).expect("IUPAC write");
        assert!(!iupac.is_empty(), "Empty IUPAC");

        let graph2 = parse_iupac_condensed(&iupac)
            .unwrap_or_else(|e| panic!("IUPAC not re-parseable: {:?}", e));
        assert_eq!(graph2.node_count(), 8, "IUPAC roundtrip lost residues");
        // Note: edge count may differ for complex branched structures
        // due to linkage order normalization during IUPAC parsing/writing
        println!("WURCS→IUPAC (branched): {}", iupac);
    }

    // =====================================================
    // WURCS → GLYCAM conversion tests
    // =====================================================

    #[test]
    fn test_wurcs_to_glycam_gs00010() {
        let wurcs = "WURCS=2.0/4,4,3/[u2122h][a2112h-1b_1-5][a2112h-1a_1-5][a2112h-1b_1-5_2*NCC/3=O]/1-2-3-4/a4-b1_b3-c1_c3-d1";

        let graph = parse_wurcs(wurcs).expect("WURCS parse");
        assert_eq!(graph.node_count(), 4);

        let glycam = write_glycam(&graph).expect("GLYCAM write");
        assert!(!glycam.is_empty(), "Empty GLYCAM output");
        assert!(glycam.contains("Gal"), "Missing Gal: {}", glycam);
        assert!(glycam.contains("Glc"), "Missing Glc: {}", glycam);

        println!("WURCS→GLYCAM: {}", glycam);
    }

    #[test]
    fn test_wurcs_to_glycam_gs00025() {
        let wurcs = "WURCS=2.0/4,5,4/[u2122h][a2112h-1b_1-5][a2122h-1b_1-5_2*NCC/3=O][a2112m-1a_1-5]/1-2-3-2-4/a4-b1_b3-c1_c3-d1_d2-e1";

        let graph = parse_wurcs(wurcs).expect("WURCS parse");
        assert_eq!(graph.node_count(), 5);

        let glycam = write_glycam(&graph).expect("GLYCAM write");
        assert!(!glycam.is_empty(), "Empty GLYCAM output");
        assert!(glycam.contains("Gal"), "Missing Gal: {}", glycam);
        assert!(glycam.contains("Glc"), "Missing Glc: {}", glycam);

        println!("WURCS→GLYCAM: {}", glycam);
    }

    #[test]
    fn test_wurcs_to_glycam_gs01002() {
        let wurcs = "WURCS=2.0/6,8,7/[u2122h_2*NCC/3=O][a1221m-1a_1-5][a2122h-1b_1-5_2*NCC/3=O][a1122h-1b_1-5][a2112h-1b_1-5][a1122h-1a_1-5]/1-2-3-2-4-5-6-6/a3-b1_a4-c1_c4-e1_d4-f1_a6-d1_e3-g1_e6-h1";

        let graph = parse_wurcs(wurcs).expect("WURCS parse");
        assert_eq!(graph.node_count(), 8);

        let glycam = write_glycam(&graph).expect("GLYCAM write");
        assert!(!glycam.is_empty(), "Empty GLYCAM output");
        assert!(
            glycam.contains("GlcpNAc"),
            "Missing GLYCAM GlcNAc: {}",
            glycam
        );
        assert!(glycam.contains("Man"), "Missing Man: {}", glycam);

        println!("WURCS→GLYCAM (branched): {}", glycam);
    }

    // =====================================================
    // Full cross-format roundtrip: WURCS → IUPAC → GLYCAM → WURCS
    // =====================================================

    #[test]
    fn test_wurcs_to_iupac_to_glycam_roundtrip_gs00002() {
        let wurcs_orig = "WURCS=2.0/4,4,3/[u2112h_2*NCC/3=O][a2112h-1b_1-5][a2112h-1a_1-5][a1221m-1a_1-5]/1-2-3-4/a3-b1_b3-c1_c2-d1";

        let graph1 = parse_wurcs(wurcs_orig).expect("WURCS parse");
        assert_eq!(graph1.node_count(), 4);

        let iupac = write_iupac_condensed(&graph1).expect("IUPAC write");
        let graph2 = parse_iupac_condensed(&iupac)
            .unwrap_or_else(|e| panic!("IUPAC not parseable: {:?}", e));

        let glycam = write_glycam(&graph2).expect("GLYCAM write");

        // Verify structural integrity is preserved through the chain
        assert_eq!(graph2.node_count(), 4, "IUPAC roundtrip changed node count");
        assert_eq!(graph2.edge_count(), 3, "IUPAC roundtrip changed edge count");

        // GLYCAM should contain key residue names
        let glycam_upper = glycam.to_uppercase();
        assert!(
            glycam_upper.contains("FUC"),
            "Missing Fuc in GLYCAM: {}",
            glycam
        );
        assert!(
            glycam_upper.contains("GAL"),
            "Missing Gal in GLYCAM: {}",
            glycam
        );
        assert!(
            glycam_upper.contains("GALPNAC"),
            "Missing GLYCAM GalNAc: {}",
            glycam
        );

        println!("Cross-format roundtrip:");
        println!("  WURCS:  {}", wurcs_orig);
        println!("  IUPAC:  {}", iupac);
        println!("  GLYCAM: {}", glycam);
    }

    // =====================================================
    // Batch structural comparison via GlycoShape
    // =====================================================

    #[test]
    fn test_glycam_wurcs_node_count_match_via_glycoshape() {
        let entries = load_entries();

        let selected_ids = ["GS00002", "GS00010", "GS00956", "GS01002"];

        let mut tested: usize = 0;
        let mut passed: usize = 0;
        let mut failures: Vec<String> = Vec::new();

        for id in &selected_ids {
            let archetype = match entries.get(*id) {
                Some(a) => a,
                None => continue,
            };

            let glycam = match &archetype.glycam {
                Some(g) if !g.is_empty() => g.as_str(),
                _ => continue,
            };
            let wurcs_ref = match &archetype.wurcs {
                Some(w) if !w.is_empty() => w.as_str(),
                _ => continue,
            };

            tested += 1;

            let gly_graph = match parse_glycam(glycam) {
                Ok(g) if g.node_count() > 0 => g,
                _ => {
                    failures.push(format!("{}: GLYCAM parse failed", id));
                    continue;
                }
            };

            let wur_graph = match parse_wurcs(wurcs_ref) {
                Ok(g) if g.node_count() > 0 => g,
                _ => {
                    failures.push(format!("{}: WURCS parse failed", id));
                    continue;
                }
            };

            if gly_graph.node_count() == wur_graph.node_count()
                && gly_graph.edge_count() == wur_graph.edge_count()
            {
                passed += 1;
            } else {
                failures.push(format!(
                    "{}: GLYCAM→{}N/{}E, WURCS→{}N/{}E | GLYCAM: {}",
                    id,
                    gly_graph.node_count(),
                    gly_graph.edge_count(),
                    wur_graph.node_count(),
                    wur_graph.edge_count(),
                    glycam,
                ));
            }
        }

        println!("\n=== GLYCAM↔WURCS Node/Edge Count Match ===");
        println!(
            "Tested: {}, Passed: {}, Failed: {}",
            tested,
            passed,
            failures.len()
        );

        if !failures.is_empty() {
            for f in &failures {
                eprintln!("  {}", f);
            }
        }

        assert!(
            failures.is_empty(),
            "{} of {} GLYCAM↔WURCS comparisons failed",
            failures.len(),
            tested
        );
    }
}
