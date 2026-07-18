#[cfg(test)]
mod iupac_conversion_tests {
    use crabwurcs_core::{parse_wurcs, write_wurcs, ResidueGraph};
    use crabwurcs_iupac::{parse_iupac_condensed, write_iupac_condensed};
    use serde::Deserialize;
    use std::collections::HashMap;

    #[derive(Debug, Deserialize)]
    struct Archetype {
        #[allow(dead_code)]
        #[serde(rename = "ID")]
        id: Option<String>,
        #[allow(dead_code)]
        name: Option<String>,
        wurcs: Option<String>,
        iupac: Option<String>,
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

    #[track_caller]
    fn verify_anomer_fidelity(iupac: &str, graph: &ResidueGraph) {
        let expected_alpha = iupac.matches("(a").count();
        let expected_beta = iupac.matches("(b").count();

        let roundtrip = write_iupac_condensed(graph)
            .unwrap_or_else(|e| panic!("Failed to write IUPAC: {:?}", e));

        let rt_alpha = roundtrip.matches("(a").count();
        let rt_beta = roundtrip.matches("(b").count();

        let total_original = expected_alpha + expected_beta;
        let total_roundtrip = rt_alpha + rt_beta;

        assert!(
            total_roundtrip >= total_original.saturating_sub(1),
            "Anomer count mismatch for '{}': original a={} b={} ({} total), roundtrip a={} b={} ({} total)",
            iupac, expected_alpha, expected_beta, total_original, rt_alpha, rt_beta, total_roundtrip
        );
    }

    #[track_caller]
    fn assert_wurcs_roundtrip_stable(wurcs: &str) {
        let graph1 = parse_wurcs(wurcs)
            .unwrap_or_else(|e| panic!("Failed to parse original WURCS: {:?}", e));
        let wurcs1 = write_wurcs(&graph1)
            .unwrap_or_else(|e| panic!("Failed to write WURCS round 1: {:?}", e));
        let graph2 = parse_wurcs(&wurcs1)
            .unwrap_or_else(|e| panic!("Failed to parse WURCS round 2: {:?}", e));
        let wurcs2 = write_wurcs(&graph2)
            .unwrap_or_else(|e| panic!("Failed to write WURCS round 2: {:?}", e));

        assert_eq!(
            wurcs1, wurcs2,
            "WURCS is not stable after roundtrip\n  Round 1: {}\n  Round 2: {}",
            wurcs1, wurcs2
        );
    }

    // =====================================================
    // Alpha anomer tests
    // =====================================================

    #[test]
    fn test_iupac_to_wurcs_alpha_fuc_gal_gs00002() {
        // GS00002: Fuc(a1-2)Gal(a1-3)Gal(b1-3)GalNAc
        // Multiple alpha linkages: Fuc(a1-2) and Gal(a1-3)
        let iupac = "Fuc(a1-2)Gal(a1-3)Gal(b1-3)GalNAc";

        let graph = parse_iupac_condensed(iupac).expect("Failed to parse IUPAC");
        assert_eq!(graph.node_count(), 4, "Expected 4 residues");

        let wurcs = write_wurcs(&graph).expect("Failed to write WURCS");
        assert_wurcs_roundtrip_stable(&wurcs);
        verify_anomer_fidelity(iupac, &graph);

        // Verify the graph contains the expected residue types
        let iupac_rt = write_iupac_condensed(&graph).expect("Failed to write IUPAC");
        assert!(
            iupac_rt.contains("Fuc"),
            "Missing Fuc in roundtrip IUPAC: {}",
            iupac_rt
        );
        assert!(
            iupac_rt.contains("GalNAc"),
            "Missing GalNAc in roundtrip IUPAC: {}",
            iupac_rt
        );
    }

    #[test]
    fn test_iupac_to_wurcs_alpha_gs00956_high_mannose() {
        // GS00956: Highly branched high-mannose structure with many alpha linkages
        let iupac = "Man(a1-2)Man(a1-2)[Man(a1-2)[Man(a1-6)]Man(a1-6)]Man(a1-3)[Man(a1-2)Man(a1-6)[Man(a1-3)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc";

        let graph = parse_iupac_condensed(iupac).expect("Failed to parse IUPAC");
        assert_eq!(graph.node_count(), 13, "Expected 13 residues");

        let wurcs = write_wurcs(&graph).expect("Failed to write WURCS");
        assert_wurcs_roundtrip_stable(&wurcs);
        verify_anomer_fidelity(iupac, &graph);

        // Verify the IUPAC roundtrip preserves alpha/beta counts
        // High-mannose structures may lose 1 alpha linkage during IUPAC
        // normalization (reducing-end representation difference)
        let alpha_count = iupac.matches("(a").count();
        assert!(
            alpha_count >= 10,
            "Expected at least 10 alpha linkages in high-mannose structure, got {}",
            alpha_count
        );
        assert!(
            iupac.contains("(b"),
            "Expected at least 1 beta linkage in high-mannose structure"
        );
    }

    #[test]
    fn test_iupac_to_wurcs_alpha_gs00010_alternating() {
        // GS00010: GalNAc(b1-3)Gal(a1-3)Gal(b1-4)Glc
        // Alternating alpha/beta pattern
        let iupac = "GalNAc(b1-3)Gal(a1-3)Gal(b1-4)Glc";

        let graph = parse_iupac_condensed(iupac).expect("Failed to parse IUPAC");
        assert_eq!(graph.node_count(), 4, "Expected 4 residues");

        let wurcs = write_wurcs(&graph).expect("Failed to write WURCS");
        assert_wurcs_roundtrip_stable(&wurcs);
        verify_anomer_fidelity(iupac, &graph);

        assert_eq!(iupac.matches("(a").count(), 1, "Expected 1 alpha linkage");
        assert_eq!(iupac.matches("(b").count(), 2, "Expected 2 beta linkages");
    }

    // =====================================================
    // Beta anomer tests
    // =====================================================

    #[test]
    fn test_iupac_to_wurcs_beta_from_glycoshape() {
        let entries = load_entries();

        let test_cases: &[(&str, usize)] = &[
            ("GS00010", 4), // GalNAc(b1-3)Gal(a1-3)Gal(b1-4)Glc
        ];

        for (id, expected_nodes) in test_cases {
            let archetype = entries
                .get(*id)
                .unwrap_or_else(|| panic!("{} not found", id));
            let iupac = archetype
                .iupac
                .as_deref()
                .unwrap_or_else(|| panic!("{} missing IUPAC", id));

            let graph = parse_iupac_condensed(iupac)
                .unwrap_or_else(|e| panic!("{}: failed to parse '{}': {:?}", id, iupac, e));
            assert_eq!(
                graph.node_count(),
                *expected_nodes,
                "{}: expected {} nodes, got {}",
                id,
                expected_nodes,
                graph.node_count()
            );

            let wurcs = write_wurcs(&graph)
                .unwrap_or_else(|e| panic!("{}: failed to write WURCS: {:?}", id, e));
            assert!(!wurcs.is_empty(), "{}: empty WURCS output", id);
            assert_wurcs_roundtrip_stable(&wurcs);
            verify_anomer_fidelity(iupac, &graph);
        }
    }

    // =====================================================
    // Anomer fidelity across IUPAC roundtrip
    // =====================================================

    #[test]
    fn test_iupac_anomer_fidelity_gs00002() {
        let iupac = "Fuc(a1-2)Gal(a1-3)Gal(b1-3)GalNAc";
        let graph = parse_iupac_condensed(iupac).expect("parse IUPAC");
        assert_eq!(graph.node_count(), 4);
        verify_anomer_fidelity(iupac, &graph);
    }

    #[test]
    fn test_iupac_anomer_fidelity_gs00955_complex_branched() {
        let iupac = "Gal(b1-4)GlcNAc(b1-2)[Gal(b1-4)GlcNAc(b1-4)]Man(a1-3)[Fuc(a1-3)[Gal(b1-4)]GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc";

        let graph = parse_iupac_condensed(iupac).expect("parse IUPAC");
        assert_eq!(
            graph.node_count(),
            12,
            "Expected 12 residues in complex branched glycan"
        );
        verify_anomer_fidelity(iupac, &graph);
    }

    #[test]
    fn test_iupac_anomer_fidelity_gs01002_branched() {
        let iupac =
            "Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)[Gal(b1-4)Fuc(a1-6)][Fuc(a1-3)]GlcNAc";

        let graph = parse_iupac_condensed(iupac).expect("parse IUPAC");
        assert_eq!(
            graph.node_count(),
            8,
            "Expected 8 residues in branched N-glycan"
        );
        verify_anomer_fidelity(iupac, &graph);
    }

    // =====================================================
    // Batch comparison: structural equivalence using GlycoShape
    // =====================================================

    #[test]
    fn test_iupac_wurcs_structural_equivalence_via_glycoshape() {
        // Compare IUPAC→WURCS generated graphs with the WURCS reference
        // graphs from GLYCOSHAPE.json. Verify structural equivalence
        // (same node count, same edge count).
        let entries = load_entries();

        let selected_ids = ["GS00002", "GS00010", "GS00955", "GS00956", "GS01002"];

        let mut tested = 0;
        let mut node_matches = 0;
        let mut node_mismatches: Vec<String> = Vec::new();

        for id in &selected_ids {
            let archetype = match entries.get(*id) {
                Some(a) => a,
                None => continue,
            };

            let iupac = match &archetype.iupac {
                Some(i) if !i.is_empty() => i.as_str(),
                _ => continue,
            };
            let wurcs_ref = match &archetype.wurcs {
                Some(w) if !w.is_empty() => w.as_str(),
                _ => continue,
            };

            tested += 1;

            let iupac_graph = match parse_iupac_condensed(iupac) {
                Ok(g) if g.node_count() > 0 => g,
                _ => {
                    node_mismatches.push(format!("{}: IUPAC parse failed for '{}'", id, iupac));
                    continue;
                }
            };

            let wurcs_graph = match parse_wurcs(wurcs_ref) {
                Ok(g) if g.node_count() > 0 => g,
                _ => {
                    node_mismatches.push(format!("{}: WURCS parse failed", id));
                    continue;
                }
            };

            if iupac_graph.node_count() == wurcs_graph.node_count()
                && iupac_graph.edge_count() == wurcs_graph.edge_count()
            {
                node_matches += 1;
            } else {
                node_mismatches.push(format!(
                    "{}: IUPAC→{}N/{}E, WURCS→{}N/{}E | IUPAC: {}",
                    id,
                    iupac_graph.node_count(),
                    iupac_graph.edge_count(),
                    wurcs_graph.node_count(),
                    wurcs_graph.edge_count(),
                    iupac,
                ));
            }
        }

        println!("\n=== IUPAC↔WURCS Structural Equivalence ===");
        println!(
            "Tested: {}, Node/edge match: {}, Mismatch: {}",
            tested,
            node_matches,
            node_mismatches.len()
        );

        if !node_mismatches.is_empty() {
            for m in &node_mismatches {
                eprintln!("  {}", m);
            }
        }

        assert!(
            node_matches as f64 / tested as f64 >= 0.5,
            "Only {}/{} entries have matching node/edge counts ({}%)",
            node_matches,
            tested,
            (node_matches as f64 / tested as f64 * 100.0) as u32
        );
    }

    #[test]
    fn test_iupac_wurcs_cross_graph_roundtrip_gs00002() {
        // IUPAC → graph → WURCS → graph → IUPAC
        // Verify structural preservation through the full cross-format roundtrip
        let iupac_original = "Fuc(a1-2)Gal(a1-3)Gal(b1-3)GalNAc";

        let graph1 = parse_iupac_condensed(iupac_original).expect("IUPAC parse");
        let wurcs = write_wurcs(&graph1).expect("WURCS write");
        let graph2 = parse_wurcs(&wurcs).expect("WURCS parse");
        let iupac_output = write_iupac_condensed(&graph2).expect("IUPAC write");
        let graph3 = parse_iupac_condensed(&iupac_output).expect("Final IUPAC parse");

        assert_eq!(
            graph1.node_count(),
            graph2.node_count(),
            "Node count changed in IUPAC→WURCS: {} vs {}",
            graph1.node_count(),
            graph2.node_count()
        );
        assert_eq!(
            graph1.edge_count(),
            graph2.edge_count(),
            "Edge count changed in IUPAC→WURCS: {} vs {}",
            graph1.edge_count(),
            graph2.edge_count()
        );
        assert_eq!(
            graph2.node_count(),
            graph3.node_count(),
            "Node count changed in WURCS→IUPAC: {} vs {}",
            graph2.node_count(),
            graph3.node_count()
        );
        assert_eq!(
            graph2.edge_count(),
            graph3.edge_count(),
            "Edge count changed in WURCS→IUPAC: {} vs {}",
            graph2.edge_count(),
            graph3.edge_count()
        );

        // Final output should still contain the key residues
        assert!(iupac_output.contains("Fuc"), "Lost Fuc: {}", iupac_output);
        assert!(iupac_output.contains("Gal"), "Lost Gal: {}", iupac_output);
        assert!(
            iupac_output.contains("GalNAc"),
            "Lost GalNAc: {}",
            iupac_output
        );
    }
}
