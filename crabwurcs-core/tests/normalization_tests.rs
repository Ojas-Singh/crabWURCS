use crabwurcs_core::{normalize_wurcs, parse_wurcs, write_wurcs_canonical};

#[test]
fn glycanformatconverter_golden_cases_match_and_are_idempotent() {
    for line in include_str!("fixtures/wurcs_normalization.tsv").lines() {
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let (input, expected) = line.split_once('\t').expect("fixture must be TSV");
        let actual = normalize_wurcs(input).expect("fixture must normalize");
        assert_eq!(actual, expected, "input: {input}");
        assert_eq!(normalize_wurcs(&actual).unwrap(), actual);

        let reparsed = parse_wurcs(&actual).unwrap();
        assert_eq!(write_wurcs_canonical(&reparsed).unwrap(), actual);
    }
}

#[test]
fn modeled_advanced_constructs_normalize_idempotently() {
    let cases = [
        // Composition.
        "WURCS=2.0/4,15,0+/[AUd21122h_5*NCC/3=O][uxxxxh_2*NCC/3=O][uxxxxh][u1221m]/1-2-2-2-2-2-2-3-3-3-4-4-4-4-4/",
        // Cycle.
        "WURCS=2.0/1,7,7/[a2122h-1a_1-5_2*OC_3*OC_6*N]/1-1-1-1-1-1-1/a1-g4_a4-b1_b4-c1_c4-d1_d4-e1_e4-f1_f4-g1",
        // Variable repeat.
        "WURCS=2.0/4,5,5/[a2112h-1b_1-5_2*NCC/3=O][a2122A-1b_1-5][a2112h-1a_1-5][a2122h-1b_1-5]/1-2-1-3-4/a3-b1_b4-c1_c4-d1_d3-e1_a1-e4~n",
        // Undefined substituent attachment.
        "WURCS=2.0/2,2,1/[u2122h][u2112h]/1-2/a?|b?}*OCC/3=O",
        // Alternative position and endpoint probability.
        "WURCS=2.0/1,2,1/[a2122h-1b_1-5]/1-1/b1-a4%.55%",
        // Indexed MAP bridge.
        "WURCS=2.0/2,2,1/[hxh][a2122h-1b_1-5]/1-2/a3n2-b1n1*1NCCOP^XO*2/6O/6=O",
        // Open-chain backbone.
        "WURCS=2.0/1,1,0/[o261m_3*CO]/1/",
    ];

    for input in cases {
        let before = parse_wurcs(input).unwrap();
        let normalized = normalize_wurcs(input).unwrap();
        assert_eq!(normalize_wurcs(&normalized).unwrap(), normalized);
        let after = parse_wurcs(&normalized).unwrap();
        assert_eq!(before.node_count(), after.node_count(), "{input}");
        assert_eq!(before.edge_count(), after.edge_count(), "{input}");
    }
}
