use crabwurcs_core::{parse_wurcs, write_wurcs};
use crabwurcs_iupac::{
    parse_glycam, parse_iupac_condensed, parse_iupac_extended, write_glycam, write_iupac_condensed,
    write_iupac_extended,
};

#[test]
fn condensed_composition_is_editable_and_reaches_wurcs() {
    let input = "{GlcNAc}2,{Man}3,{Fuc}1";
    let mut graph = parse_iupac_condensed(input).unwrap();
    assert!(graph.is_composition());
    assert_eq!(graph.node_count(), 6);
    assert_eq!(write_iupac_condensed(&graph).unwrap(), input);

    let _ = graph.inner_mut();
    assert_eq!(write_iupac_condensed(&graph).unwrap(), input);
    let wurcs = write_wurcs(&graph).unwrap();
    assert!(wurcs.contains("/3,6,0+/"), "{wurcs}");
    let reparsed = parse_wurcs(&wurcs).unwrap();
    assert!(reparsed.is_composition());
    assert_eq!(reparsed.node_count(), 6);
}

#[test]
fn extended_and_glycam_compositions_roundtrip() {
    let extended = "{D-GlcpNAc}2,{D-Manp}3,{L-Fucp}1";
    let mut graph = parse_iupac_extended(extended).unwrap();
    assert!(graph.is_composition());
    assert_eq!(graph.node_count(), 6);
    assert_eq!(write_iupac_extended(&graph).unwrap(), extended);
    let _ = graph.inner_mut();
    let regenerated = write_iupac_extended(&graph).unwrap();
    assert!(regenerated.contains("{D-GlcpNAc}2"), "{regenerated}");
    assert!(regenerated.contains("{D-Manp}3"), "{regenerated}");
    assert!(regenerated.contains("{L-Fucp}1"), "{regenerated}");

    let glycam = "{DGlcpNAc}2,{DManp}3,{LFucp}1";
    let mut graph = parse_glycam(glycam).unwrap();
    assert_eq!(write_glycam(&graph).unwrap(), glycam);
    let _ = graph.inner_mut();
    let regenerated = write_glycam(&graph).unwrap();
    assert!(regenerated.contains("{DGlcpNAc}2"), "{regenerated}");
    assert!(regenerated.contains("{DManp}3"), "{regenerated}");
    assert!(regenerated.contains("{LFucp}1"), "{regenerated}");
}

#[test]
fn wurcs_repeat_and_cycle_export_closure_annotations() {
    let repeat = "WURCS=2.0/1,3,3/[a2122h-1a_1-5]/1-1-1/a4-b1_b4-c1_a1-c4~2-5";
    let mut graph = parse_wurcs(repeat).unwrap();
    let _ = graph.inner_mut();
    let iupac = write_iupac_condensed(&graph).unwrap();
    assert!(iupac.starts_with("[4)"), "{iupac}");
    assert!(iupac.ends_with("(a1-]2-5"), "{iupac}");
    let reparsed = parse_iupac_condensed(&iupac).unwrap();
    assert!(reparsed.inner().edge_weights().any(|linkage| {
        linkage.repeat
            == Some(crabwurcs_core::RepeatCount::Range {
                min: Some(2),
                max: Some(5),
            })
    }));
    let extended = write_iupac_extended(&graph).unwrap();
    assert!(extended.starts_with("[4)-"), "{extended}");
    assert!(extended.ends_with("-(1→]2-5"), "{extended}");
    let reparsed_extended = parse_iupac_extended(&extended).unwrap();
    assert!(reparsed_extended
        .inner()
        .edge_weights()
        .any(|linkage| linkage.repeat.is_some()));

    let cyclic = "WURCS=2.0/1,3,3/[a2122h-1a_1-5]/1-1-1/a1-c4_a4-b1_b4-c1";
    let mut graph = parse_wurcs(cyclic).unwrap();
    let _ = graph.inner_mut();
    let iupac = write_iupac_condensed(&graph).unwrap();
    assert!(iupac.starts_with("4)"), "{iupac}");
    assert!(iupac.ends_with("(a1-"), "{iupac}");
    let reparsed = parse_iupac_condensed(&iupac).unwrap();
    assert_eq!(
        reparsed
            .inner()
            .edge_weights()
            .filter(|linkage| linkage.cyclic)
            .count(),
        1
    );
}

#[test]
fn undefined_fragment_exports_reference_style_anchor_labels() {
    let wurcs = "WURCS=2.0/6,11,10/[a2122h-1x_1-5_2*NCC/3=O][a2122h-1b_1-5_2*NCC/3=O][a1122h-1b_1-5][a1122h-1a_1-5][a2112h-1b_1-5][a1221m-1a_1-5]/1-2-3-4-2-5-4-2-6-2-5/a4-b1_a6-i1_b4-c1_c3-d1_c6-g1_d2-e1_e4-f1_g2-h1_j4-k1_j1-d4|d6|g4|g6}";
    let mut graph = parse_wurcs(wurcs).unwrap();
    assert_eq!(graph.edge_count(), 9);
    assert_eq!(graph.undefined_linkages().len(), 1);
    let _ = graph.inner_mut();

    let iupac = write_iupac_condensed(&graph).unwrap();
    assert!(iupac.contains("=1$,"), "{iupac}");
    assert!(iupac.matches("1$").count() >= 3, "{iupac}");
    assert!(iupac.contains("-4/6)"), "{iupac}");

    let reparsed = parse_iupac_condensed(&iupac).unwrap();
    assert_eq!(reparsed.node_count(), 11);
    assert_eq!(reparsed.edge_count(), 9);
    assert_eq!(reparsed.undefined_linkages().len(), 1);
    assert_eq!(reparsed.undefined_linkages()[0].parents.len(), 2);

    let extended = write_iupac_extended(&graph).unwrap();
    assert!(extended.contains("=1$,"), "{extended}");
    assert!(extended.contains('→'), "{extended}");
    let reparsed = parse_iupac_extended(&extended).unwrap();
    assert_eq!(reparsed.node_count(), 11);
    assert_eq!(reparsed.edge_count(), 9);
    assert_eq!(reparsed.undefined_linkages().len(), 1);
}

#[test]
fn linkage_probability_roundtrips_between_wurcs_and_iupac() {
    let wurcs = "WURCS=2.0/2,2,1/[a2122h-1x_1-5][a2112h-1b_1-5]/1-2/b1-a4%.55%";
    let mut graph = parse_wurcs(wurcs).unwrap();
    let _ = graph.inner_mut();
    let condensed = write_iupac_condensed(&graph).unwrap();
    assert!(condensed.contains("-55%4)"), "{condensed}");
    let reparsed = parse_iupac_condensed(&condensed).unwrap();
    assert_eq!(
        reparsed
            .inner()
            .edge_weights()
            .next()
            .unwrap()
            .parent_probability,
        graph
            .inner()
            .edge_weights()
            .next()
            .unwrap()
            .parent_probability
    );

    let extended = write_iupac_extended(&graph).unwrap();
    assert!(extended.contains("→55%4"), "{extended}");
    let reparsed = parse_iupac_extended(&extended).unwrap();
    assert!(reparsed
        .inner()
        .edge_weights()
        .next()
        .unwrap()
        .parent_probability
        .is_some());
}

#[test]
fn substituent_probability_uses_reference_iupac_annotation() {
    let unknown = "WURCS=2.0/1,1,0/[a2122A-1a_1-5_6%?%*OC]/1/";
    let mut graph = parse_wurcs(unknown).unwrap();
    let _ = graph.inner_mut();
    let condensed = write_iupac_condensed(&graph).unwrap();
    assert_eq!(condensed, "GlcA6(?%)Me");
    let reparsed = parse_iupac_condensed(&condensed).unwrap();
    let modification = &reparsed
        .residue(reparsed.root().unwrap())
        .unwrap()
        .modifications[0];
    assert_eq!(
        modification.probability,
        Some(crabwurcs_core::Probability {
            lower: crabwurcs_core::ProbabilityValue::Unknown,
            upper: crabwurcs_core::ProbabilityValue::Unknown,
        })
    );

    let exact = "WURCS=2.0/1,1,0/[a2211m-1b_1-5_2%.5%*OCC/3=O]/1/";
    let mut graph = parse_wurcs(exact).unwrap();
    let _ = graph.inner_mut();
    let condensed = write_iupac_condensed(&graph).unwrap();
    assert_eq!(condensed, "Rha2(50%)Ac");
    let reparsed = parse_iupac_condensed(&condensed).unwrap();
    assert_eq!(
        reparsed
            .residue(reparsed.root().unwrap())
            .unwrap()
            .modifications[0]
            .probability,
        Some(crabwurcs_core::Probability {
            lower: crabwurcs_core::ProbabilityValue::Known(5000),
            upper: crabwurcs_core::ProbabilityValue::Known(5000),
        })
    );
}

#[test]
fn map_bridge_roundtrips_through_condensed_iupac() {
    let wurcs = "WURCS=2.0/2,2,1/[hxh][a2122h-1b_1-5]/1-2/a3n2-b1n1*1NCCOP^XO*2/6O/6=O";
    let mut graph = parse_wurcs(wurcs).unwrap();
    let _ = graph.inner_mut();
    let condensed = write_iupac_condensed(&graph).unwrap();
    assert!(condensed.contains("(b1-1PEtn2-3)"), "{condensed}");

    let parsed = parse_iupac_condensed(&condensed).unwrap();
    let bridge = parsed
        .inner()
        .edge_weights()
        .find(|linkage| linkage.map_code.is_some())
        .unwrap();
    assert_eq!(bridge.map_code.as_deref(), Some("*1NCCOP^XO*2/6O/6=O"));
    assert_eq!(bridge.child_modification_position, Some(1));
    assert_eq!(bridge.parent_modification_position, Some(2));

    let generated = write_wurcs(&parsed).unwrap();
    assert!(
        generated.contains("n2-") && generated.contains("n1*1NCCOP^XO*2/6O/6=O"),
        "{generated}"
    );

    let extended = write_iupac_extended(&graph).unwrap();
    assert!(extended.contains("(1-1PEtn2→3)"), "{extended}");
    let reparsed_extended = parse_iupac_extended(&extended).unwrap();
    assert!(reparsed_extended
        .inner()
        .edge_weights()
        .any(|linkage| { linkage.map_code.as_deref() == Some("*1NCCOP^XO*2/6O/6=O") }));

    let glycam = write_glycam(&graph).unwrap();
    assert!(glycam.contains("b1-1PEtn2-3"), "{glycam}");
    let reparsed_glycam = parse_glycam(&glycam).unwrap();
    assert!(reparsed_glycam
        .inner()
        .edge_weights()
        .any(|linkage| { linkage.map_code.as_deref() == Some("*1NCCOP^XO*2/6O/6=O") }));

    let phosphate = "WURCS=2.0/2,2,1/[a2122h-1x_1-5][a2112h-1b_1-5]/1-2/a3-b1*OPO*/3O/3=O";
    let mut graph = parse_wurcs(phosphate).unwrap();
    let _ = graph.inner_mut();
    let condensed = write_iupac_condensed(&graph).unwrap();
    assert!(condensed.contains("(b1-P-3)"), "{condensed}");
}
