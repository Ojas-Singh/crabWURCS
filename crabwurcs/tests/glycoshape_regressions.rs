//! Regression tests for defects found while auditing the GlycoShape database
//! (entries named by GlycoShape ID). They use only the public API so they
//! can be run against earlier releases.

use crabwurcs::snfg::{Shape, colour, symbol_for};
use crabwurcs::{Format, ResidueGraph, convert, normalize_wurcs, parse_notation, write_notation};

fn wurcs_of(notation: &str, format: Format) -> String {
    let graph = parse_notation(notation, format)
        .unwrap_or_else(|error| panic!("cannot read {notation}: {error}"));
    crabwurcs::core::write_wurcs(&graph).unwrap()
}

fn root_residue(graph: &ResidueGraph) -> &crabwurcs::core::Monosaccharide {
    graph.residue(graph.root().unwrap()).unwrap()
}

/// GLYCAM substituent lists (`[9A]`, `[2A]`, `[6PC]`, `[2Ac,4Ac]`) are read.
///
/// Expected WURCS: GlyTouCan (GlyCosmos `sparqlist/gtcid2seqs`, looked up
/// 2026-10-14) for the archetype accession of each entry; GlyGen holds the
/// same WURCS for G92493ZW and G35061DQ and has no record of the others.
#[test]
fn glycam_substituent_lists_are_read() {
    for (entry, name, accession, glytoucan) in [
        (
            "GS00237",
            "DNeup5Ac[9A]a2-6DGalpb1-4DGlcpNAc",
            "G92493ZW",
            "WURCS=2.0/3,3,2/[u2122h_2*NCC/3=O][a2112h-1b_1-5][Aad21122h-2a_2-6_5*NCC/3=O_9*OCC/3=O]/1-2-3/a4-b1_b6-c2",
        ),
        (
            "GS00238",
            "DNeup5Ac[9A]a2-3DGalpb1-4DGlcpNAc",
            "G35061DQ",
            "WURCS=2.0/3,3,2/[u2122h_2*NCC/3=O][a2112h-1b_1-5][Aad21122h-2a_2-6_5*NCC/3=O_9*OCC/3=O]/1-2-3/a4-b1_b3-c2",
        ),
        (
            "GS00314",
            "DNeup5Ac[9A]a2-3DGalpb1-4[LFucpa1-3]DGlcpNAc",
            "G65204RX",
            "WURCS=2.0/4,4,3/[u2122h_2*NCC/3=O][a1221m-1a_1-5][a2112h-1b_1-5][Aad21122h-2a_2-6_5*NCC/3=O_9*OCC/3=O]/1-2-3-4/a3-b1_a4-c1_c3-d2",
        ),
        (
            "GS00386",
            "DGalpNAcb1-4[LFucpa1-3]DGlcpNAc[6PC]b1-2DManpa1-3[DManpa1-6]DManpb1-4DGlcpNAcb1-4[LFucpa1-6]DGlcpNAc",
            "G20461UY",
            "WURCS=2.0/7,9,8/[u2122h_2*NCC/3=O][a2122h-1b_1-5_2*NCC/3=O][a1122h-1b_1-5][a1122h-1a_1-5][a2122h-1b_1-5_2*NCC/3=O_6*OP^XOCCNC/7C/7C/3O/3=O][a1221m-1a_1-5][a2112h-1b_1-5_2*NCC/3=O]/1-2-3-4-5-6-7-4-6/a4-b1_a6-i1_b4-c1_c3-d1_c6-h1_d2-e1_e3-f1_e4-g1",
        ),
        (
            "GS00667",
            "DGlcpNAcb1-2[DGlcpa1-3]LRhapa1-2LRhapa1-3LRhap[2A]",
            "G61219EA",
            "WURCS=2.0/4,5,4/[u2211m_2*OCC/3=O][a2211m-1a_1-5][a2122h-1b_1-5_2*NCC/3=O][a2122h-1a_1-5]/1-2-2-3-4/a3-b1_b2-c1_c2-d1_c3-e1",
        ),
        (
            // The archetype name leaves the reducing end open; the beta name
            // is the structure of the accession.
            "GS00954",
            "DGlcpb1-3[DGalpNAca1-4DGalpNAca1-4]DGalpNAca1-4DGalpNAca1-4DGalpNAca1-3DBacp[2Ac,4Ac]b1-OH",
            "G66397DB",
            "WURCS=2.0/3,7,6/[a2122m-1b_1-5_2*NCC/3=O_4*NCC/3=O][a2112h-1a_1-5_2*NCC/3=O][a2122h-1b_1-5]/1-2-2-2-3-2-2/a3-b1_b4-c1_c4-d1_d3-e1_d4-f1_f4-g1",
        ),
    ] {
        assert_eq!(
            wurcs_of(name, Format::Glycam),
            normalize_wurcs(glytoucan).unwrap(),
            "{entry} ({accession}): {name}"
        );
        // The alpha name (every reducing end here is an aldose) keeps the
        // substituent too; the stored alpha/beta WURCS of these entries lack it.
        let alpha = format!("{}a1-OH", name.trim_end_matches("b1-OH"));
        let wurcs = wurcs_of(&alpha, Format::Glycam);
        for substituent in ["9*OCC/3=O", "6*OP^XOCCNC", "2*OCC/3=O_", "4*NCC/3=O"] {
            if glytoucan.contains(substituent) {
                assert!(wurcs.contains(substituent), "{entry}: {alpha} -> {wurcs}");
            }
        }
    }
    // An unknown GLYCAM code is an error, not a silently dropped substituent.
    assert!(parse_notation("DGlcpNAc[6X]b1-OH", Format::Glycam).is_err());
}

/// Residue names carry the configuration and the ring size where they
/// differ from the default of the name.
#[test]
fn residue_names_carry_configuration_and_ring_size() {
    for (wurcs, expected) in [
        ("WURCS=2.0/1,1,0/[u2112m]/1/", "D-Fuc"),
        ("WURCS=2.0/1,1,0/[u1221m]/1/", "Fuc"),
        ("WURCS=2.0/1,1,0/[u122h]/1/", "D-Ara"),
        ("WURCS=2.0/1,1,0/[u2211h]/1/", "L-Man"),
        ("WURCS=2.0/1,1,0/[u1122m]/1/", "D-Rha"),
        ("WURCS=2.0/1,1,0/[u1121h]/1/", "L-Gul"),
        ("WURCS=2.0/1,1,0/[u2221h]/1/", "L-Tal"),
        (
            "WURCS=2.0/2,2,1/[u2122h][a211h-1b_1-4]/1-2/a4-b1",
            "Araf(b1-4)Glc",
        ),
        (
            "WURCS=2.0/2,2,1/[u2122h][a211h-1b_1-5]/1-2/a4-b1",
            "Ara(b1-4)Glc",
        ),
        (
            "WURCS=2.0/2,2,1/[u2122h][a122h-1b_1-4]/1-2/a4-b1",
            "D-Araf(b1-4)Glc",
        ),
        (
            "WURCS=2.0/2,2,1/[u2122h][a2112h-1b_1-4]/1-2/a4-b1",
            "Galf(b1-4)Glc",
        ),
        (
            "WURCS=2.0/2,2,1/[u2122h][ha122h-2b_2-5]/1-2/a1-b2",
            "Fruf(b2-1)Glc",
        ),
        (
            "WURCS=2.0/2,2,1/[u2122h][a2112m-1a_1-5]/1-2/a2-b1",
            "D-Fuc(a1-2)Glc",
        ),
    ] {
        assert_eq!(
            convert(wurcs, Format::Wurcs, Format::IupacCondensed).unwrap(),
            expected,
            "{wurcs}"
        );
    }
}

/// SNFG: D-Fuc and D-Rha are drawn with the Fuc and Rha triangles and a
/// configuration letter, not with the Gal and Man circles (GS00025, GS00657).
#[test]
fn snfg_draws_rare_configurations_with_their_own_symbols() {
    for (wurcs, shape, fill) in [
        (
            "WURCS=2.0/1,1,0/[a2112m-1a_1-5]/1/",
            Shape::Triangle,
            colour::RED,
        ),
        (
            "WURCS=2.0/1,1,0/[a1122m-1a_1-5]/1/",
            Shape::Triangle,
            colour::GREEN,
        ),
        (
            "WURCS=2.0/1,1,0/[a1211m-1a_1-5]/1/",
            Shape::Triangle,
            colour::BLUE,
        ),
    ] {
        let graph = parse_notation(wurcs, Format::Wurcs).unwrap();
        let symbol = symbol_for(root_residue(&graph)).unwrap();
        assert_eq!((symbol.shape, symbol.fill), (shape, fill), "{wurcs}");
        let svg = crabwurcs::render_svg(&graph).unwrap();
        assert!(svg.contains("class=\"config-label\">"), "{wurcs}");
    }
}

/// The 28 GlycoShape glycans whose written names did not return unchanged
/// (21 with a furanose, 7 with rare residues): the WURCS of each name, and of
/// its alpha and beta forms, returns from GLYCAM, IUPAC condensed and IUPAC
/// extended. The reducing-end anomer is written too.
#[test]
fn glycoshape_names_return_unchanged_through_every_written_notation() {
    const NAMES: &[(&str, &str, &str, &str)] = &[
        (
            "GS00025",
            "DFucpa1-2DGalpb1-3DGlcpNAcb1-3DGalpb1-4DGlcp",
            "DFucpa1-2DGalpb1-3DGlcpNAcb1-3DGalpb1-4DGlcpa1-OH",
            "DFucpa1-2DGalpb1-3DGlcpNAcb1-3DGalpb1-4DGlcpb1-OH",
        ),
        (
            "GS00430",
            "DGalpb1-4DGalpa1-3[2,4-diacetimido-2,4,6-trideoxyhexose]",
            "DGalpb1-4DGalpa1-3[2,4-diacetimido-2,4,6-trideoxyhexose]aNone-OH",
            "DGalpb1-4DGalpa1-3[2,4-diacetimido-2,4,6-trideoxyhexose]bNone-OH",
        ),
        (
            "GS00472",
            "LArafb1-2LArafb1-2LAraf",
            "LArafb1-2LArafb1-2LArafa1-OH",
            "LArafb1-2LArafb1-2LArafb1-OH",
        ),
        (
            "GS00473",
            "LArafb1-2LAraf",
            "LArafb1-2LArafa1-OH",
            "LArafb1-2LArafb1-OH",
        ),
        (
            "GS00475",
            "LArafa1-3LArafb1-2LArafb1-2LAraf",
            "LArafa1-3LArafb1-2LArafb1-2LArafa1-OH",
            "LArafa1-3LArafb1-2LArafb1-2LArafb1-OH",
        ),
        (
            "GS00480",
            "LArafa1-2LArafa1-3LArafb1-2LArafb1-2LAraf",
            "LArafa1-2LArafa1-3LArafb1-2LArafb1-2LArafa1-OH",
            "LArafa1-2LArafa1-3LArafb1-2LArafb1-2LArafb1-OH",
        ),
        (
            "GS00482",
            "LArafa1-3LArafb1-2LArafb1-2[LArafa1-3]LAraf",
            "LArafa1-3LArafb1-2LArafb1-2[LArafa1-3]LArafa1-OH",
            "LArafa1-3LArafb1-2LArafb1-2[LArafa1-3]LArafb1-OH",
        ),
        (
            "GS00483",
            "LArafb1-3LArafb1-2LArafb1-2LAraf",
            "LArafb1-3LArafb1-2LArafb1-2LArafa1-OH",
            "LArafb1-3LArafb1-2LArafb1-2LArafb1-OH",
        ),
        (
            "GS00502",
            "DFrufb2-1DFrufb2-1DFruf",
            "DFrufb2-1DFrufb2-1DFrufa2-OH",
            "DFrufb2-1DFrufb2-1DFrufb2-OH",
        ),
        (
            "GS00534",
            "DGalfb1-2LArafb1-2LAraf",
            "DGalfb1-2LArafb1-2LArafa1-OH",
            "DGalfb1-2LArafb1-2LArafb1-OH",
        ),
        (
            "GS00556",
            "LArafa1-3DGalfb1-2LArafb1-2LAraf",
            "LArafa1-3DGalfb1-2LArafb1-2LArafa1-OH",
            "LArafa1-3DGalfb1-2LArafb1-2LArafb1-OH",
        ),
        (
            "GS00557",
            "LArafb1-2LArafb1-2[LArafa1-3]LAraf",
            "LArafb1-2LArafb1-2[LArafa1-3]LArafa1-OH",
            "LArafb1-2LArafb1-2[LArafa1-3]LArafb1-OH",
        ),
        (
            "GS00558",
            "LArafb1-2LArafb1-2LArafb1-2LAraf",
            "LArafb1-2LArafb1-2LArafb1-2LArafa1-OH",
            "LArafb1-2LArafb1-2LArafb1-2LArafb1-OH",
        ),
        (
            "GS00562",
            "DArafb1-2[DArafa1-3]DAraf",
            "DArafb1-2[DArafa1-3]DArafa1-OH",
            "DArafb1-2[DArafa1-3]DArafb1-OH",
        ),
        (
            "GS00569",
            "LArafb1-2LArafb1-3LArafb1-2LArafb1-2LAraf",
            "LArafb1-2LArafb1-3LArafb1-2LArafb1-2LArafa1-OH",
            "LArafb1-2LArafb1-3LArafb1-2LArafb1-2LArafb1-OH",
        ),
        (
            "GS00572",
            "LArafa1-2LArafa1-3DGalfb1-2LArafb1-2LAraf",
            "LArafa1-2LArafa1-3DGalfb1-2LArafb1-2LArafa1-OH",
            "LArafa1-2LArafa1-3DGalfb1-2LArafb1-2LArafb1-OH",
        ),
        (
            "GS00616",
            "DGalfb1-2LArafb1-2[LArafa1-3]LAraf",
            "DGalfb1-2LArafb1-2[LArafa1-3]LArafa1-OH",
            "DGalfb1-2LArafb1-2[LArafa1-3]LArafb1-OH",
        ),
        (
            "GS00617",
            "DGalfb1-2[LArafa1-5]LArafb1-2LAraf",
            "DGalfb1-2[LArafa1-5]LArafb1-2LArafa1-OH",
            "DGalfb1-2[LArafa1-5]LArafb1-2LArafb1-OH",
        ),
        (
            "GS00621",
            "LArafa1-3DGalfb1-2LArafb1-2[LArafa1-3]LAraf",
            "LArafa1-3DGalfb1-2LArafb1-2[LArafa1-3]LArafa1-OH",
            "LArafa1-3DGalfb1-2LArafb1-2[LArafa1-3]LArafb1-OH",
        ),
        (
            "GS00632",
            "DGalfb1-2[LArafa1-5]LArafb1-2[LArafa1-3]LAraf",
            "DGalfb1-2[LArafa1-5]LArafb1-2[LArafa1-3]LArafa1-OH",
            "DGalfb1-2[LArafa1-5]LArafb1-2[LArafa1-3]LArafb1-OH",
        ),
        (
            "GS00644",
            "DFrufb2-4DFruf",
            "DFrufb2-4DFrufa2-OH",
            "DFrufb2-4DFrufb2-OH",
        ),
        (
            "GS00650",
            "LManpa1-2LGulp",
            "LManpa1-2LGulpa1-OH",
            "LManpa1-2LGulpb1-OH",
        ),
        (
            "GS00657",
            "DGlcpNAcb1-3DRhapa1-3DGlcpa1-6DGlcpNAcb1-2DGalp",
            "DGlcpNAcb1-3DRhapa1-3DGlcpa1-6DGlcpNAcb1-2DGalpa1-OH",
            "DGlcpNAcb1-3DRhapa1-3DGlcpa1-6DGlcpNAcb1-2DGalpb1-OH",
        ),
        (
            "GS00664",
            "DGalpa1-2[DManpa1-3DRhapa1-3][LRhap[2Me,3Me]a1-2[DArapb1-3]LRhapb1-4DXylpb1-4]LFucpa1-3[DXylpb1-4]DGlcp",
            "DGalpa1-2[DManpa1-3DRhapa1-3][LRhap[2Me,3Me]a1-2[DArapb1-3]LRhapb1-4DXylpb1-4]LFucpa1-3[DXylpb1-4]DGlcpa1-OH",
            "DGalpa1-2[DManpa1-3DRhapa1-3][LRhap[2Me,3Me]a1-2[DArapb1-3]LRhapb1-4DXylpb1-4]LFucpa1-3[DXylpb1-4]DGlcpb1-OH",
        ),
        (
            "GS00675",
            "DGalpa1-2[DManpa1-3DRhapa1-3][LRhap[2Me,3Me]a1-2[LArafb1-3]LRhapb1-4DXylpb1-4]LFucpa1-3[DXylpb1-4]DGlcp",
            "DGalpa1-2[DManpa1-3DRhapa1-3][LRhap[2Me,3Me]a1-2[LArafb1-3]LRhapb1-4DXylpb1-4]LFucpa1-3[DXylpb1-4]DGlcpa1-OH",
            "DGalpa1-2[DManpa1-3DRhapa1-3][LRhap[2Me,3Me]a1-2[LArafb1-3]LRhapb1-4DXylpb1-4]LFucpa1-3[DXylpb1-4]DGlcpb1-OH",
        ),
        (
            "GS00890",
            "DArafb1-2DArafa1-5DArafa1-5DAraf",
            "DArafb1-2DArafa1-5DArafa1-5DArafa1-OH",
            "DArafb1-2DArafa1-5DArafa1-5DArafb1-OH",
        ),
        ("GS00908", "LGulp", "LGulpa1-OH", "LGulpb1-OH"),
        ("GS00912", "LTalp", "LTalpa1-OH", "LTalpb1-OH"),
    ];
    assert_eq!(NAMES.len(), 28);
    for (entry, archetype, alpha, beta) in NAMES {
        for name in [archetype, alpha, beta] {
            let wurcs = wurcs_of(name, Format::Glycam);
            let graph = parse_notation(&wurcs, Format::Wurcs).unwrap();
            for format in [
                Format::Glycam,
                Format::IupacCondensed,
                Format::IupacExtended,
            ] {
                let written = write_notation(&graph, format).unwrap_or_else(|error| {
                    panic!("{entry}: cannot write {wurcs} as {format:?}: {error}")
                });
                assert_eq!(
                    wurcs_of(&written, format),
                    wurcs,
                    "{entry} {name}: {format:?} {written}"
                );
            }
        }
    }
    // The reducing-end anomer is part of every written name.
    let alpha = wurcs_of("LArafb1-2LArafa1-OH", Format::Glycam);
    assert_eq!(
        convert(&alpha, Format::Wurcs, Format::Glycam).unwrap(),
        "LArafb1-2LArafa1-OH"
    );
    assert_eq!(
        convert(&alpha, Format::Wurcs, Format::IupacCondensed).unwrap(),
        "Araf(b1-2)Araf(a1-"
    );
    assert_eq!(
        convert(&alpha, Format::Wurcs, Format::IupacExtended).unwrap(),
        "β-L-Araf-(1→2)-α-L-Araf-(1→"
    );
}

/// A furanose ring closes at O4 of an aldopentose or aldohexose (WURCS
/// `1-4`) and at O5 of a 2-ketohexose (`2-5`), from every notation; the
/// pyranose closures `1-5` / `2-6` stored for 25 GlycoShape furanose entries
/// are a different ring.
#[test]
fn furanose_rings_close_at_the_right_oxygen() {
    for (notation, format, ring) in [
        ("LArafa1-OH", Format::Glycam, "[a211h-1a_1-4]"),
        ("Araf(a1-", Format::IupacCondensed, "[a211h-1a_1-4]"),
        ("α-L-Araf-(1→", Format::IupacExtended, "[a211h-1a_1-4]"),
        ("DGalfb1-OH", Format::Glycam, "[a2112h-1b_1-4]"),
        ("DRibfa1-OH", Format::Glycam, "[a222h-1a_1-4]"),
        ("DFrufb2-OH", Format::Glycam, "[ha122h-2b_2-5]"),
        ("Fruf(b2-", Format::IupacCondensed, "[ha122h-2b_2-5]"),
        ("DPsifa2-OH", Format::Glycam, "[ha222h-2a_2-5]"),
    ] {
        let wurcs = wurcs_of(notation, format);
        assert!(wurcs.contains(ring), "{notation} -> {wurcs}");
        // and the ring survives a trip through IUPAC condensed
        let iupac = convert(&wurcs, Format::Wurcs, Format::IupacCondensed).unwrap();
        assert_eq!(wurcs_of(&iupac, Format::IupacCondensed), wurcs, "{iupac}");
    }
}

/// GlycoShape's multi-model PDB files open with free-text REMARK lines and
/// repeat every atom in each model; the first model is read.
#[test]
fn pdb_files_with_free_text_remarks_and_several_models_are_read() {
    let atom = |serial: u32, name: &str, residue: &str, sequence: u32, x: f64, y: f64, z: f64| {
        let element = &name[..1];
        format!(
            "ATOM  {serial:5} {name:<4} {residue:3} X{sequence:4}    {x:8.3}{y:8.3}{z:8.3}  1.00  0.00          {element:>2}"
        )
    };
    let model = [
        atom(1, "C1", "NAG", 1, 0.0, 0.0, 0.0),
        atom(2, "C2", "NAG", 1, 1.5, 0.0, 0.0),
        atom(3, "C3", "NAG", 1, 2.0, 1.4, 0.0),
        atom(4, "C4", "NAG", 1, 1.3, 2.6, 0.0),
        atom(5, "O4", "NAG", 1, 1.8, 3.9, 0.0),
        atom(6, "C5", "NAG", 1, -0.2, 2.5, 0.0),
        atom(7, "O5", "NAG", 1, -0.6, 1.2, 0.0),
        atom(8, "C1", "GAL", 2, 1.4, 5.3, 0.0),
        atom(9, "C2", "GAL", 2, 2.6, 6.2, 0.0),
        atom(10, "O5", "GAL", 2, 0.3, 5.9, 0.0),
    ]
    .join("\n");
    let bonds = "CONECT    5    4    8\nCONECT    8    5    9   10\n";
    let pdb = format!(
        "REMARK     GENERATED BY GlycanAnalysisPipeline from GlycoShape\n\
         REMARK        ____ _                ____  _\n\
         REMARK    IUPAC Name: Gal(b1-4)GlcNAc(b1-\n\
         REMARK   GLYCAN MULTI-MODEL PDB  LEVEL LEVEL_1  ANOMER BETA\n\
         MODEL        1\n{model}\n{bonds}ENDMDL\nMODEL        2\n{model}\n{bonds}ENDMDL\nEND\n"
    );
    let glycans =
        crabwurcs::extract_glycans_from_str(&pdb, false).unwrap_or_else(|error| panic!("{error}"));
    assert_eq!(glycans.len(), 1);
    assert_eq!(
        crabwurcs::core::write_wurcs(&glycans[0].graph).unwrap(),
        "WURCS=2.0/2,2,1/[a2122h-1b_1-5_2*NCC/3=O][a2112h-1b_1-5]/1-2/a4-b1"
    );

    let directory = std::env::temp_dir().join(format!("crabwurcs-pdb-{}", std::process::id()));
    std::fs::create_dir_all(&directory).unwrap();
    let path = directory.join("beta.pdb");
    std::fs::write(&path, &pdb).unwrap();
    assert_eq!(
        crabwurcs::extract_glycans_from_file(&path).unwrap().len(),
        1
    );
    std::fs::remove_dir_all(directory).unwrap();
}

/// Opt-in check on the Level 1 files of a GlycoShape static build, as
/// distributed (free-text REMARK banner, five models, every residue an ATOM
/// record, GLYCAM codes such as `0FA` for D-Fuc):
///
/// ```text
/// CRABWURCS_GLYCOSHAPE_STATIC=/path/to/static \
///   cargo test -p crabwurcs --test glycoshape_regressions -- --ignored
/// ```
#[test]
#[ignore = "needs a GlycoShape static build (CRABWURCS_GLYCOSHAPE_STATIC)"]
fn glycoshape_level1_pdb_files_are_read_as_their_names() {
    let root = std::path::PathBuf::from(
        std::env::var("CRABWURCS_GLYCOSHAPE_STATIC")
            .expect("set CRABWURCS_GLYCOSHAPE_STATIC to a GlycoShape static build"),
    );
    for (entry, variant, name) in [
        (
            "GS00002",
            "alpha",
            "LFucpa1-2DGalpa1-3DGalpb1-3DGalpNAca1-OH",
        ),
        (
            "GS00025",
            "beta",
            "DFucpa1-2DGalpb1-3DGlcpNAcb1-3DGalpb1-4DGlcpb1-OH",
        ),
        ("GS00237", "alpha", "DNeup5Ac[9A]a2-6DGalpb1-4DGlcpNAca1-OH"),
        (
            "GS00386",
            "beta",
            "DGalpNAcb1-4[LFucpa1-3]DGlcpNAc[6PC]b1-2DManpa1-3[DManpa1-6]DManpb1-4DGlcpNAcb1-4[LFucpa1-6]DGlcpNAcb1-OH",
        ),
        ("GS00472", "beta", "LArafb1-2LArafb1-2LArafb1-OH"),
        ("GS00502", "alpha", "DFrufb2-1DFrufb2-1DFrufa2-OH"),
        ("GS00893", "alpha", "DAllpa1-OH"),
    ] {
        let path = root
            .join(entry)
            .join("output/level_1/PDB")
            .join(format!("{variant}.pdb"));
        let glycans = crabwurcs::extract_glycans_from_file(&path)
            .unwrap_or_else(|error| panic!("{}: {error}", path.display()));
        assert_eq!(glycans.len(), 1, "{entry} {variant}");
        assert_eq!(
            crabwurcs::core::write_wurcs(&glycans[0].graph).unwrap(),
            wurcs_of(name, Format::Glycam),
            "{entry} {variant}"
        );
    }
}

/// `convert` computes: inputs that are records of the bundled GlycoShape
/// table are no longer answered with the table's strings (`G60371D-N` for
/// the WURCS of GS00728, a damaged string for GS00430, the GLYCAM name for
/// GS00954).
#[test]
fn conversion_of_table_records_is_computed() {
    for (entry, wurcs, expected) in [
        (
            "GS00728",
            "WURCS=2.0/8,23,22/[u2122h_2*NCC/3=O][a2122h-1b_1-5_2*NCC/3=O][a1122h-1b_1-5][a1122h-1a_1-5][a2112h-1b_1-5][Aad21122h-2a_2-6_5*NCC/3=O][a2112h-1b_1-5_2*NCC/3=O][a1221m-1a_1-5]/1-2-3-4-2-5-6-7-2-5-6-7-2-4-2-5-6-7-2-5-6-7-8/a4-b1_a6-w1_b4-c1_c3-d1_c4-m1_c6-n1_d2-e1_d4-i1_e4-f1_f3-g2_f4-h1_i4-j1_j3-k2_j4-l1_n2-o1_n6-s1_o4-p1_p3-q2_p4-r1_s4-t1_t3-u2_t4-v1",
            "Neu5Ac(a2-3)[GalNAc(b1-4)]Gal(b1-4)GlcNAc(b1-2)[Neu5Ac(a2-3)[GalNAc(b1-4)]Gal(b1-4)GlcNAc(b1-4)]Man(a1-3)[GlcNAc(b1-4)][Neu5Ac(a2-3)[GalNAc(b1-4)]Gal(b1-4)GlcNAc(b1-2)[Neu5Ac(a2-3)[GalNAc(b1-4)]Gal(b1-4)GlcNAc(b1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)[Fuc(a1-6)]GlcNAc",
        ),
        (
            "GS00430",
            "WURCS=2.0/3,3,2/[uxxxxm_2*NCC/3=O_4*NCC/3=O][a2112h-1a_1-5][a2112h-1b_1-5]/1-2-3/a3-b1_b4-c1",
            "Gal(b1-4)Gal(a1-3)dHexNAc4NAc",
        ),
        (
            "GS00954",
            "WURCS=2.0/3,7,6/[a2122m-1b_1-5_2*NCC/3=O_4*NCC/3=O][a2112h-1a_1-5_2*NCC/3=O][a2122h-1b_1-5]/1-2-2-2-3-2-2/a3-b1_b4-c1_c4-d1_d3-e1_d4-f1_f4-g1",
            "Glc(b1-3)[GalNAc(a1-4)GalNAc(a1-4)]GalNAc(a1-4)GalNAc(a1-4)GalNAc(a1-3)Bac2Ac4Ac(b1-",
        ),
    ] {
        let iupac = convert(wurcs, Format::Wurcs, Format::IupacCondensed).unwrap();
        assert_eq!(iupac, expected, "{entry}");
        assert_eq!(
            wurcs_of(&iupac, Format::IupacCondensed),
            normalize_wurcs(wurcs).unwrap(),
            "{entry}: {iupac}"
        );
    }
    // A GlyTouCan accession is not a notation.
    assert!(convert("G60371D-N", Format::Auto, Format::Wurcs).is_err());
}
