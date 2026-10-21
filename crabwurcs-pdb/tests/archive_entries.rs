//! wwPDB archive entries: the residue a glycan is attached to, read from
//! `_struct_conn` and `LINK` records and, when a file records no connection,
//! from distances; space-group items pdbtbx rejects; anomers deposited as two
//! residues on the same atoms; thioglycosides, a fluorinated C6 and a
//! cyclodextrin; residues IUPAC condensed has no name for.
//!
//! The fixtures (`data/*.cif`, `data/*.pdb`) are cut from archive files; each
//! says where it comes from and how to regenerate it.

use crabwurcs_core::write_wurcs;
use crabwurcs_iupac::{IupacError, parse_iupac_condensed, write_iupac_condensed};
use crabwurcs_pdb::{PdbError, PdbReadOptions, read_glycans_from_str};

const N_GLYCANS_1EO8_CIF: &str = include_str!("data/1eo8_glycosylation.cif");
const N_GLYCANS_1EO8_PDB: &str = include_str!("data/1eo8_glycosylation.pdb");
const C_MANNOSIDE_5LFR: &str = include_str!("data/5lfr_c_mannosylation.cif");
const O_GLYCAN_1G1S: &str = include_str!("data/1g1s_o_glycosylation.cif");
const CALCIUM_1K9I: &str = include_str!("data/1k9i_calcium.pdb");
const SPACE_GROUP_31NP: &str = include_str!("data/31np_space_group.cif");
const ANOMER_RESIDUES_4D52: &str = include_str!("data/4d52_anomer_residues.cif");
const ANOMER_ALTERNATIVES_1AXZ: &str = include_str!("data/1axz_anomer_alternatives.cif");
const THIODIGALACTOSIDE_1A78: &str = include_str!("data/1a78_thiodigalactoside.cif");
const THIO_LINK_5M77: &str = include_str!("data/5m77_thio_link_c6.cif");
const FLUORO_C6_9I4K: &str = include_str!("data/9i4k_fluoro_c6.cif");
const CYCLODEXTRIN_3CK7: &str = include_str!("data/3ck7_cyclodextrin.cif");
const SELENOFUCOSIDE_5H47: &str = include_str!("data/5h47_selenofucoside.cif");
const C_FUCOSIDE_5D2A: &str = include_str!("data/5d2a_c_fucoside.cif");

/// `(attachment site, WURCS)` of every glycan, in file order.
fn glycans(text: &str, is_mmcif: bool) -> Vec<(String, String)> {
    read_glycans_from_str(text, is_mmcif, &PdbReadOptions::default())
        .unwrap()
        .glycans
        .iter()
        .map(|glycan| {
            (
                glycan.attachment_site.clone().unwrap_or_default(),
                write_wurcs(&glycan.graph).unwrap(),
            )
        })
        .collect()
}

fn sites(text: &str, is_mmcif: bool) -> Vec<String> {
    glycans(text, is_mmcif)
        .into_iter()
        .map(|(site, _)| site)
        .collect()
}

/// The text without its `LINK` records, and with its `_struct_conn` loop
/// turned into a category nothing reads.
fn without_connections(text: &str) -> String {
    text.lines()
        .filter(|line| !line.starts_with("LINK"))
        .map(|line| format!("{}\n", line.replace("_struct_conn.", "_unread_conn.")))
        .collect()
}

const GLCNAC: &str = "WURCS=2.0/1,1,0/[a2122h-1b_1-5_2*NCC/3=O]/1/";
const MAN_GLCNAC2: &str =
    "WURCS=2.0/2,3,2/[a2122h-1b_1-5_2*NCC/3=O][a1122h-1b_1-5]/1-1-2/a4-b1_b4-c1";

#[test]
fn struct_conn_names_the_asparagine_of_an_n_glycan() {
    assert_eq!(
        glycans(N_GLYCANS_1EO8_CIF, true),
        vec![
            ("A/ASN/81".to_string(), GLCNAC.to_string()),
            ("A/ASN/165".to_string(), MAN_GLCNAC2.to_string()),
        ]
    );
}

#[test]
fn link_records_give_what_struct_conn_gives() {
    // The PDB fixture has LINK records and no CONECT records.
    assert!(
        !N_GLYCANS_1EO8_PDB
            .lines()
            .any(|line| line.starts_with("CONECT"))
    );
    assert_eq!(
        glycans(N_GLYCANS_1EO8_PDB, false),
        glycans(N_GLYCANS_1EO8_CIF, true)
    );
}

#[test]
fn without_connection_records_the_attachment_is_read_from_distances() {
    for (text, is_mmcif) in [(N_GLYCANS_1EO8_CIF, true), (N_GLYCANS_1EO8_PDB, false)] {
        let bare = without_connections(text);
        assert!(!bare.contains("_struct_conn.") && !bare.contains("\nLINK"));
        assert_eq!(glycans(&bare, is_mmcif), glycans(text, is_mmcif));
    }
}

#[test]
fn a_c_mannoside_is_attached_to_its_tryptophan() {
    assert_eq!(sites(C_MANNOSIDE_5LFR, true), vec!["A/TRP/22"]);
    assert_eq!(
        sites(&without_connections(C_MANNOSIDE_5LFR), true),
        vec!["A/TRP/22"]
    );
}

#[test]
fn an_o_glycan_is_attached_to_its_threonine() {
    assert_eq!(sites(O_GLYCAN_1G1S, true), vec!["C/THR/616"]);
    assert_eq!(
        sites(&without_connections(O_GLYCAN_1G1S), true),
        vec!["C/THR/616"]
    );
}

#[test]
fn a_contact_is_no_attachment_when_the_file_records_connections() {
    // The file records its connections, but not the one of this glycan: the
    // glycan stays unattached, however close the asparagine.
    let text = N_GLYCANS_1EO8_CIF
        .lines()
        .filter(|line| !line.starts_with("covale1 "))
        .map(|line| format!("{line}\n"))
        .collect::<String>();
    assert_eq!(sites(&text, true), vec!["A/NAG/440", "A/ASN/165"]);
}

#[test]
fn a_metal_coordination_names_no_attachment() {
    // LINK and CONECT records bond the calcium ion to O3 and O4 of MAN K2.
    assert_eq!(
        glycans(CALCIUM_1K9I, false),
        vec![(
            "K/MAN/1".to_string(),
            "WURCS=2.0/2,5,4/[a1122h-1a_1-5][a2122h-1b_1-5_2*NCC/3=O]/1-1-2-1-2/a3-b1_a6-d1_b2-c1_d2-e1"
                .to_string()
        )]
    );
}

#[test]
fn a_disagreeing_space_group_does_not_stop_the_reading() {
    let reading = read_glycans_from_str(SPACE_GROUP_31NP, true, &PdbReadOptions::default())
        .expect("the space-group items are not needed");
    assert_eq!(reading.glycans.len(), 1);
    assert_eq!(reading.glycans[0].graph.node_count(), 2);
    assert_eq!(
        reading.glycans[0].attachment_site.as_deref(),
        Some("B/A1KB4/1")
    );
}

#[test]
fn an_unknown_space_group_or_cell_does_not_stop_the_reading() {
    // pdbtbx panics on `?` for these items.
    for (item, value) in [
        ("_symmetry.Int_Tables_number                5", "?"),
        ("_cell.length_a                     61.482", "?"),
    ] {
        let tag = item.split_whitespace().next().unwrap();
        let text = SPACE_GROUP_31NP.replace(item, &format!("{tag} {value}"));
        assert_ne!(text, SPACE_GROUP_31NP);
        let reading = read_glycans_from_str(&text, true, &PdbReadOptions::default()).unwrap();
        assert_eq!(reading.glycans.len(), 1);
    }
}

#[test]
fn other_parse_errors_still_stop_the_reading() {
    // A record type pdbtbx refuses, in a file whose space group disagrees.
    let broken = SPACE_GROUP_31NP.replacen("\nHETATM", "\nHETATX", 1);
    match read_glycans_from_str(&broken, true, &PdbReadOptions::default()) {
        Err(PdbError::ParseError(message)) => {
            assert!(message.contains("Atom type not correct"), "{message}");
        }
        other => panic!("expected a parse error, got {other:?}"),
    }
}

/// `(WURCS, IUPAC condensed)` of the only glycan of an mmCIF text, the IUPAC
/// read back to the same WURCS.
fn single_glycan(text: &str) -> (String, String) {
    let reading = read_glycans_from_str(text, true, &PdbReadOptions::default()).unwrap();
    assert_eq!(reading.glycans.len(), 1, "{:?}", reading.glycans);
    let graph = &reading.glycans[0].graph;
    let wurcs = write_wurcs(graph).unwrap();
    let iupac = write_iupac_condensed(graph).unwrap();
    assert_eq!(
        write_wurcs(&parse_iupac_condensed(&iupac).unwrap()).unwrap(),
        wurcs,
        "{iupac}"
    );
    (wurcs, iupac)
}

#[test]
fn an_anomer_deposited_as_a_second_residue_is_read_once() {
    // 4D52: GXL (alpha, occupancy 0.2) and GIV (beta, 0.8) share the ring
    // atoms. Read together, every centre was unreadable (`Res4`).
    assert_eq!(
        single_glycan(ANOMER_RESIDUES_4D52),
        (
            "WURCS=2.0/1,1,0/[a1221h-1b_1-5]/1/".to_string(),
            "L-Gal(b1-".to_string()
        )
    );
    // 1AXZ: GLA in alternative location A, GAL in B; they read as a
    // 1-1-linked "GalA(b1-1)Gal(a1-".
    assert_eq!(
        single_glycan(ANOMER_ALTERNATIVES_1AXZ),
        (
            "WURCS=2.0/1,1,0/[a2112h-1a_1-5]/1/".to_string(),
            "Gal(a1-".to_string()
        )
    );
}

#[test]
fn a_thioglycosidic_bond_keeps_its_sulfur() {
    // Thiodigalactoside (YIO, 1-thio-beta-D-galactose, and GAL): the S
    // bridge was read as an O-glycosidic bond.
    assert_eq!(
        single_glycan(THIODIGALACTOSIDE_1A78),
        (
            "WURCS=2.0/1,2,1/[a2112h-1b_1-5]/1-1/a1-b1*S*".to_string(),
            "Gal(b1-SH-1)Gal(b1-".to_string()
        )
    );
    // 5M77: a mannose S-linked to C6 of a mannose, whose C6 (CH2-S) was
    // unreadable.
    assert_eq!(
        single_glycan(THIO_LINK_5M77),
        (
            "WURCS=2.0/1,2,1/[a1122h-1a_1-5]/1-1/a6-b1*S*".to_string(),
            "Man(a1-SH-6)Man(a1-".to_string()
        )
    );
}

#[test]
fn a_terminal_carbon_with_a_fluorine_is_a_hydroxymethyl_position() {
    let reading = read_glycans_from_str(FLUORO_C6_9I4K, true, &PdbReadOptions::default()).unwrap();
    assert!(reading.warnings.is_empty(), "{:?}", reading.warnings);
    assert_eq!(
        write_wurcs(&reading.glycans[0].graph).unwrap(),
        "WURCS=2.0/2,2,1/[a2122h-1b_1-5_1*OC_2*NCC/3=O_6*F][a2112h-1b_1-5]/1-2/a4-b1"
    );
}

#[test]
fn a_cyclodextrin_is_a_closed_ring() {
    let (wurcs, iupac) = single_glycan(CYCLODEXTRIN_3CK7);
    assert_eq!(
        wurcs,
        "WURCS=2.0/1,6,6/[a2122h-1a_1-5]/1-1-1-1-1-1/a4-b1_b4-c1_c4-d1_d4-e1_e4-f1_f4-a1"
    );
    // The closing linkage was dropped from the IUPAC.
    assert_eq!(
        iupac,
        "4)Glc(a1-4)Glc(a1-4)Glc(a1-4)Glc(a1-4)Glc(a1-4)Glc(a1-"
    );
}

#[test]
fn residues_iupac_condensed_has_no_name_for_are_an_error() {
    // FSW (a 2-Se-methyl-2-seleno-L-fucoside: C2 is not read through Se) and
    // ZDC (a C-fucoside, read as a ketose) were written `Res41Me` and
    // `Res7`.
    for text in [SELENOFUCOSIDE_5H47, C_FUCOSIDE_5D2A] {
        let reading = read_glycans_from_str(text, true, &PdbReadOptions::default()).unwrap();
        assert_eq!(reading.glycans.len(), 1);
        match write_iupac_condensed(&reading.glycans[0].graph) {
            Err(IupacError::UnrepresentableInFormat { residue, .. }) => {
                assert!(residue.contains("no IUPAC condensed name"), "{residue}");
            }
            other => panic!("expected an error, got {other:?}"),
        }
    }
}
