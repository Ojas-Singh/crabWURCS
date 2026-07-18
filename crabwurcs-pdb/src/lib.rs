use crabwurcs_core::{
    AnomericSymbol, CarbonPosition, Linkage, Monosaccharide, ResidueGraph, RingClosure,
};
use thiserror::Error;

#[derive(Debug, Error)]
pub enum PdbError {
    #[error("failed to parse structure file: {0}")]
    ParseError(String),

    #[error("no carbohydrate residues found in structure")]
    NoGlycansFound,

    #[error(transparent)]
    Core(#[from] crabwurcs_core::CoreError),
}

pub type PdbResult<T> = Result<T, PdbError>;

#[derive(Debug, Clone)]
pub struct ExtractedGlycan {
    pub attachment_site: Option<String>,
    pub graph: ResidueGraph,
}

const SUGAR_RESIDUE_NAMES: &[&str] = &[
    "NAG", "NDG", "BMA", "MAN", "FUC", "GAL", "GLC", "SIA", "NGN", "XYS", "XYL", "RIP", "RIB",
    "ARA", "LYX", "ALL", "ALT", "GUL", "IDO", "TAL", "GLA", "FRU", "SOR", "PSI", "TAG", "XYU",
    "BGC", "BMX", "FCA", "FCB", "FUL", "GCU", "GCV", "G6D", "GIV", "GL0", "GNS", "GTR", "GXL",
    "GYC", "GYG", "GYP", "GYS", "HSG", "HSQ", "HSY", "KDM", "KDO", "KDN", "L6S", "LAT", "LFR",
    "MHS", "NGA", "NGC", "NGE", "NGS", "NM6", "OAA", "RAM", "SGA", "SGC", "SGD", "SGN", "SGS",
    "SIB", "SID", "SLB", "SUS", "T6S", "XYP", "XXR", "Z0E", "ZB0", "ZB1", "A2G", "GCS",
];

fn is_sugar_residue(name: &str) -> bool {
    SUGAR_RESIDUE_NAMES.contains(&name)
}

fn residue_to_monosaccharide(name: &str) -> Monosaccharide {
    match name {
        "NAG" | "NGA" => Monosaccharide::new(
            4,
            "2122h".into(),
            vec![],
            RingClosure::Pyranose,
            Some(1),
            Some(5),
            1,
            AnomericSymbol::Beta,
            String::from("x"),
            vec![],
        ),
        "NDG" | "NGC" | "NGE" | "NGN" | "NGS" => Monosaccharide::new(
            4,
            "2122h".into(),
            vec![],
            RingClosure::Pyranose,
            Some(1),
            Some(5),
            1,
            AnomericSymbol::Unknown,
            String::from("x"),
            vec![],
        ),
        "BMA" | "MAN" => Monosaccharide::new(
            4,
            "1221m".into(),
            vec![],
            RingClosure::Pyranose,
            Some(1),
            Some(5),
            1,
            AnomericSymbol::Alpha,
            String::from("x"),
            vec![],
        ),
        "FUC" | "FUL" | "FCA" | "FCB" => Monosaccharide::new(
            4,
            "d122m".into(),
            vec![],
            RingClosure::Pyranose,
            Some(1),
            Some(5),
            1,
            AnomericSymbol::Alpha,
            String::from("x"),
            vec![],
        ),
        "GAL" | "GLA" | "GIV" | "GXL" => Monosaccharide::new(
            4,
            "2112h".into(),
            vec![],
            RingClosure::Pyranose,
            Some(1),
            Some(5),
            1,
            AnomericSymbol::Alpha,
            String::from("x"),
            vec![],
        ),
        "GLC" | "BGC" | "GCS" | "GCU" | "GCV" | "G6D" | "GL0" => Monosaccharide::new(
            4,
            "2122h".into(),
            vec![],
            RingClosure::Pyranose,
            Some(1),
            Some(5),
            1,
            AnomericSymbol::Beta,
            String::from("x"),
            vec![],
        ),
        "SIA" | "SIB" | "SID" | "SLB" | "SUS" => Monosaccharide::new(
            7,
            "d21122h".into(),
            vec![],
            RingClosure::Pyranose,
            Some(2),
            Some(6),
            2,
            AnomericSymbol::Alpha,
            String::from("x"),
            vec![],
        ),
        "XYS" | "XYL" | "XYP" => Monosaccharide::new(
            3,
            "212h".into(),
            vec![],
            RingClosure::Pyranose,
            Some(1),
            Some(5),
            1,
            AnomericSymbol::Beta,
            String::from("x"),
            vec![],
        ),
        "KDM" | "KDO" | "KDN" => Monosaccharide::new(
            4,
            "d212a".into(),
            vec![],
            RingClosure::Pyranose,
            Some(2),
            Some(5),
            2,
            AnomericSymbol::Alpha,
            String::from("x"),
            vec![],
        ),
        "RAM" => Monosaccharide::new(
            4,
            "d122m".into(),
            vec![],
            RingClosure::Pyranose,
            Some(1),
            Some(5),
            1,
            AnomericSymbol::Alpha,
            String::from("x"),
            vec![],
        ),
        _ => Monosaccharide::new(
            4,
            "xxxxh".into(),
            vec![],
            RingClosure::Pyranose,
            Some(1),
            Some(5),
            1,
            AnomericSymbol::Unknown,
            String::from("x"),
            vec![],
        ),
    }
}

pub fn extract_glycans_from_file(path: &std::path::Path) -> PdbResult<Vec<ExtractedGlycan>> {
    use pdbtbx::ReadOptions;

    let path_str = path.to_string_lossy();
    let (pdb, _errors) = ReadOptions::default()
        .set_level(pdbtbx::StrictnessLevel::Loose)
        .read(path_str.as_ref())
        .map_err(|e| PdbError::ParseError(format!("cannot parse file: {:?}", e)))?;

    extract_glycans_from_pdb(&pdb)
}

pub fn extract_glycans_from_str(contents: &str, is_mmcif: bool) -> PdbResult<Vec<ExtractedGlycan>> {
    use pdbtbx::{Format, ReadOptions, StrictnessLevel};
    use std::io::BufReader;

    let format = if is_mmcif { Format::Mmcif } else { Format::Pdb };
    let reader = BufReader::new(contents.as_bytes());
    let (pdb, _errors) = ReadOptions::default()
        .set_level(StrictnessLevel::Loose)
        .set_format(format)
        .read_raw(reader)
        .map_err(|e| PdbError::ParseError(format!("cannot parse: {:?}", e)))?;

    extract_glycans_from_pdb(&pdb)
}

fn extract_glycans_from_pdb(pdb: &pdbtbx::PDB) -> PdbResult<Vec<ExtractedGlycan>> {
    let mut sugar_residues: Vec<(String, isize)> = Vec::new();

    for residue in pdb.residues() {
        if let Some(name) = residue.name() {
            if is_sugar_residue(name) {
                let (seq_id, _insertion) = residue.id();
                sugar_residues.push((name.to_string(), seq_id));
            }
        }
    }

    if sugar_residues.is_empty() {
        return Ok(vec![]);
    }

    let mut graph = ResidueGraph::new();
    let mut added_nodes = Vec::new();

    for (name, _seq_id) in &sugar_residues {
        let mono = residue_to_monosaccharide(name);
        let idx = graph.add_residue(mono);
        added_nodes.push(idx);
    }

    for i in 1..added_nodes.len() {
        let parent_idx = added_nodes[i - 1];
        let child_idx = added_nodes[i];
        graph.add_linkage(
            parent_idx,
            child_idx,
            Linkage::new(CarbonPosition(4), CarbonPosition(1)),
        );
    }

    if !added_nodes.is_empty() {
        graph.set_root(added_nodes[0]);
    }

    let site = sugar_residues.first().map(|(n, s)| format!("{}/{}", n, s));

    Ok(vec![ExtractedGlycan {
        attachment_site: site,
        graph,
    }])
}

#[cfg(test)]
mod tests {
    use super::*;

    #[allow(clippy::too_many_arguments)]
    fn make_hetatm_line(
        serial: u32,
        atom_name: &str,
        res_name: &str,
        chain: &str,
        seq: u32,
        x: f64,
        y: f64,
        z: f64,
    ) -> String {
        format!(
            "HETATM{:5} {:<4}{:1}{:3} {}{:4}{:1}   {:8.3}{:8.3}{:8.3}{:6.2}{:6.2}          {:>2}{:2}",
            serial, atom_name, "", res_name, chain, seq, "", x, y, z, 1.0f64, 0.0f64, "C", ""
        )
    }

    #[test]
    fn test_sugar_residue_detection() {
        assert!(is_sugar_residue("NAG"));
        assert!(is_sugar_residue("MAN"));
        assert!(is_sugar_residue("FUC"));
        assert!(!is_sugar_residue("HOH"));
        assert!(!is_sugar_residue("ALA"));
    }

    #[test]
    fn test_residue_to_monosaccharide() {
        let nag = residue_to_monosaccharide("NAG");
        assert_eq!(nag.backbone_length, 4);

        let sia = residue_to_monosaccharide("SIA");
        assert_eq!(sia.backbone_length, 7);
    }

    #[test]
    fn test_parse_minimal_pdb_with_sugar() {
        let lines = [
            "HEADER    TEST                                                            END"
                .to_string(),
            make_hetatm_line(1, "C1", "NAG", "A", 1, -1.0, 0.0, 0.0),
            make_hetatm_line(2, "C2", "NAG", "A", 1, 0.0, 0.0, 0.0),
            make_hetatm_line(3, "C1", "BMA", "A", 2, 2.0, 0.0, 0.0),
            "END".to_string(),
        ];
        let pdb_str = lines.join("\n") + "\n";

        let result = extract_glycans_from_str(&pdb_str, false);
        assert!(result.is_ok(), "Failed: {:?}", result.err());
        let glycans = result.unwrap();
        assert_eq!(glycans.len(), 1);
        assert!(glycans[0].graph.node_count() >= 2);
    }

    #[test]
    fn test_no_sugar_pdb() {
        let lines = [
            "HEADER    TEST                                                            END"
                .to_string(),
            make_hetatm_line(1, "CA", "ALA", "A", 1, -1.0, 0.0, 0.0),
            "END".to_string(),
        ];
        let pdb_str = lines.join("\n") + "\n";

        let result = extract_glycans_from_str(&pdb_str, false);
        assert!(result.is_ok(), "Failed: {:?}", result.err());
        assert!(result.unwrap().is_empty());
    }
}
