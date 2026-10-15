//! crabWURCS: a pure-Rust toolkit for glycan notation conversion,
//! chemical-structure interop, and SNFG rendering.
//!
//! This crate is a thin facade — it does no work itself, it just
//! re-exports the workspace's other crates under one name so a consumer
//! can `use crabwurcs::...` instead of depending on each sub-crate
//! individually. All real logic lives in:
//! - [`core`] — the WURCS grammar and the shared [`core::ResidueGraph`] model
//! - [`iupac`] — IUPAC condensed/extended notation
//! - [`mol`] — MOL/SDF/SMILES chemical structure interop (MolWURCS equivalent)
//! - [`pdb`] — glycan extraction from PDB/mmCIF structures
//! - [`snfg`] — SNFG SVG rendering (seq2snfg equivalent)
//!
//! Conversions are computed from the parsed residue graph. The facade also
//! carries a table of the notations of 942 GlycoShape records (938
//! molecularly specified structures plus four notation-only edge cases); it is
//! consulted only when asked, through [`convert_with_table`], either as a
//! fallback for inputs that cannot be converted or as a check of a computed
//! structure. SMILES is written and read by the pure-Rust molecular backend.
//!
//! The sections below — generated from the renderer's own output — document
//! every supported monosaccharide symbol, the SNFG rendering rules, and motif
//! highlighting.

#![doc = include_str!("../docs/supported-monosaccharides.md")]
#![doc = include_str!("../docs/rendering.md")]
#![doc = include_str!("../docs/motif-highlighting.md")]

use std::collections::HashMap;
use std::sync::OnceLock;

pub use crabwurcs_core as core;
pub use crabwurcs_iupac as iupac;
pub use crabwurcs_mol as mol;
pub use crabwurcs_pdb as pdb;
pub use crabwurcs_snfg as snfg;

// Convenience re-exports of the most commonly needed types, so simple
// consumers don't need to reach into `core::`.
pub use crabwurcs_core::{
    CoreError, CoreResult, MotifError, MotifMatch, ResidueGraph, ResidueKind, classify_residue,
    find_motif_matches, normalize_wurcs, residue_from_kind, standardize_wurcs,
    write_wurcs_canonical,
};
pub use crabwurcs_iupac::write_iupac_condensed_canonical;
pub use crabwurcs_pdb::{
    ExtractedGlycan, ExtractedGlycanWithProvenance, NamePolicy, PdbReadOptions, PdbReading,
    PdbResidueReference, PdbWarning, PdbWarningKind, ResidueReading, StereoCentre,
    extract_glycans_from_file, extract_glycans_from_str, extract_glycans_with_provenance_from_file,
    extract_glycans_with_provenance_from_str, read_glycans_from_file, read_glycans_from_str,
};
pub use crabwurcs_snfg::{
    HighlightSelection, RenderOptions, SnfgError, SnfgResult, SourceNotation, render_png,
    render_png_with_motifs, render_png_with_options, render_png_with_selection, render_svg,
    render_svg_with_motifs, render_svg_with_options, render_svg_with_selection, render_symbol_svg,
};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Format {
    Auto,
    Wurcs,
    IupacCondensed,
    IupacExtended,
    Glycam,
    Smiles,
}

#[derive(Debug, thiserror::Error)]
pub enum ConversionError {
    #[error(transparent)]
    Core(#[from] core::CoreError),
    #[error(transparent)]
    Iupac(#[from] iupac::IupacError),
    #[error(transparent)]
    Molecule(#[from] mol::MolError),
    #[error("automatic format detection cannot be used as an output format")]
    AutoOutput,
}

pub type ConversionResult<T> = Result<T, ConversionError>;

/// One exact GlycoShape archetype record. Keeping this as a tab-separated
/// compile-time table avoids parsing the much larger metadata JSON at runtime.
#[derive(Clone, Copy)]
struct CorpusRecord<'a> {
    wurcs: &'a str,
    iupac: &'a str,
    iupac_extended: &'a str,
    glycam: &'a str,
    smiles: &'a str,
}

impl CorpusRecord<'static> {
    fn notation(self, format: Format) -> Option<&'static str> {
        let value = match format {
            Format::Wurcs => self.wurcs,
            Format::IupacCondensed => self.iupac,
            Format::IupacExtended => self.iupac_extended,
            Format::Glycam => self.glycam,
            Format::Smiles => self.smiles,
            Format::Auto => return None,
        };
        (!value.is_empty()).then_some(value)
    }
}

fn corpus_records() -> impl Iterator<Item = CorpusRecord<'static>> {
    const DATA: &str = include_str!("../data/glycoshape_notations.tsv");
    const DERIVED: &str = include_str!("../data/glycoshape_derived_notations.tsv");
    const NOTATION_ONLY: &str = include_str!("../data/glycoshape_notation_only.tsv");
    DATA.lines()
        .chain(DERIVED.lines())
        .chain(NOTATION_ONLY.lines())
        .filter_map(|line| {
            let mut fields = line.splitn(5, '\t');
            Some(CorpusRecord {
                wurcs: fields.next()?,
                iupac: fields.next()?,
                iupac_extended: fields.next()?,
                glycam: fields.next()?,
                smiles: fields.next()?,
            })
        })
}

fn corpus_wurcs_normalizations() -> &'static HashMap<&'static str, String> {
    static NORMALIZED: OnceLock<HashMap<&'static str, String>> = OnceLock::new();
    NORMALIZED.get_or_init(|| {
        corpus_records()
            .filter_map(|record| {
                let graph = core::parse_wurcs(record.wurcs).ok()?;
                let normalized = core::write_wurcs(&graph).ok()?;
                Some((record.wurcs, normalized))
            })
            .collect()
    })
}

fn normalized_wurcs_corpus() -> &'static HashMap<String, CorpusRecord<'static>> {
    static NORMALIZED: OnceLock<HashMap<String, CorpusRecord<'static>>> = OnceLock::new();
    NORMALIZED.get_or_init(|| {
        corpus_records()
            .filter_map(|record| {
                let normalized = corpus_wurcs_normalizations().get(record.wurcs)?;
                Some((normalized.clone(), record))
            })
            .collect()
    })
}

fn corpus_record_for_input(input: &str, format: Format) -> Option<CorpusRecord<'static>> {
    let input = input.trim();
    if input.is_empty() {
        return None;
    }
    if format == Format::Auto {
        // Any column: a fallback must not depend on detecting the format.
        if input.starts_with("WURCS=") {
            return corpus_record_for_input(input, Format::Wurcs);
        }
        return corpus_records().find(|record| {
            [
                record.iupac,
                record.iupac_extended,
                record.glycam,
                record.smiles,
            ]
            .contains(&input)
        });
    }
    if format == Format::Wurcs {
        if let Some(record) = corpus_records().find(|record| record.wurcs == input) {
            return Some(record);
        }
        let normalized = core::normalize_wurcs(input).ok()?;
        return normalized_wurcs_corpus().get(&normalized).copied();
    }
    corpus_records().find(|record| match format {
        Format::Wurcs => unreachable!("handled above"),
        Format::IupacCondensed => record.iupac == input,
        Format::IupacExtended => record.iupac_extended == input,
        Format::Glycam => record.glycam == input,
        Format::Smiles => record.smiles == input,
        Format::Auto => false,
    })
}

/// Whether `value` holds a GLYCAM linkage written without parentheses
/// (`a1-4`, `b2-1`, `a2-3`).
fn glycam_linkage(value: &str) -> bool {
    value.as_bytes().windows(4).any(|window| {
        matches!(window[0], b'a' | b'b')
            && window[1].is_ascii_digit()
            && window[2] == b'-'
            && (window[3].is_ascii_digit() || window[3] == b'?' || window[3] == b'O')
    })
}

/// Guess the notation of `input` from its syntax.
pub fn detect_format(input: &str) -> Format {
    let value = input.trim();
    if value.starts_with("WURCS=") {
        Format::Wurcs
    } else if value.contains('→')
        || value.contains("α-")
        || value.contains("β-")
        || (!value.contains('(') && (value.starts_with("D-") || value.starts_with("L-")))
    {
        Format::IupacExtended
    } else if value.contains('@')
        || value.contains('#')
        || value.contains("C=")
        || value.starts_with("OC[")
        || value.starts_with("C[")
    {
        Format::Smiles
    } else if !value.contains('(')
        && (glycam_linkage(value) || value.starts_with('D') || value.starts_with('L'))
    {
        Format::Glycam
    } else {
        Format::IupacCondensed
    }
}

pub fn parse_notation(input: &str, format: Format) -> ConversionResult<ResidueGraph> {
    let format = if format == Format::Auto {
        detect_format(input)
    } else {
        format
    };
    Ok(match format {
        Format::Auto => unreachable!(),
        Format::Wurcs => core::parse_wurcs(input)?,
        Format::IupacCondensed => iupac::parse_iupac_condensed(input)?,
        Format::IupacExtended => iupac::parse_iupac_extended(input)?,
        Format::Glycam => iupac::parse_glycam(input)?,
        Format::Smiles => mol::wurcs_from_molecule(input, mol::ChemFormat::Smiles)?,
    })
}

pub fn write_notation(graph: &ResidueGraph, format: Format) -> ConversionResult<String> {
    Ok(match format {
        Format::Auto => return Err(ConversionError::AutoOutput),
        Format::Wurcs => core::write_wurcs(graph)?,
        Format::IupacCondensed => iupac::write_iupac_condensed(graph)?,
        Format::IupacExtended => iupac::write_iupac_extended(graph)?,
        Format::Glycam => iupac::write_glycam(graph)?,
        Format::Smiles => mol::molecule_from_wurcs(graph, mol::ChemFormat::Smiles)?,
    })
}

/// Convert `input` from one notation to another by parsing it into the
/// residue graph and writing that graph. The bundled GlycoShape table is not
/// consulted; see [`convert_with_table`].
pub fn convert(input: &str, from: Format, to: Format) -> ConversionResult<String> {
    if to == Format::Auto {
        return Err(ConversionError::AutoOutput);
    }
    write_notation(&parse_notation(input, from)?, to)
}

/// How [`convert_with_table`] uses the bundled table of GlycoShape notations.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum TableMode {
    /// Compute the conversion; never consult the table (as [`convert`]).
    #[default]
    Off,
    /// Compute the conversion; only when the input cannot be converted,
    /// answer from the table record that holds the input, if there is one.
    Fallback,
    /// Compute the conversion and require the structure read from the input
    /// to be the one of the table record that holds the input, if any.
    Check,
}

#[derive(Debug, thiserror::Error)]
pub enum TableConversionError {
    #[error(transparent)]
    Conversion(#[from] ConversionError),
    #[error(
        "the input reads as {computed}, but the bundled GlycoShape record for it holds {table}"
    )]
    Disagreement { computed: String, table: String },
}

/// The notation in `to` of the bundled GlycoShape record that holds `input`
/// in `from` (`Format::Auto`: in any column), without converting anything.
/// WURCS answers are normalized.
pub fn table_lookup(input: &str, from: Format, to: Format) -> Option<String> {
    let record = corpus_record_for_input(input, from)?;
    let value = record.notation(to)?;
    if to == Format::Wurcs {
        return corpus_wurcs_normalizations()
            .get(value)
            .cloned()
            .or_else(|| Some(value.to_owned()));
    }
    Some(value.to_owned())
}

/// [`convert`], with the bundled GlycoShape table used as `mode` says.
pub fn convert_with_table(
    input: &str,
    from: Format,
    to: Format,
    mode: TableMode,
) -> Result<String, TableConversionError> {
    if to == Format::Auto {
        return Err(ConversionError::AutoOutput.into());
    }
    match mode {
        TableMode::Off => Ok(convert(input, from, to)?),
        TableMode::Fallback => match convert(input, from, to) {
            Ok(output) => Ok(output),
            Err(error) => table_lookup(input, from, to).ok_or(error.into()),
        },
        TableMode::Check => {
            let graph = parse_notation(input, from)?;
            if let Some(table) = table_lookup(input, from, Format::Wurcs) {
                let computed = core::write_wurcs(&graph).map_err(ConversionError::from)?;
                let table = core::normalize_wurcs(&table).unwrap_or(table);
                if computed != table {
                    return Err(TableConversionError::Disagreement { computed, table });
                }
            }
            Ok(write_notation(&graph, to)?)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const WURCS: &str = "WURCS=2.0/4,4,3/[u2112h_2*NCC/3=O][a2112h-1b_1-5][a2112h-1a_1-5][a1221m-1a_1-5]/1-2-3-4/a3-b1_b3-c1_c2-d1";
    const IUPAC: &str = "Fuc(a1-2)Gal(a1-3)Gal(b1-3)GalNAc";

    #[test]
    fn corpus_iupac_to_smiles_and_back_to_wurcs() {
        let smiles = convert(IUPAC, Format::IupacCondensed, Format::Smiles).unwrap();
        assert!(smiles.contains('@'));
        assert_eq!(
            convert(&smiles, Format::Smiles, Format::Wurcs).unwrap(),
            WURCS
        );
    }

    #[test]
    fn wurcs_conversion_normalizes_reordered_input() {
        let reordered = "WURCS=2.0/2,2,1/[a2112h-1b_1-5][u2122h]/1-2/b3-a1";
        assert_eq!(
            convert(reordered, Format::Wurcs, Format::Wurcs).unwrap(),
            "WURCS=2.0/2,2,1/[u2122h][a2112h-1b_1-5]/1-2/a3-b1"
        );
    }

    #[test]
    fn unknown_smiles_uses_the_chemistry_backend() {
        let result = parse_notation("C1CC1", Format::Smiles);
        assert!(matches!(
            result,
            Err(ConversionError::Molecule(mol::MolError::NoGlycanFound))
        ));
    }

    /// The only table notations that do not compute to their record's WURCS.
    /// Each is a defect of the stored notation, not of the readers.
    const TABLE_RECORDS_THAT_DISAGREE: &[(&str, &str)] = &[
        // GlyTouCan accessions stored where a name belongs (GS00722, GS00728):
        // identifiers, which crabWURCS does not look up.
        ("glycam", "G50508SG"),
        // GS00430: a damaged IUPAC string.
        (
            "iupac",
            "Galb1-4)Gala1-3(?24-diacetimido-246-trideoxyhexose)",
        ),
        ("iupac", "G60371D-N"),
        ("glycam", "G60371DN"),
        // GS00954: the archetype name leaves the reducing-end anomer open,
        // the record's WURCS (GlyTouCan G66397DB) declares it beta.
        (
            "glycam",
            "DGlcpb1-3[DGalpNAca1-4DGalpNAca1-4]DGalpNAca1-4DGalpNAca1-4DGalpNAca1-3DBacp[2Ac,4Ac]",
        ),
    ];

    /// Table notations that only the table can resolve.
    fn only_in_the_table(notation: &str) -> bool {
        TABLE_RECORDS_THAT_DISAGREE[..4]
            .iter()
            .any(|(_, value)| *value == notation)
    }

    #[test]
    fn autodetection_parses_every_nonempty_corpus_notation() {
        for record in corpus_records() {
            for notation in [
                record.wurcs,
                record.iupac,
                record.iupac_extended,
                record.glycam,
                record.smiles,
            ] {
                if notation.is_empty() || only_in_the_table(notation) {
                    continue;
                }
                let graph = parse_notation(notation, Format::Auto).unwrap_or_else(|error| {
                    panic!("failed to autodetect/parse {notation}: {error}")
                });
                assert!(graph.node_count() > 0, "{notation}");
            }
        }
    }

    #[test]
    fn notation_only_glycoshape_rows_convert_by_computation() {
        let kdo = convert("D-KDOp", Format::Auto, Format::Wurcs).unwrap();
        assert_eq!(kdo, "WURCS=2.0/1,1,0/[AUd1122h]/1/");

        let bac_iupac = "DGlcpb1-3[DGalpNAca1-4DGalpNAca1-4]DGalpNAca1-4DGalpNAca1-4DGalpNAca1-3DBacp[2Ac,4Ac]b1-OH";
        let bac_glycam = convert(bac_iupac, Format::Auto, Format::Glycam).unwrap();
        assert_eq!(bac_glycam, bac_iupac);
        assert!(
            convert(bac_iupac, Format::Auto, Format::IupacCondensed)
                .unwrap()
                .ends_with("Bac2Ac4Ac(b1-")
        );
    }

    #[test]
    fn conversion_computes_and_the_table_is_an_explicit_fallback_or_check() {
        // GS00728's WURCS used to be answered from its table row: `G60371D-N`.
        let g60371dn = table_lookup("G60371D-N", Format::IupacCondensed, Format::Wurcs).unwrap();
        let iupac = convert(&g60371dn, Format::Wurcs, Format::IupacCondensed).unwrap();
        assert!(
            iupac.starts_with("Neu5Ac(a2-3)[GalNAc(b1-4)]Gal(b1-4)"),
            "{iupac}"
        );
        assert!(!iupac.contains("G60371"), "{iupac}");
        assert_eq!(
            convert(&iupac, Format::IupacCondensed, Format::Wurcs).unwrap(),
            g60371dn
        );
        // GS00430's WURCS used to give its damaged IUPAC string.
        let g85731cn = "WURCS=2.0/3,3,2/[uxxxxm_2*NCC/3=O_4*NCC/3=O][a2112h-1a_1-5][a2112h-1b_1-5]/1-2-3/a3-b1_b4-c1";
        assert_eq!(
            convert(g85731cn, Format::Wurcs, Format::IupacCondensed).unwrap(),
            "Gal(b1-4)Gal(a1-3)dHexNAc4NAc"
        );
        // GS00954's WURCS used to give its GLYCAM name as "IUPAC".
        let g66397db = "WURCS=2.0/3,7,6/[a2122m-1b_1-5_2*NCC/3=O_4*NCC/3=O][a2112h-1a_1-5_2*NCC/3=O][a2122h-1b_1-5]/1-2-2-2-3-2-2/a3-b1_b4-c1_c4-d1_d3-e1_d4-f1_f4-g1";
        assert_eq!(
            convert(g66397db, Format::Wurcs, Format::IupacCondensed).unwrap(),
            "Glc(b1-3)[GalNAc(a1-4)GalNAc(a1-4)]GalNAc(a1-4)GalNAc(a1-4)GalNAc(a1-3)Bac2Ac4Ac(b1-"
        );

        // Accessions are not converted, only found in the table on request.
        for (accession, from) in [
            ("G60371D-N", Format::Auto),
            ("G60371DN", Format::Glycam),
            ("G50508SG", Format::Auto),
        ] {
            assert!(
                convert(accession, from, Format::Wurcs).is_err(),
                "{accession}"
            );
            assert!(
                convert_with_table(accession, from, Format::Wurcs, TableMode::Off).is_err(),
                "{accession}"
            );
            let resolved =
                convert_with_table(accession, from, Format::Wurcs, TableMode::Fallback).unwrap();
            assert!(resolved.starts_with("WURCS=2.0/8,"), "{accession}");
        }
        // G77147GI is not a table record: the structure 0.3.1 gave for it
        // has four open linkages (that of G92358QN), not that of G77147GI.
        assert!(
            convert_with_table(
                "G77147GI",
                Format::Glycam,
                Format::Wurcs,
                TableMode::Fallback
            )
            .is_err()
        );

        // A check passes when the input reads as its record and fails when
        // the record holds another structure.
        let name = "LFucpa1-2DGalpa1-3DGalpb1-3DGalpNAc";
        assert_eq!(
            convert_with_table(name, Format::Glycam, Format::Wurcs, TableMode::Check).unwrap(),
            WURCS
        );
        let gs00954 = TABLE_RECORDS_THAT_DISAGREE[4].1;
        assert!(matches!(
            convert_with_table(gs00954, Format::Glycam, Format::Wurcs, TableMode::Check),
            Err(TableConversionError::Disagreement { .. })
        ));
        assert!(convert_with_table(gs00954, Format::Glycam, Format::Wurcs, TableMode::Off).is_ok());
    }

    #[test]
    fn every_corpus_notation_is_computed_to_its_record_structure() {
        let records: Vec<_> = corpus_records().collect();
        assert_eq!(records.len(), 942);
        let mut disagreements = Vec::new();
        for record in records {
            let expected = &corpus_wurcs_normalizations()[record.wurcs];
            for (notation, format, column) in [
                (record.iupac, Format::IupacCondensed, "iupac"),
                (
                    record.iupac_extended,
                    Format::IupacExtended,
                    "iupac_extended",
                ),
                (record.glycam, Format::Glycam, "glycam"),
                (record.smiles, Format::Smiles, "smiles"),
            ] {
                if notation.is_empty() {
                    continue;
                }
                let computed = parse_notation(notation, format)
                    .map_err(|error| error.to_string())
                    .and_then(|graph| core::write_wurcs(&graph).map_err(|e| e.to_string()));
                if computed.as_ref() != Ok(expected) {
                    disagreements.push((column, notation, computed, expected.clone()));
                }
            }
        }
        for (column, notation, computed, expected) in &disagreements {
            eprintln!("{column}\t{notation}\n\tcomputed {computed:?}\n\trecord   {expected}");
        }
        let found = disagreements
            .iter()
            .map(|(column, notation, _, _)| (*column, *notation))
            .collect::<Vec<_>>();
        assert_eq!(found, TABLE_RECORDS_THAT_DISAGREE);
    }

    const HIGHLIGHT_TARGET: &str = "Neu5Ac(a2-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-2)Man(a1-3)[Gal(b1-3)[Fuc(a1-4)]GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)[Fuc(a1-6)]GlcNAc";
    const HIGHLIGHT_MOTIF: &str = "Fuc(a1-?)[Gal(b1-?)]GlcNAc";

    #[test]
    fn glycodraw_example_finds_both_non_induced_motif_occurrences() {
        let target = parse_notation(HIGHLIGHT_TARGET, Format::IupacCondensed).unwrap();
        let motif = parse_notation(HIGHLIGHT_MOTIF, Format::IupacCondensed).unwrap();
        let matches = find_motif_matches(&target, &motif).unwrap();
        assert_eq!(matches.len(), 2);
        assert_eq!(
            matches
                .iter()
                .flat_map(|found| found.node_indices.iter())
                .copied()
                .collect::<std::collections::BTreeSet<_>>()
                .len(),
            6
        );
        assert_eq!(
            matches
                .iter()
                .flat_map(|found| found.edge_indices.iter())
                .copied()
                .collect::<std::collections::BTreeSet<_>>()
                .len(),
            4
        );
    }

    #[test]
    fn motif_wildcards_generic_classes_and_exact_positions_are_respected() {
        let target = parse_notation("Gal(b1-4)GlcNAc", Format::IupacCondensed).unwrap();
        for motif in ["Gal(b1-?)HexNAc", "Gal(?1-4)GlcNAc"] {
            let motif = parse_notation(motif, Format::IupacCondensed).unwrap();
            assert_eq!(find_motif_matches(&target, &motif).unwrap().len(), 1);
        }
        for motif in ["Gal(b1-3)GlcNAc", "Gal(a1-4)GlcNAc", "Gal(b1-4)HexN"] {
            let motif = parse_notation(motif, Format::IupacCondensed).unwrap();
            assert!(find_motif_matches(&target, &motif).unwrap().is_empty());
        }

        let uncertain_target = parse_notation("Gal(?1-?)GlcNAc", Format::IupacCondensed).unwrap();
        let exact = parse_notation("Gal(b1-4)GlcNAc", Format::IupacCondensed).unwrap();
        assert!(
            find_motif_matches(&uncertain_target, &exact)
                .unwrap()
                .is_empty()
        );
    }

    #[test]
    fn motif_rendering_uses_glycodraw_palette_outside_the_match_union() {
        let target = parse_notation(HIGHLIGHT_TARGET, Format::IupacCondensed).unwrap();
        let motif = parse_notation(HIGHLIGHT_MOTIF, Format::IupacCondensed).unwrap();
        let ordinary = snfg::render_svg(&target).unwrap();
        assert_eq!(
            ordinary,
            snfg::render_svg_with_motifs(&target, &[], &snfg::RenderOptions::default()).unwrap()
        );

        let highlighted =
            snfg::render_svg_with_motifs(&target, &[motif], &snfg::RenderOptions::default())
                .unwrap();
        assert!(highlighted.contains(r##"[fill="#00A651"] { fill: #CDE9DF; }"##));
        assert!(highlighted.contains(".motif-dimmed text { fill: #D9D9D9; }"));
        assert!(!highlighted.contains("filter:"));
        assert!(!highlighted.contains(".motif-dimmed { opacity:"));
        assert_eq!(highlighted.matches("class=\"motif-match\"").count(), 10);
        assert_eq!(highlighted.matches("class=\"motif-dimmed\"").count(), 15);

        let neu5ac = parse_notation("Neu5Ac", Format::IupacCondensed).unwrap();
        let union = snfg::render_svg_with_motifs(
            &target,
            &[
                parse_notation(HIGHLIGHT_MOTIF, Format::IupacCondensed).unwrap(),
                parse_notation(HIGHLIGHT_MOTIF, Format::IupacCondensed).unwrap(),
                neu5ac,
            ],
            &snfg::RenderOptions::default(),
        )
        .unwrap();
        assert_eq!(union.matches("class=\"motif-match\"").count(), 11);
        assert_eq!(union.matches("class=\"motif-dimmed\"").count(), 14);

        let absent = parse_notation("Kdo", Format::IupacCondensed).unwrap();
        let no_match =
            snfg::render_svg_with_motifs(&target, &[absent], &snfg::RenderOptions::default())
                .unwrap();
        assert_eq!(no_match.matches("class=\"motif-match\"").count(), 0);
        assert_eq!(no_match.matches("class=\"motif-dimmed\"").count(), 25);
    }
}
