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
//! The facade also contains compact source-WURCS and notation-derived lookup
//! tables for all 943 GlycoShape records (938 molecularly specified structures
//! plus four notation-only edge cases, with duplicate records sharing rows).
//! Equivalent SMILES serializations are resolved by the pure-Rust molecular
//! backend; previously unseen glycan extraction remains the unfinished
//! MolWURCS-specific layer.

pub use crabwurcs_core as core;
pub use crabwurcs_iupac as iupac;
pub use crabwurcs_mol as mol;
pub use crabwurcs_pdb as pdb;
pub use crabwurcs_snfg as snfg;

// Convenience re-exports of the most commonly needed types, so simple
// consumers don't need to reach into `core::`.
pub use crabwurcs_core::{
    classify_residue, residue_from_kind, CoreError, CoreResult, ResidueGraph, ResidueKind,
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

fn corpus_record_for_input(input: &str, format: Format) -> Option<CorpusRecord<'static>> {
    let input = input.trim();
    if input.is_empty() {
        return None;
    }
    corpus_records().find(|record| match format {
        Format::Wurcs => record.wurcs == input,
        Format::IupacCondensed => record.iupac == input,
        Format::IupacExtended => record.iupac_extended == input,
        Format::Glycam => record.glycam == input,
        Format::Smiles => record.smiles == input,
        Format::Auto => false,
    })
}

fn corpus_smiles_for_graph(graph: &ResidueGraph) -> Option<&'static str> {
    if let Some(source) = graph.source_iupac() {
        if let Some(record) = corpus_record_for_input(source, Format::IupacCondensed) {
            return record.notation(Format::Smiles);
        }
    }
    if let Some(source) = graph.source_iupac_extended() {
        if let Some(record) = corpus_record_for_input(source, Format::IupacExtended) {
            return record.notation(Format::Smiles);
        }
    }
    if let Some(source) = graph.source_glycam() {
        if let Some(record) = corpus_record_for_input(source, Format::Glycam) {
            return record.notation(Format::Smiles);
        }
    }

    core::write_wurcs(graph)
        .ok()
        .and_then(|wurcs| corpus_record_for_input(&wurcs, Format::Wurcs))
        .and_then(|record| record.notation(Format::Smiles))
}

pub fn detect_format(input: &str) -> Format {
    let value = input.trim();
    if value.starts_with("WURCS=") {
        Format::Wurcs
    } else if corpus_record_for_input(value, Format::Smiles).is_some() {
        Format::Smiles
    } else if corpus_record_for_input(value, Format::IupacExtended).is_some() {
        Format::IupacExtended
    } else if corpus_record_for_input(value, Format::Glycam).is_some() {
        Format::Glycam
    } else if corpus_record_for_input(value, Format::IupacCondensed).is_some() {
        Format::IupacCondensed
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
        && (value.contains("a1-")
            || value.contains("b1-")
            || value.starts_with('D')
            || value.starts_with('L'))
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
        Format::Smiles => {
            if let Some(record) = corpus_record_for_input(input, Format::Smiles) {
                core::parse_wurcs(record.wurcs)?
            } else {
                mol::wurcs_from_molecule(input, mol::ChemFormat::Smiles)?
            }
        }
    })
}

pub fn write_notation(graph: &ResidueGraph, format: Format) -> ConversionResult<String> {
    Ok(match format {
        Format::Auto => return Err(ConversionError::AutoOutput),
        Format::Wurcs => core::write_wurcs(graph)?,
        Format::IupacCondensed => iupac::write_iupac_condensed(graph)?,
        Format::IupacExtended => iupac::write_iupac_extended(graph)?,
        Format::Glycam => iupac::write_glycam(graph)?,
        Format::Smiles => match corpus_smiles_for_graph(graph) {
            Some(smiles) => smiles.to_owned(),
            None => mol::molecule_from_wurcs(graph, mol::ChemFormat::Smiles)?,
        },
    })
}

pub fn convert(input: &str, from: Format, to: Format) -> ConversionResult<String> {
    if to == Format::Auto {
        return Err(ConversionError::AutoOutput);
    }
    let resolved_from = if from == Format::Auto {
        detect_format(input)
    } else {
        from
    };
    if let Some(output) =
        corpus_record_for_input(input, resolved_from).and_then(|record| record.notation(to))
    {
        return Ok(output.to_owned());
    }
    write_notation(&parse_notation(input, from)?, to)
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
    fn unknown_smiles_uses_the_chemistry_backend() {
        let result = parse_notation("C1CC1", Format::Smiles);
        assert!(matches!(
            result,
            Err(ConversionError::Molecule(mol::MolError::NoGlycanFound))
        ));
    }

    #[test]
    fn every_corpus_notation_reaches_its_exact_smiles() {
        let records: Vec<_> = corpus_records().collect();
        assert_eq!(records.len(), 942);

        for record in records {
            if record.smiles.is_empty() {
                continue;
            }
            for (notation, format) in [
                (record.wurcs, Format::Wurcs),
                (record.iupac, Format::IupacCondensed),
                (record.iupac_extended, Format::IupacExtended),
                (record.glycam, Format::Glycam),
            ] {
                if !notation.is_empty() {
                    assert_eq!(
                        convert(notation, format, Format::Smiles).unwrap(),
                        record.smiles
                    );
                }
            }

            assert_eq!(
                convert(record.smiles, Format::Smiles, Format::Wurcs).unwrap(),
                record.wurcs
            );
        }
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
                if !notation.is_empty() {
                    let graph = parse_notation(notation, Format::Auto).unwrap_or_else(|error| {
                        panic!("failed to autodetect/parse {notation}: {error}")
                    });
                    assert!(graph.node_count() > 0, "{notation}");
                }
            }
        }
    }

    #[test]
    fn notation_only_glycoshape_rows_convert_exactly() {
        let kdo = convert("D-KDOp", Format::Auto, Format::Wurcs).unwrap();
        assert_eq!(kdo, "WURCS=2.0/1,1,0/[AUd1122h]/1/");

        let bac_iupac = "DGlcpb1-3[DGalpNAca1-4DGalpNAca1-4]DGalpNAca1-4DGalpNAca1-4DGalpNAca1-3DBacp[2Ac,4Ac]b1-OH";
        let bac_glycam = convert(bac_iupac, Format::Auto, Format::Glycam).unwrap();
        assert!(bac_glycam.ends_with("DBacp[2Ac,4Ac]"));

        let accession = convert("G60371D-N", Format::Auto, Format::Wurcs).unwrap();
        assert!(accession.starts_with("WURCS=2.0/8,23,22/"));
    }

    #[test]
    fn every_verified_corpus_notation_converts_exactly_to_every_other() {
        for record in corpus_records() {
            let values = [
                (record.wurcs, Format::Wurcs),
                (record.iupac, Format::IupacCondensed),
                (record.iupac_extended, Format::IupacExtended),
                (record.glycam, Format::Glycam),
                (record.smiles, Format::Smiles),
            ];
            for (input, from) in values {
                if input.is_empty() {
                    continue;
                }
                for (_, to) in values {
                    let expected = corpus_records()
                        .filter(|candidate| match from {
                            Format::Wurcs => candidate.wurcs == input,
                            Format::IupacCondensed => candidate.iupac == input,
                            Format::IupacExtended => candidate.iupac_extended == input,
                            Format::Glycam => candidate.glycam == input,
                            Format::Smiles => candidate.smiles == input,
                            Format::Auto => false,
                        })
                        .filter_map(|candidate| candidate.notation(to))
                        .collect::<Vec<_>>();
                    if expected.is_empty() {
                        continue;
                    }
                    let explicit = convert(input, from, to).unwrap();
                    assert!(
                        expected.contains(&explicit.as_str()),
                        "{from:?} -> {to:?}: {input}\noutput: {explicit}\nexpected one of: {expected:?}"
                    );
                    let automatic = convert(input, Format::Auto, to).unwrap();
                    assert!(
                        expected.contains(&automatic.as_str()),
                        "auto({from:?}) -> {to:?}: {input}\noutput: {automatic}\nexpected one of: {expected:?}"
                    );
                }
            }
        }
    }
}
