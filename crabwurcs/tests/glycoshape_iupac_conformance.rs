//! Offline conformance checks against the IUPAC notations bundled from
//! GlycoShape. The WURCS column is treated as the reference structure; both
//! sides are normalized before they are compared so source traversal order
//! does not affect the result.

use crabwurcs::{core, iupac};

#[derive(Debug, Clone, Copy)]
struct Record<'a> {
    wurcs: &'a str,
    condensed: &'a str,
    extended: &'a str,
}

fn records() -> impl Iterator<Item = Record<'static>> {
    const DATA: &str = include_str!("../data/glycoshape_notations.tsv");
    const DERIVED: &str = include_str!("../data/glycoshape_derived_notations.tsv");
    const NOTATION_ONLY: &str = include_str!("../data/glycoshape_notation_only.tsv");

    DATA.lines()
        .chain(DERIVED.lines())
        .chain(NOTATION_ONLY.lines())
        .filter_map(|line| {
            let mut fields = line.splitn(5, '\t');
            Some(Record {
                wurcs: fields.next()?,
                condensed: fields.next()?,
                extended: fields.next()?,
            })
        })
}

/// Strings stored in the IUPAC column that are not IUPAC notations: a
/// GlyTouCan accession (GlycoShape GS00728) and a damaged string (GS00430).
/// Only the table fallback resolves them.
const NOT_IUPAC: &[&str] = &[
    "G60371D-N",
    "Galb1-4)Gala1-3(?24-diacetimido-246-trideoxyhexose)",
];

#[test]
fn glycoshape_iupac_notations_match_normalized_wurcs() {
    let mut checked = 0usize;
    let mut skipped = Vec::new();
    let mut mismatches = Vec::new();

    for record in records() {
        if record.wurcs.is_empty() || NOT_IUPAC.contains(&record.condensed) {
            continue;
        }

        let expected = match core::normalize_wurcs(record.wurcs) {
            Ok(value) => value,
            Err(error) => {
                skipped.push(format!("reference WURCS: {} ({error})", record.wurcs));
                continue;
            }
        };
        for (kind, notation, parse) in [
            (
                "condensed",
                record.condensed,
                iupac::parse_iupac_condensed as fn(&str) -> Result<_, _>,
            ),
            (
                "extended",
                record.extended,
                iupac::parse_iupac_extended as fn(&str) -> Result<_, _>,
            ),
        ] {
            if notation.is_empty() {
                continue;
            }
            let graph = match parse(notation) {
                Ok(graph) => graph,
                Err(error) => {
                    skipped.push(format!("{kind} `{notation}` ({error})"));
                    continue;
                }
            };
            let actual = match core::write_wurcs(&graph) {
                Ok(value) => value,
                Err(error) => {
                    skipped.push(format!("write {kind} `{notation}` ({error})"));
                    continue;
                }
            };
            checked += 1;
            if actual != expected {
                mismatches.push((kind, notation, expected.clone(), actual));
            }
        }
    }

    println!(
        "GlycoShape IUPAC conformance: checked={}, skipped={}, mismatches={}",
        checked,
        skipped.len(),
        mismatches.len()
    );
    if !skipped.is_empty() {
        println!("Skipped unsupported records (first 20):");
        for entry in skipped.iter().take(20) {
            println!("  {entry}");
        }
    }
    if !mismatches.is_empty() {
        println!("Mismatches (first 20):");
        for (kind, notation, expected, actual) in mismatches.iter().take(20) {
            println!("  {kind}: {notation}\n    expected: {expected}\n    actual:   {actual}");
        }
    }

    assert!(checked > 0, "no GlycoShape IUPAC records were checked");
    assert!(
        skipped.is_empty(),
        "{} GlycoShape IUPAC records could not be parsed or written",
        skipped.len()
    );
    assert!(
        mismatches.is_empty(),
        "{} GlycoShape IUPAC records produced a different normalized WURCS",
        mismatches.len()
    );
}
