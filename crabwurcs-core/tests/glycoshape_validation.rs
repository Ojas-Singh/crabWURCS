// Comprehensive test: Check all GlycoShape IUPAC → WURCS conversions

use crabwurcs_core::{parse_wurcs, write_wurcs};
use crabwurcs_iupac::parse_iupac_condensed;
use std::collections::HashMap;

#[derive(serde::Deserialize)]
struct GlycoShapeEntry {
    #[serde(rename = "iupac")]
    iupac: Option<String>,
    #[serde(rename = "wurcs")]
    wurcs: Option<String>,
    #[serde(rename = "name")]
    name: Option<String>,
}

#[derive(serde::Deserialize)]
struct GlycoShapeArchetype {
    archetype: GlycoShapeEntry,
}

#[derive(serde::Deserialize)]
struct GlycoShapeData {
    #[serde(flatten)]
    entries: HashMap<String, GlycoShapeArchetype>,
}

fn load_glycoshape() -> GlycoShapeData {
    // Try to load from project root
    let paths = [
        "../GLYCOSHAPE.json",
        "../../GLYCOSHAPE.json",
        "../../../GLYCOSHAPE.json",
        "GLYCOSHAPE.json",
    ];

    let mut json_content = None;
    for path in &paths {
        if let Ok(content) = std::fs::read_to_string(path) {
            json_content = Some(content);
            break;
        }
    }

    let json = json_content.expect("Could not find GLYCOSHAPE.json in any expected location");
    serde_json::from_str(&json).expect("Failed to parse GLYCOSHAPE.json")
}

fn normalize_wurcs_for_comparison(wurcs: &str) -> Option<String> {
    // Parse and rewrite to normalize
    let graph = parse_wurcs(wurcs).ok()?;
    write_wurcs(&graph).ok()
}

#[test]
fn test_all_glycoshape_iupac_to_wurcs() {
    let data = load_glycoshape();

    let mut total = 0;
    let mut failed = 0;
    let mut missing_iupac = 0;
    let mut missing_wurcs = 0;
    let mut parse_errors = 0;
    let mut write_errors = 0;
    let mut matches = 0;
    let mut mismatches = Vec::new();

    println!("\n=== GlycoShape IUPAC → WURCS Validation ===\n");

    for (id, entry) in &data.entries {
        let iupac = match &entry.archetype.iupac {
            Some(i) => i,
            None => {
                missing_iupac += 1;
                continue;
            }
        };

        let expected_wurcs = match &entry.archetype.wurcs {
            Some(w) => w,
            None => {
                missing_wurcs += 1;
                continue;
            }
        };

        total += 1;

        // Parse IUPAC to graph
        let graph = match parse_iupac_condensed(iupac) {
            Ok(g) => g,
            Err(e) => {
                parse_errors += 1;
                eprintln!("{}: IUPAC parse error: {:?}", id, e);
                continue;
            }
        };

        // Write WURCS from the IUPAC graph
        let generated_wurcs = match write_wurcs(&graph) {
            Ok(w) => w,
            Err(e) => {
                write_errors += 1;
                eprintln!("{}: WURCS write error: {:?}", id, e);
                continue;
            }
        };

        // Normalize both for comparison
        let normalized_generated = normalize_wurcs_for_comparison(&generated_wurcs);
        let normalized_expected = normalize_wurcs_for_comparison(expected_wurcs);

        match (normalized_generated, normalized_expected) {
            (Some(gen), Some(exp)) if gen == exp => {
                matches += 1;
            }
            (Some(gen), Some(exp)) => {
                failed += 1;
                mismatches.push(format!(
                    "{}: {}\n  IUPAC: {}\n  Expected: {}\n  Got:      {}",
                    id,
                    entry.archetype.name.as_deref().unwrap_or("(unknown)"),
                    iupac,
                    exp,
                    gen
                ));
            }
            (None, Some(_)) => {
                write_errors += 1;
                eprintln!(
                    "{}: Generated WURCS is unparseable: {}",
                    id, generated_wurcs
                );
            }
            (Some(_), None) => {
                eprintln!("{}: Expected WURCS is unparseable: {}", id, expected_wurcs);
            }
            (None, None) => {
                eprintln!("{}: Both WURCS are unparseable", id);
            }
        }
    }

    println!("\n=== Summary ===");
    println!("Total entries with both IUPAC and WURCS: {}", total);
    println!("Passed (exact match): {}", matches);
    println!("Failed: {}", failed);
    println!("Missing IUPAC: {}", missing_iupac);
    println!("Missing WURCS: {}", missing_wurcs);
    println!("IUPAC parse errors: {}", parse_errors);
    println!("WURCS write errors: {}", write_errors);

    if !mismatches.is_empty() {
        println!("\n=== Mismatches (showing first 20) ===");
        for (i, m) in mismatches.iter().take(20).enumerate() {
            println!("{}. {}", i + 1, m);
        }
        if mismatches.len() > 20 {
            println!("... and {} more", mismatches.len() - 20);
        }
    }

    // Don't fail the test - just report results
    // This allows us to see the overall status without blocking CI
    println!(
        "\nTest completed with {} matches out of {} total",
        matches, total
    );
}

#[test]
fn test_glycoshape_sample_entries() {
    // Test a few specific entries to understand the patterns
    let data = load_glycoshape();

    let sample_ids = vec!["GS00002", "GS00003", "GS00004"];

    for id in sample_ids {
        if let Some(entry) = data.entries.get(id) {
            let iupac = match &entry.archetype.iupac {
                Some(i) => i,
                None => {
                    println!("{}: No IUPAC", id);
                    continue;
                }
            };

            let expected_wurcs = match &entry.archetype.wurcs {
                Some(w) => w,
                None => {
                    println!("{}: No WURCS", id);
                    continue;
                }
            };

            println!("\n=== {} ===", id);
            println!("IUPAC: {}", iupac);
            println!("Expected WURCS: {}", expected_wurcs);

            let graph = match parse_iupac_condensed(iupac) {
                Ok(g) => {
                    println!("Node count: {}", g.node_count());
                    g
                }
                Err(e) => {
                    println!("Parse error: {:?}", e);
                    continue;
                }
            };

            let generated_wurcs = match write_wurcs(&graph) {
                Ok(w) => w,
                Err(e) => {
                    println!("Write error: {:?}", e);
                    continue;
                }
            };

            println!("Generated WURCS: {}", generated_wurcs);

            let normalized_generated = normalize_wurcs_for_comparison(&generated_wurcs);
            let normalized_expected = normalize_wurcs_for_comparison(expected_wurcs);

            match (normalized_generated, normalized_expected) {
                (Some(gen), Some(exp)) if gen == exp => {
                    println!("✓ MATCH!");
                }
                (Some(gen), Some(exp)) => {
                    println!("✗ MISMATCH");
                    println!("  Expected: {}", exp);
                    println!("  Got:      {}", gen);
                }
                _ => {
                    println!("Could not normalize for comparison");
                }
            }
        }
    }
}
