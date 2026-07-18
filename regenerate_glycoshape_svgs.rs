// Regenerate all GlycoShape SVG files from WURCS data in GLYCOSHAPE.json

use crabwurcs_core::parse_wurcs;
use crabwurcs_snfg::render_svg;
use serde::Deserialize;
use std::collections::HashMap;
use std::fs::{self, File};
use std::io::Write;

#[derive(Debug, Deserialize)]
struct Archetype {
    #[serde(rename = "ID")]
    id: String,
    #[allow(dead_code)]
    name: Option<String>,
    wurcs: Option<String>,
    #[allow(dead_code)]
    iupac: Option<String>,
    #[allow(dead_code)]
    iupac_extended: Option<String>,
    #[allow(dead_code)]
    glytoucan: Option<String>,
}

#[derive(Debug, Deserialize)]
struct Entry {
    archetype: Archetype,
}

fn load_entries() -> Vec<(String, Archetype)> {
    let path = "GLYCOSHAPE.json";
    let data = std::fs::read_to_string(path).expect("Cannot read GLYCOSHAPE.json");
    let raw: HashMap<String, Entry> =
        serde_json::from_str(&data).expect("Cannot parse GLYCOSHAPE.json");
    raw.into_iter().map(|(k, v)| (k, v.archetype)).collect()
}

fn main() {
    let entries = load_entries();
    println!("Loaded {} entries from GLYCOSHAPE.json", entries.len());

    let mut success_count = 0;
    let mut failure_count = 0;
    let mut skipped_count = 0;

    for (key, archetype) in &entries {
        let wurcs = match &archetype.wurcs {
            Some(w) => w,
            None => {
                skipped_count += 1;
                continue;
            }
        };

        // Parse WURCS
        let graph = match parse_wurcs(wurcs) {
            Ok(g) => g,
            Err(e) => {
                eprintln!("Failed to parse WURCS for {}: {:?}", key, e);
                failure_count += 1;
                continue;
            }
        };

        // Render SNFG
        let svg = match render_svg(&graph) {
            Ok(s) => s,
            Err(e) => {
                eprintln!("Failed to render SNFG for {}: {:?}", key, e);
                failure_count += 1;
                continue;
            }
        };

        // Write SVG file
        let filename = format!("{}.snfg.svg", key);
        match File::create(&filename) {
            Ok(mut file) => {
                if let Err(e) = file.write_all(svg.as_bytes()) {
                    eprintln!("Failed to write SVG file {}: {:?}", filename, e);
                    failure_count += 1;
                } else {
                    success_count += 1;
                }
            }
            Err(e) => {
                eprintln!("Failed to create SVG file {}: {:?}", filename, e);
                failure_count += 1;
            }
        }
    }

    println!("\n--- Summary ---");
    println!("Total entries: {}", entries.len());
    println!("Successfully generated: {} SVG files", success_count);
    println!("Failed: {}", failure_count);
    println!("Skipped (no WURCS): {}", skipped_count);
}
