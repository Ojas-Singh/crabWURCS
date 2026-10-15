//! Opt-in checks of the coordinate reading on the GlycoShape database.
//!
//! ```text
//! CRABWURCS_GLYCOSHAPE_STATIC=/path/to/static \
//! CRABWURCS_ENTRY_AUDIT=/path/to/entry-audit \
//! CRABWURCS_STEREO_CHECK=/path/to/stereo-check \
//! CRABWURCS_GLYCOSHAPE_ARCHIVE=/path/to/archive \
//! CRABWURCS_STEREO_LEVEL1_RINGS=/path/to/level1_rings.jsonl \
//!   cargo test -p crabwurcs-pdb --test glycoshape_coordinates -- --ignored --nocapture
//! ```
//!
//! * `CRABWURCS_GLYCOSHAPE_STATIC`: a static build (`GLYCOSHAPE.json`,
//!   `GSxxxxx/output/level_1/PDB/{alpha,beta}.pdb`).
//! * `CRABWURCS_ENTRY_AUDIT`: GlycanAnalysisPipeline's entry audit
//!   (`GSxxxxx/issues.json`).
//! * `CRABWURCS_STEREO_CHECK`: an independent, name-blind reading of ring
//!   configurations (`results/GSxxxxx.json`, one record per ring with the
//!   signs of its signed volumes), here compared on the raw topologies of
//!   `CRABWURCS_GLYCOSHAPE_ARCHIVE` (`<archive_folder>/topology.pdb`).
//! * `CRABWURCS_STEREO_LEVEL1_RINGS` (optional): the same independent
//!   reader run on every Level 1 file, one JSON record per ring.

use crabwurcs_core::write_wurcs;
use crabwurcs_pdb::{
    NamePolicy, PdbReadOptions, PdbReading, PdbWarningKind, ResidueReading, read_glycans_from_file,
};
use serde_json::Value;
use std::collections::{BTreeMap, HashMap, HashSet};
use std::path::{Path, PathBuf};

fn directory(variable: &str) -> PathBuf {
    PathBuf::from(std::env::var(variable).unwrap_or_else(|_| panic!("set {variable}")))
}

fn json(path: &Path) -> Value {
    serde_json::from_str(
        &std::fs::read_to_string(path)
            .unwrap_or_else(|error| panic!("{}: {error}", path.display())),
    )
    .unwrap()
}

/// Issue codes of every audited entry.
fn audit_issues(audit: &Path) -> HashMap<String, HashSet<String>> {
    let mut issues = HashMap::new();
    for entry in std::fs::read_dir(audit).unwrap().flatten() {
        let path = entry.path().join("issues.json");
        if !path.exists() {
            continue;
        }
        let value = json(&path);
        let id = value["glycoshape_id"].as_str().unwrap().to_string();
        let codes = value["issues"]
            .as_array()
            .unwrap()
            .iter()
            .filter_map(|issue| issue["code"].as_str().map(str::to_string))
            .collect();
        issues.insert(id, codes);
    }
    issues
}

fn read(path: &Path) -> PdbReading {
    read_glycans_from_file(path, &PdbReadOptions::with_policy(NamePolicy::Coordinates))
        .unwrap_or_else(|error| panic!("{}: {error}", path.display()))
}

#[test]
#[ignore = "needs a GlycoShape static build and its entry audit (see the module documentation)"]
fn level1_files_read_as_their_names_except_where_the_audit_finds_otherwise() {
    let root = directory("CRABWURCS_GLYCOSHAPE_STATIC");
    let issues = audit_issues(&directory("CRABWURCS_ENTRY_AUDIT"));
    let entries = json(&root.join("GLYCOSHAPE.json"));
    let mut files = 0;
    let mut matching = 0;
    let mut name_unreadable = 0;
    let mut explained: BTreeMap<&str, Vec<String>> = BTreeMap::new();
    let mut unexplained = Vec::new();
    let mut contradiction_files = 0;
    let mut contradiction_warnings = 0;
    let mut other_warnings = 0;
    for (id, entry) in entries.as_object().unwrap() {
        for form in ["alpha", "beta"] {
            let Some(name) = entry[form]["glycam"].as_str() else {
                continue;
            };
            let path = root
                .join(id)
                .join("output/level_1/PDB")
                .join(format!("{form}.pdb"));
            if !path.exists() {
                continue;
            }
            files += 1;
            let reading = read(&path);
            let contradictions = reading
                .warnings
                .iter()
                .filter(|warning| warning.kind == PdbWarningKind::NameContradictsCoordinates)
                .count();
            contradiction_warnings += contradictions;
            other_warnings += reading.warnings.len() - contradictions;
            if contradictions > 0 {
                contradiction_files += 1;
            }
            let Ok(expected) = crabwurcs_iupac::parse_glycam(name)
                .and_then(|graph| write_wurcs(&graph).map_err(crabwurcs_iupac::IupacError::from))
            else {
                name_unreadable += 1;
                continue;
            };
            let read_back = reading
                .glycans
                .iter()
                .map(|glycan| write_wurcs(&glycan.graph).unwrap())
                .collect::<Vec<_>>();
            if read_back == [expected.clone()] {
                matching += 1;
                continue;
            }
            let codes = issues.get(id).cloned().unwrap_or_default();
            let reason = [
                "structure-contradicts-name",
                "level1-file-anomer",
                "mass-mismatch",
            ]
            .into_iter()
            .find(|code| codes.contains(*code));
            match reason {
                Some(reason) => explained
                    .entry(reason)
                    .or_default()
                    .push(format!("{id} {form}")),
                None => unexplained.push(format!("{id} {form}: {read_back:?} (name {expected})")),
            }
        }
    }
    println!("Level 1 files read: {files}");
    println!("  read as the entry name: {matching}");
    println!("  entry name not a GLYCAM sequence: {name_unreadable}");
    for (reason, list) in &explained {
        println!(
            "  differ, audit issue {reason}: {} ({})",
            list.len(),
            list.join(", ")
        );
    }
    println!("  differ, unexplained: {}", unexplained.len());
    println!(
        "  files with a name-contradicts-coordinates warning: {contradiction_files} ({contradiction_warnings} warnings); other warnings: {other_warnings}"
    );
    assert!(unexplained.is_empty(), "{unexplained:#?}");
}

/// The Fischer side (`'1'` left, `'2'` right) that a centre of the
/// independent reader's signature stands for, given its label and the ring
/// size. That reader orders the three substituents of the ring-closing
/// carbon (`R5` of a pyranose, `R4` of a furanose) and of the exocyclic
/// centre `X` the other way round from this crate.
fn fischer_side(label: &str, sign: char, ring_size: usize) -> char {
    let flipped = label == "X" || label == format!("R{}", ring_size - 1);
    match (sign, flipped) {
        ('+', false) | ('-', true) => '1',
        _ => '2',
    }
}

/// One ring of the independent reader.
struct Ring<'a> {
    residue_name: &'a str,
    residue_number: isize,
    size: usize,
    centres: Vec<String>,
    signature: &'a str,
}

#[derive(Default)]
struct Agreement {
    rings: usize,
    rings_agreeing: usize,
    centres: usize,
    centres_agreeing: usize,
    missing: Vec<String>,
    disagreeing: Vec<String>,
}

impl Agreement {
    fn compare(&mut self, context: &str, residues: &[ResidueReading], ring: &Ring) {
        self.rings += 1;
        let (residue_name, residue_number) = (ring.residue_name, ring.residue_number);
        let Some(residue) = residues.iter().find(|residue| {
            residue.residue_name == residue_name && residue.sequence_number == residue_number
        }) else {
            self.missing
                .push(format!("{context} {residue_name}{residue_number}"));
            return;
        };
        let mut agree = true;
        for (atom, item) in ring.centres.iter().zip(ring.signature.split_whitespace()) {
            let (label, sign) = item.split_at(item.len() - 1);
            let label = label.trim_end_matches('k');
            let expected = fischer_side(label, sign.chars().next().unwrap(), ring.size);
            // The ketose anomeric centre (`R1k`) is ordered (ring O, C3, O2)
            // by both readers.
            self.centres += 1;
            match residue.centres.iter().find(|centre| &centre.atom == atom) {
                Some(centre) if centre.code == expected => self.centres_agreeing += 1,
                found => {
                    agree = false;
                    self.disagreeing.push(format!(
                        "{context} {residue_name}{residue_number} {atom}: independent {item}, crabWURCS {:?}",
                        found.map(|centre| (centre.code, centre.volume))
                    ));
                }
            }
        }
        if agree {
            self.rings_agreeing += 1;
        }
    }

    fn report(&self, what: &str) {
        println!(
            "{what}: {} of {} rings and {} of {} centres agree; {} rings not found",
            self.rings_agreeing,
            self.rings,
            self.centres_agreeing,
            self.centres,
            self.missing.len()
        );
    }
}

#[test]
#[ignore = "needs the independent stereo check and the GlycoShape archive (see the module documentation)"]
fn coordinates_agree_with_an_independent_reader_ring_by_ring() {
    let stereo = directory("CRABWURCS_STEREO_CHECK");
    let archive = directory("CRABWURCS_GLYCOSHAPE_ARCHIVE");
    let mut agreement = Agreement::default();
    let mut entries = 0;
    let mut paths = std::fs::read_dir(stereo.join("results"))
        .unwrap()
        .flatten()
        .map(|entry| entry.path())
        .collect::<Vec<_>>();
    paths.sort();
    for path in paths {
        let result = json(&path);
        let id = result["glycoshape_id"].as_str().unwrap();
        // Entries found to hold the other reducing-end anomer were renamed
        // in the archive after the check (`...a1-OH` to `...b1-OH`).
        let folder = result["archive_folder"].as_str().unwrap();
        let renamed = folder
            .strip_suffix("a1-OH")
            .map(|stem| format!("{stem}b1-OH"))
            .or_else(|| {
                folder
                    .strip_suffix("b1-OH")
                    .map(|stem| format!("{stem}a1-OH"))
            });
        let Some(topology) = std::iter::once(folder.to_string())
            .chain(renamed)
            .map(|folder| archive.join(folder).join("topology.pdb"))
            .find(|path| path.exists())
        else {
            continue;
        };
        entries += 1;
        let reading = read(&topology);
        for ring in result["rings"].as_array().unwrap() {
            // `QYS2`, `0fA7`: a three-character residue name, then its number.
            let (name, number) = ring["residue"].as_str().unwrap().split_at(3);
            let centres = ring["centre_atom_names"]
                .as_array()
                .unwrap()
                .iter()
                .map(|value| value.as_str().unwrap().to_string())
                .collect::<Vec<_>>();
            agreement.compare(
                id,
                &reading.residues,
                &Ring {
                    residue_name: name,
                    residue_number: number.parse().unwrap(),
                    size: ring["ring_atom_names"].as_array().unwrap().len(),
                    centres,
                    signature: ring["signature_topology"].as_str().unwrap(),
                },
            );
        }
    }
    agreement.report(&format!("raw topologies of {entries} entries"));
    assert!(entries > 0);
    assert!(agreement.missing.is_empty(), "{:#?}", agreement.missing);
    assert!(
        agreement.disagreeing.is_empty(),
        "{:#?}",
        agreement.disagreeing
    );

    let Ok(level1) = std::env::var("CRABWURCS_STEREO_LEVEL1_RINGS") else {
        return;
    };
    let mut agreement = Agreement::default();
    let mut readings: HashMap<String, PdbReading> = HashMap::new();
    for line in std::fs::read_to_string(level1).unwrap().lines() {
        let ring: Value = serde_json::from_str(line).unwrap();
        let file = ring["file"].as_str().unwrap().to_string();
        let reading = readings
            .entry(file.clone())
            .or_insert_with(|| read(Path::new(&file)));
        let centres = ring["centres"]
            .as_array()
            .unwrap()
            .iter()
            .map(|value| value.as_str().unwrap().to_string())
            .collect::<Vec<_>>();
        agreement.compare(
            &file,
            &reading.residues,
            &Ring {
                residue_name: ring["residue_name"].as_str().unwrap(),
                residue_number: ring["residue_number"].as_i64().unwrap() as isize,
                size: ring["ring_size"].as_u64().unwrap() as usize,
                centres,
                signature: ring["signature"].as_str().unwrap(),
            },
        );
    }
    agreement.report(&format!("Level 1 files ({})", readings.len()));
    assert!(agreement.missing.is_empty(), "{:#?}", agreement.missing);
    assert!(
        agreement.disagreeing.is_empty(),
        "{:#?}",
        agreement.disagreeing
    );
}
