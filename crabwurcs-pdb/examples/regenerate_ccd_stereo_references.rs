//! Regenerate the wwPDB CCD coordinate references of the stereochemistry
//! tests (`tests/data/ccd_stereo_references.tsv`).
//!
//! Usage:
//! `cargo run -p crabwurcs-pdb --example regenerate_ccd_stereo_references -- DIRECTORY DOWNLOADED > crabwurcs-pdb/tests/data/ccd_stereo_references.tsv`
//!
//! `DIRECTORY` holds one wwPDB Chemical Component Dictionary file per
//! component, `<ID>.cif`, as served at
//! `https://files.wwpdb.org/pub/pdb/refdata/chem_comp/<last character>/<ID>/<ID>.cif`;
//! `DOWNLOADED` is the date (or dates) they were downloaded. The command
//! reads local files only. It writes the heavy atoms of each component's
//! ideal (`pdbx_model_Cartn_*_ideal`) and model (`model_Cartn_*`)
//! coordinates, with the component's name and the SHA-256 of its file.

use std::path::Path;
use std::process::Command;

fn sha256(path: &Path) -> Result<String, String> {
    let output = Command::new("shasum")
        .args(["-a", "256"])
        .arg(path)
        .output()
        .map_err(|error| format!("cannot run shasum: {error}"))?;
    if !output.status.success() {
        return Err(format!("shasum failed with {}", output.status));
    }
    String::from_utf8(output.stdout)
        .map_err(|error| error.to_string())?
        .split_whitespace()
        .next()
        .map(str::to_string)
        .ok_or_else(|| "shasum produced no checksum".into())
}

/// Split one CIF data line into values (single or double quotes allowed).
fn values(line: &str) -> Vec<String> {
    let mut values = Vec::new();
    let mut characters = line.chars().peekable();
    while let Some(&character) = characters.peek() {
        if character.is_whitespace() {
            characters.next();
            continue;
        }
        let mut value = String::new();
        if character == '"' || character == '\'' {
            characters.next();
            while let Some(next) = characters.next() {
                if next == character && characters.peek().is_none_or(|after| after.is_whitespace())
                {
                    break;
                }
                value.push(next);
            }
        } else {
            while let Some(&next) = characters.peek() {
                if next.is_whitespace() {
                    break;
                }
                value.push(next);
                characters.next();
            }
        }
        values.push(value);
    }
    values
}

struct Component {
    id: String,
    name: String,
    atoms: Vec<Vec<(String, String)>>,
}

fn read_component(path: &Path) -> Result<Component, String> {
    let text = std::fs::read_to_string(path).map_err(|error| error.to_string())?;
    let mut id = None;
    let mut name = None;
    let lines = text.lines().collect::<Vec<_>>();
    for (index, line) in lines.iter().enumerate() {
        if let Some(rest) = line.strip_prefix("_chem_comp.id ") {
            id = values(rest).into_iter().next();
        }
        if let Some(rest) = line.strip_prefix("_chem_comp.name ") {
            name = values(rest).into_iter().next().or_else(|| {
                // A multi-line value: `;text` on the next line.
                lines
                    .get(index + 1)
                    .and_then(|next| next.strip_prefix(';'))
                    .map(str::to_string)
            });
        }
    }
    let mut atoms = Vec::new();
    let mut index = 0;
    while index < lines.len() {
        if lines[index].trim() == "loop_"
            && lines
                .get(index + 1)
                .is_some_and(|line| line.starts_with("_chem_comp_atom."))
        {
            index += 1;
            let mut keys = Vec::new();
            while let Some(key) = lines[index].strip_prefix("_chem_comp_atom.") {
                keys.push(key.trim().to_string());
                index += 1;
            }
            while index < lines.len()
                && !lines[index].starts_with('#')
                && !lines[index].starts_with('_')
                && !lines[index].starts_with("loop_")
            {
                let row = values(lines[index]);
                if row.len() == keys.len() {
                    atoms.push(keys.iter().cloned().zip(row).collect());
                }
                index += 1;
            }
            break;
        }
        index += 1;
    }
    Ok(Component {
        id: id.ok_or_else(|| format!("{}: no _chem_comp.id", path.display()))?,
        name: name.unwrap_or_default(),
        atoms,
    })
}

fn main() -> Result<(), String> {
    let mut arguments = std::env::args().skip(1);
    let directory = arguments
        .next()
        .ok_or("usage: regenerate_ccd_stereo_references DIRECTORY DOWNLOADED")?;
    let downloaded = arguments
        .next()
        .ok_or("usage: regenerate_ccd_stereo_references DIRECTORY DOWNLOADED")?;
    let mut paths = std::fs::read_dir(&directory)
        .map_err(|error| format!("{directory}: {error}"))?
        .filter_map(|entry| entry.ok().map(|entry| entry.path()))
        .filter(|path| path.extension().is_some_and(|extension| extension == "cif"))
        .collect::<Vec<_>>();
    paths.sort();

    println!(
        "# crabWURCS stereochemistry references: wwPDB Chemical Component Dictionary coordinates"
    );
    println!(
        "# Source: https://files.wwpdb.org/pub/pdb/refdata/chem_comp/<last character>/<ID>/<ID>.cif"
    );
    println!("# Downloaded: {downloaded}");
    println!(
        "# Regenerate with: cargo run -p crabwurcs-pdb --example regenerate_ccd_stereo_references -- DIRECTORY '{downloaded}' > crabwurcs-pdb/tests/data/ccd_stereo_references.tsv"
    );
    println!(
        "# Heavy atoms of the ideal (pdbx_model_Cartn_*_ideal) and model (model_Cartn_*) coordinates."
    );
    println!("# Component lines: ## ID<TAB>SHA-256 of ID.cif<TAB>_chem_comp.name");
    println!("# Columns: component, coordinates (ideal|model), atom, element, x, y, z");
    for path in paths {
        let component = read_component(&path)?;
        println!(
            "## {}\t{}\t{}",
            component.id,
            sha256(&path)?,
            component.name
        );
        for (set, prefix, suffix) in [
            ("ideal", "pdbx_model_Cartn_", "_ideal"),
            ("model", "model_Cartn_", ""),
        ] {
            for atom in &component.atoms {
                let field = |key: &str| {
                    atom.iter()
                        .find(|(name, _)| name == key)
                        .map(|(_, value)| value.as_str())
                        .unwrap_or("?")
                };
                let element = field("type_symbol");
                if matches!(element, "H" | "D") {
                    continue;
                }
                let coordinates = ["x", "y", "z"]
                    .iter()
                    .map(|axis| field(&format!("{prefix}{axis}{suffix}")))
                    .collect::<Vec<_>>();
                if coordinates.contains(&"?") {
                    continue;
                }
                println!(
                    "{}\t{set}\t{}\t{element}\t{}\t{}\t{}",
                    component.id,
                    field("atom_id"),
                    coordinates[0],
                    coordinates[1],
                    coordinates[2]
                );
            }
        }
    }
    Ok(())
}
