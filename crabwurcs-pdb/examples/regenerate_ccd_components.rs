//! Deterministically regenerate the bundled CCD carbohydrate lookup table.
//!
//! Usage:
//! `cargo run -p crabwurcs-pdb --example regenerate_ccd_components -- components.cif.gz [output.tsv]`
//!
//! The command reads a local, pinned wwPDB CCD snapshot. It does not download
//! data. Gzip input is decoded with the platform `gzip` command; uncompressed
//! CIF is read directly. Only components declared as saccharides and accepted
//! by the molecular recognizer are emitted.

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use std::process::Command;

use crabwurcs_core::{AnomericSymbol, RingClosure, classify_residue};
use crabwurcs_mol::{ChemFormat, wurcs_from_molecule};

fn read_snapshot(path: &Path) -> Result<String, String> {
    if path.extension().is_some_and(|extension| extension == "gz") {
        let output = Command::new("gzip")
            .args(["-dc"])
            .arg(path)
            .output()
            .map_err(|error| format!("cannot run gzip: {error}"))?;
        if !output.status.success() {
            return Err(format!("gzip failed with {}", output.status));
        }
        String::from_utf8(output.stdout).map_err(|error| format!("CCD is not UTF-8: {error}"))
    } else {
        std::fs::read_to_string(path).map_err(|error| format!("cannot read CCD: {error}"))
    }
}

fn snapshot_sha256(path: &Path) -> Result<String, String> {
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

fn cif_tokens(input: &str) -> Vec<String> {
    let bytes = input.as_bytes();
    let mut tokens = Vec::new();
    let mut cursor = 0;
    let mut line_start = true;
    while cursor < bytes.len() {
        if bytes[cursor].is_ascii_whitespace() {
            line_start = bytes[cursor] == b'\n';
            cursor += 1;
            continue;
        }
        if bytes[cursor] == b'#' {
            cursor += 1;
            while cursor < bytes.len() && bytes[cursor] != b'\n' {
                cursor += 1;
            }
            continue;
        }
        if line_start && bytes[cursor] == b';' {
            cursor += 1;
            let start = cursor;
            while cursor + 1 < bytes.len() && !(bytes[cursor] == b'\n' && bytes[cursor + 1] == b';')
            {
                cursor += 1;
            }
            tokens.push(String::from_utf8_lossy(&bytes[start..cursor]).into_owned());
            cursor = (cursor + 2).min(bytes.len());
            line_start = false;
            continue;
        }
        line_start = false;
        if matches!(bytes[cursor], b'\'' | b'"') {
            let quote = bytes[cursor];
            cursor += 1;
            let start = cursor;
            while cursor < bytes.len() && bytes[cursor] != quote {
                cursor += 1;
            }
            tokens.push(String::from_utf8_lossy(&bytes[start..cursor]).into_owned());
            cursor = (cursor + 1).min(bytes.len());
            continue;
        }
        let start = cursor;
        while cursor < bytes.len() && !bytes[cursor].is_ascii_whitespace() {
            cursor += 1;
        }
        tokens.push(String::from_utf8_lossy(&bytes[start..cursor]).into_owned());
    }
    tokens
}

#[derive(Default)]
struct Component {
    kind: String,
    descriptors: Vec<(String, String, String)>,
}

fn parse_components(input: &str) -> BTreeMap<String, Component> {
    let tokens = cif_tokens(input);
    let mut components = BTreeMap::<String, Component>::new();
    let mut component_id = String::new();
    let mut cursor = 0;
    while cursor < tokens.len() {
        let token = &tokens[cursor];
        if let Some(id) = token.strip_prefix("data_") {
            component_id = id.to_ascii_uppercase();
            components.entry(component_id.clone()).or_default();
            cursor += 1;
            continue;
        }
        if token == "_chem_comp.type" && cursor + 1 < tokens.len() {
            components.entry(component_id.clone()).or_default().kind = tokens[cursor + 1].clone();
            cursor += 2;
            continue;
        }
        if token != "loop_" {
            cursor += 1;
            continue;
        }
        cursor += 1;
        let header_start = cursor;
        while cursor < tokens.len() && tokens[cursor].starts_with('_') {
            cursor += 1;
        }
        let headers = &tokens[header_start..cursor];
        if headers.is_empty() {
            continue;
        }
        let descriptor = |name: &str| headers.iter().position(|header| header == name);
        let type_column = descriptor("_pdbx_chem_comp_descriptor.type");
        let program_column = descriptor("_pdbx_chem_comp_descriptor.program");
        let value_column = descriptor("_pdbx_chem_comp_descriptor.descriptor");
        while cursor + headers.len() <= tokens.len()
            && !tokens[cursor].starts_with('_')
            && tokens[cursor] != "loop_"
            && !tokens[cursor].starts_with("data_")
            && tokens[cursor] != "stop_"
        {
            let row = &tokens[cursor..cursor + headers.len()];
            if let (Some(kind), Some(program), Some(value)) =
                (type_column, program_column, value_column)
            {
                components
                    .entry(component_id.clone())
                    .or_default()
                    .descriptors
                    .push((row[kind].clone(), row[program].clone(), row[value].clone()));
            }
            cursor += headers.len();
        }
        if cursor < tokens.len() && tokens[cursor] == "stop_" {
            cursor += 1;
        }
    }
    components
}

fn run() -> Result<(), String> {
    let mut arguments = std::env::args_os().skip(1);
    let snapshot = PathBuf::from(
        arguments
            .next()
            .ok_or("usage: regenerate_ccd_components <components.cif[.gz]> [output.tsv]")?,
    );
    let output = arguments
        .next()
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from("crabwurcs-pdb/data/pdb_carbohydrate_components.tsv"));
    let contents = read_snapshot(&snapshot)?;
    let checksum = snapshot_sha256(&snapshot)?;
    let mut rows = BTreeMap::new();
    for (id, component) in parse_components(&contents) {
        if !component.kind.to_ascii_lowercase().contains("saccharide") {
            continue;
        }
        let smiles = component
            .descriptors
            .iter()
            .filter(|(kind, _, _)| kind == "SMILES_CANONICAL")
            .min_by_key(|(_, program, _)| (!program.eq_ignore_ascii_case("OpenEye"), program))
            .map(|(_, _, value)| value);
        let Some(smiles) = smiles else {
            continue;
        };
        let Ok(graph) = wurcs_from_molecule(smiles, ChemFormat::Smiles) else {
            continue;
        };
        if graph.node_count() != 1 {
            continue;
        }
        let residue = graph
            .root()
            .and_then(|root| graph.residue(root))
            .ok_or_else(|| format!("{id}: recognizer returned a rootless graph"))?;
        let Some(kind) = classify_residue(residue) else {
            continue;
        };
        let anomer = match residue.anomeric_symbol {
            AnomericSymbol::Alpha => "a",
            AnomericSymbol::Beta => "b",
            _ => "x",
        };
        let ring = match residue.ring {
            RingClosure::Pyranose => "p",
            RingClosure::Furanose => "f",
            _ => "x",
        };
        let row = format!("{}\t{anomer}\t{ring}\tn", kind.canonical_name());
        if rows.insert(id.clone(), row).is_some() {
            return Err(format!("conflicting component mapping for {id}"));
        }
    }
    let mut generated = format!(
        "# crabWURCS bundled wwPDB carbohydrate component map\n\
         # Source file: {}\n\
         # Source SHA-256: {checksum}\n\
         # Generated by crabwurcs-pdb/examples/regenerate_ccd_components.rs\n\
         # Columns: component_id, crabwurcs_residue_kind, anomer(a|b|x), ring(p|f|x), mirror(y|n)\n",
        snapshot.display()
    );
    for (id, row) in rows {
        generated.push_str(&format!("{id}\t{row}\n"));
    }
    std::fs::write(&output, generated)
        .map_err(|error| format!("cannot write {}: {error}", output.display()))?;
    eprintln!("wrote {}", output.display());
    Ok(())
}

fn main() {
    if let Err(error) = run() {
        eprintln!("error: {error}");
        std::process::exit(1);
    }
}
