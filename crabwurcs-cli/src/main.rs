//! `crabwurcs` command-line interface.
//!
//! Covers the use cases of GlycanFormatConverter-cli (notation
//! conversion), plus CLI entry points for MolWURCS-equivalent chemistry
//! conversion, PDB glycan extraction, and SNFG rendering — one binary
//! instead of four separate tools.
//!
//! Parsing, notation conversion, corpus-backed SMILES conversion, and SNFG
//! rendering are implemented in the workspace crates. The bundled molecular
//! corpus uses a pure-Rust backend; general glycan extraction from previously
//! unseen molecules remains the MolWURCS-specific work in progress.

use clap::{Parser, Subcommand, ValueEnum};
use std::io::Read;
use std::path::PathBuf;
use std::process::ExitCode;

#[derive(Parser)]
#[command(
    name = "crabwurcs",
    version,
    about = "Convert, inspect, and render glycan structures",
    long_about = "crabwurcs is a pure-Rust toolkit that replaces \
                  GlycanFormatConverter(-cli), MolWURCS, and seq2snfg \
                  with a single binary."
)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Clone, Copy, Debug, ValueEnum)]
enum TextFormat {
    Auto,
    Wurcs,
    IupacCondensed,
    IupacExtended,
    Glycam,
    Smiles,
}

impl From<TextFormat> for crabwurcs::Format {
    fn from(value: TextFormat) -> Self {
        match value {
            TextFormat::Auto => Self::Auto,
            TextFormat::Wurcs => Self::Wurcs,
            TextFormat::IupacCondensed => Self::IupacCondensed,
            TextFormat::IupacExtended => Self::IupacExtended,
            TextFormat::Glycam => Self::Glycam,
            TextFormat::Smiles => Self::Smiles,
        }
    }
}

#[derive(Clone, Copy, Debug, ValueEnum)]
enum ChemFormatArg {
    Mol,
    Sdf,
    Smiles,
}

impl From<ChemFormatArg> for crabwurcs::mol::ChemFormat {
    fn from(value: ChemFormatArg) -> Self {
        match value {
            ChemFormatArg::Mol => crabwurcs::mol::ChemFormat::Mol,
            ChemFormatArg::Sdf => crabwurcs::mol::ChemFormat::Sdf,
            ChemFormatArg::Smiles => crabwurcs::mol::ChemFormat::Smiles,
        }
    }
}

#[derive(Subcommand)]
enum Command {
    /// Convert between WURCS, IUPAC, GLYCAM, and SMILES.
    ///
    /// Equivalent to GlycanFormatConverter-cli. GlycoCT is intentionally
    /// not supported — see the workspace README for why.
    Convert {
        #[arg(long, value_enum, default_value = "auto")]
        from: TextFormat,
        #[arg(long, value_enum)]
        to: TextFormat,
        /// Input string, or a file path if `--input-file` is set. Reads
        /// stdin if omitted.
        input: Option<String>,
        #[arg(long)]
        input_file: bool,
    },

    /// Extract WURCS from a chemical structure (MOL/SDF/SMILES).
    /// Equivalent to `MolWURCS --wurcs-from-molecules`.
    MolToWurcs {
        #[arg(long, value_enum)]
        format: ChemFormatArg,
        input: Option<PathBuf>,
    },

    /// Render a chemical structure (MOL/SDF/SMILES) from WURCS.
    /// Equivalent to `MolWURCS --molecules-from-wurcs`.
    WurcsToMol {
        #[arg(long, value_enum)]
        format: ChemFormatArg,
        input: Option<PathBuf>,
    },

    /// Extract glycan structure(s) from a PDB or mmCIF file (one line per
    /// glycan found). WURCS is the default; IUPAC and GLYCAM are available
    /// through `--to`.
    PdbToWurcs {
        input: PathBuf,
        #[arg(long, value_enum, default_value = "wurcs")]
        to: TextFormat,
    },

    /// Render a glycan as an SNFG SVG figure. Equivalent to seq2snfg.
    Render {
        /// WURCS input string, or a file path if `--input-file` is set.
        /// Reads stdin if omitted.
        input: Option<String>,
        #[arg(long)]
        input_file: bool,
        #[arg(short, long)]
        output: Option<PathBuf>,
        #[arg(long, value_enum, default_value = "auto")]
        from: TextFormat,
    },
}

fn main() -> ExitCode {
    let cli = Cli::parse();

    let result = match cli.command {
        Command::Convert {
            from,
            to,
            input,
            input_file,
        } => run_convert(from, to, input, input_file),
        Command::MolToWurcs { format, input } => run_mol_to_wurcs(format, input),
        Command::WurcsToMol { format, input } => run_wurcs_to_mol(format, input),
        Command::PdbToWurcs { input, to } => run_pdb_convert(input, to),
        Command::Render {
            input,
            input_file,
            output,
            from,
        } => run_render(input, input_file, output, from),
    };

    match result {
        Ok(()) => ExitCode::SUCCESS,
        Err(err) => {
            eprintln!("error: {err}");
            ExitCode::FAILURE
        }
    }
}

/// Reads a positional-arg-or-stdin string input, matching the convention
/// used by `convert` and `render`: if `input` is `Some` and `input_file`
/// is set, treat it as a file path; if `input` is `Some` and not a file,
/// treat it as the literal string; if `input` is `None`, read stdin.
fn read_text_input(input: Option<String>, input_file: bool) -> anyhow::Result<String> {
    match input {
        Some(value) if input_file => Ok(std::fs::read_to_string(value)?),
        Some(value) => Ok(value),
        None => {
            let mut buf = String::new();
            std::io::stdin().read_to_string(&mut buf)?;
            Ok(buf)
        }
    }
}

fn run_convert(
    from: TextFormat,
    to: TextFormat,
    input: Option<String>,
    input_file: bool,
) -> anyhow::Result<()> {
    let text = read_text_input(input, input_file)?;

    let output = crabwurcs::convert(&text, from.into(), to.into())?;

    println!("{output}");
    Ok(())
}

fn run_mol_to_wurcs(format: ChemFormatArg, input: Option<PathBuf>) -> anyhow::Result<()> {
    let text = match input {
        Some(path) => std::fs::read_to_string(path)?,
        None => {
            let mut buf = String::new();
            std::io::stdin().read_to_string(&mut buf)?;
            buf
        }
    };

    let graphs = match format {
        ChemFormatArg::Smiles => {
            vec![crabwurcs::parse_notation(&text, crabwurcs::Format::Smiles)?]
        }
        _ => crabwurcs::mol::wurcs_from_molecules(&text, format.into())?,
    };
    for graph in graphs {
        println!("{}", crabwurcs::core::write_wurcs(&graph)?);
    }
    Ok(())
}

fn run_wurcs_to_mol(format: ChemFormatArg, input: Option<PathBuf>) -> anyhow::Result<()> {
    let text = match input {
        Some(path) => std::fs::read_to_string(path)?,
        None => {
            let mut buf = String::new();
            std::io::stdin().read_to_string(&mut buf)?;
            buf
        }
    };

    let graph = crabwurcs::core::parse_wurcs(&text)?;
    let output = match format {
        ChemFormatArg::Smiles => crabwurcs::write_notation(&graph, crabwurcs::Format::Smiles)?,
        _ => crabwurcs::mol::molecule_from_wurcs(&graph, format.into())?,
    };
    println!("{output}");
    Ok(())
}

fn run_pdb_convert(input: PathBuf, to: TextFormat) -> anyhow::Result<()> {
    if matches!(to, TextFormat::Auto) {
        anyhow::bail!("--to auto is not valid for PDB output");
    }
    let glycans = crabwurcs::pdb::extract_glycans_from_file(&input)?;
    for glycan in glycans {
        let notation = crabwurcs::write_notation(&glycan.graph, to.into())?;
        match glycan.attachment_site {
            Some(site) => println!("{site}\t{notation}"),
            None => println!("{notation}"),
        }
    }
    Ok(())
}

fn run_render(
    input: Option<String>,
    input_file: bool,
    output: Option<PathBuf>,
    from_format: TextFormat,
) -> anyhow::Result<()> {
    let text = read_text_input(input, input_file)?;
    let graph = crabwurcs::parse_notation(&text, from_format.into())?;
    let svg = crabwurcs::snfg::render_svg(&graph)?;

    match output {
        Some(path) => std::fs::write(path, svg)?,
        None => println!("{svg}"),
    }
    Ok(())
}
