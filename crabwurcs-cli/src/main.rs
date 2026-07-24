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
enum MotifTextFormat {
    Auto,
    Wurcs,
    IupacCondensed,
    IupacExtended,
    Glycam,
}

impl From<MotifTextFormat> for crabwurcs::Format {
    fn from(value: MotifTextFormat) -> Self {
        match value {
            MotifTextFormat::Auto => Self::Auto,
            MotifTextFormat::Wurcs => Self::Wurcs,
            MotifTextFormat::IupacCondensed => Self::IupacCondensed,
            MotifTextFormat::IupacExtended => Self::IupacExtended,
            MotifTextFormat::Glycam => Self::Glycam,
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

    /// Render a glycan as an SNFG SVG or PNG figure. Equivalent to seq2snfg.
    Render {
        /// Glycan notation string, or a file path if `--input-file` is set.
        /// Reads stdin if omitted.
        input: Option<String>,
        #[arg(long)]
        input_file: bool,
        #[arg(short, long)]
        output: Option<PathBuf>,
        #[arg(long, value_enum, default_value = "auto")]
        from: TextFormat,
        /// Motif pattern to highlight. May be repeated; every occurrence of
        /// every motif is highlighted.
        #[arg(long = "highlight-motif")]
        highlight_motifs: Vec<String>,
        /// Notation used by every --highlight-motif value.
        #[arg(long, value_enum, default_value = "auto")]
        motif_from: MotifTextFormat,
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
            highlight_motifs,
            motif_from,
        } => run_render(
            input,
            input_file,
            output,
            from,
            highlight_motifs,
            motif_from,
        ),
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
    highlight_motifs: Vec<String>,
    motif_from: MotifTextFormat,
) -> anyhow::Result<()> {
    let text = read_text_input(input, input_file)?;
    let requested_format: crabwurcs::Format = from_format.into();
    let actual_format = if requested_format == crabwurcs::Format::Auto {
        crabwurcs::detect_format(&text)
    } else {
        requested_format
    };
    let graph = crabwurcs::parse_notation(&text, actual_format)?;
    let render_options = crabwurcs::snfg::RenderOptions {
        source_notation: Some(crabwurcs::snfg::SourceNotation::new(
            format_name(actual_format),
            text.trim(),
        )),
        ..crabwurcs::snfg::RenderOptions::default()
    };
    let output_kind = output
        .as_ref()
        .map(|path| render_output_kind(path))
        .transpose()?;
    let motifs = highlight_motifs
        .iter()
        .enumerate()
        .map(|(index, motif)| {
            crabwurcs::parse_notation(motif, motif_from.into()).map_err(|error| {
                anyhow::anyhow!(
                    "failed to parse highlight motif {} ({motif:?}): {error}",
                    index + 1
                )
            })
        })
        .collect::<anyhow::Result<Vec<_>>>()?;
    for (index, motif) in motifs.iter().enumerate() {
        crabwurcs::find_motif_matches(&graph, motif).map_err(|error| {
            anyhow::anyhow!(
                "unsupported highlight motif {} ({:?}): {error}",
                index + 1,
                highlight_motifs[index]
            )
        })?;
    }
    match (output, output_kind) {
        (None, None) => {
            let svg = if motifs.is_empty() {
                crabwurcs::snfg::render_svg_with_options(&graph, &render_options)?
            } else {
                crabwurcs::snfg::render_svg_with_motifs(&graph, &motifs, &render_options)?
            };
            println!("{svg}");
        }
        (Some(path), Some(RenderOutputKind::Svg)) => {
            let svg = if motifs.is_empty() {
                crabwurcs::snfg::render_svg_with_options(&graph, &render_options)?
            } else {
                crabwurcs::snfg::render_svg_with_motifs(&graph, &motifs, &render_options)?
            };
            std::fs::write(path, svg)?;
        }
        (Some(path), Some(RenderOutputKind::Png)) => {
            let png = if motifs.is_empty() {
                crabwurcs::snfg::render_png_with_options(&graph, &render_options)?
            } else {
                crabwurcs::snfg::render_png_with_motifs(&graph, &motifs, &render_options)?
            };
            std::fs::write(path, png)?;
        }
        _ => unreachable!("output path and inferred kind must be present together"),
    }
    Ok(())
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum RenderOutputKind {
    Svg,
    Png,
}

fn render_output_kind(path: &std::path::Path) -> anyhow::Result<RenderOutputKind> {
    let extension = path
        .extension()
        .and_then(|extension| extension.to_str())
        .ok_or_else(|| {
            anyhow::anyhow!(
                "render output path must end in .svg or .png: {}",
                path.display()
            )
        })?;
    if extension.eq_ignore_ascii_case("svg") {
        Ok(RenderOutputKind::Svg)
    } else if extension.eq_ignore_ascii_case("png") {
        Ok(RenderOutputKind::Png)
    } else {
        anyhow::bail!("unsupported render output extension .{extension}; expected .svg or .png")
    }
}

fn format_name(format: crabwurcs::Format) -> &'static str {
    match format {
        crabwurcs::Format::Auto => "auto",
        crabwurcs::Format::Wurcs => "wurcs",
        crabwurcs::Format::IupacCondensed => "iupac-condensed",
        crabwurcs::Format::IupacExtended => "iupac-extended",
        crabwurcs::Format::Glycam => "glycam",
        crabwurcs::Format::Smiles => "smiles",
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn render_accepts_repeatable_motifs_and_an_explicit_motif_format() {
        let cli = Cli::try_parse_from([
            "crabwurcs",
            "render",
            "--from",
            "iupac-condensed",
            "--highlight-motif",
            "Fuc(a1-?)GlcNAc",
            "--highlight-motif",
            "Gal(b1-?)GlcNAc",
            "--motif-from",
            "iupac-condensed",
            "Gal(b1-4)GlcNAc",
        ])
        .unwrap();
        let Command::Render {
            highlight_motifs,
            motif_from,
            ..
        } = cli.command
        else {
            panic!("expected render command");
        };
        assert_eq!(highlight_motifs.len(), 2);
        assert!(matches!(motif_from, MotifTextFormat::IupacCondensed));
    }

    #[test]
    fn smiles_is_not_an_accepted_motif_format() {
        assert!(
            Cli::try_parse_from([
                "crabwurcs",
                "render",
                "--highlight-motif",
                "C1CC1",
                "--motif-from",
                "smiles",
                "Glc",
            ])
            .is_err()
        );
    }

    #[test]
    fn render_output_format_is_inferred_case_insensitively() {
        assert_eq!(
            render_output_kind(std::path::Path::new("figure.svg")).unwrap(),
            RenderOutputKind::Svg
        );
        assert_eq!(
            render_output_kind(std::path::Path::new("figure.SVG")).unwrap(),
            RenderOutputKind::Svg
        );
        assert_eq!(
            render_output_kind(std::path::Path::new("figure.png")).unwrap(),
            RenderOutputKind::Png
        );
        assert_eq!(
            render_output_kind(std::path::Path::new("figure.PNG")).unwrap(),
            RenderOutputKind::Png
        );
        assert!(render_output_kind(std::path::Path::new("figure")).is_err());
        assert!(render_output_kind(std::path::Path::new("figure.jpg")).is_err());
    }

    #[test]
    fn render_writes_and_overwrites_svg_and_png_files() {
        let directory =
            std::env::temp_dir().join(format!("crabwurcs-render-test-{}", std::process::id()));
        std::fs::create_dir_all(&directory).unwrap();
        let svg_path = directory.join("glycan.SVG");
        let png_path = directory.join("glycan.PNG");
        std::fs::write(&svg_path, "old").unwrap();
        std::fs::write(&png_path, "old").unwrap();

        run_render(
            Some("Glc".into()),
            false,
            Some(svg_path.clone()),
            TextFormat::IupacCondensed,
            Vec::new(),
            MotifTextFormat::Auto,
        )
        .unwrap();
        run_render(
            Some("Glc".into()),
            false,
            Some(png_path.clone()),
            TextFormat::IupacCondensed,
            Vec::new(),
            MotifTextFormat::Auto,
        )
        .unwrap();

        let svg = std::fs::read_to_string(&svg_path).unwrap();
        assert!(svg.starts_with("<svg"));
        assert!(svg.contains("<metadata id=\"crabwurcs-notations\">"));
        let png = std::fs::read(&png_path).unwrap();
        assert_eq!(&png[..8], b"\x89PNG\r\n\x1a\n");

        std::fs::remove_file(svg_path).unwrap();
        std::fs::remove_file(png_path).unwrap();
        std::fs::remove_dir(directory).unwrap();
    }
}
