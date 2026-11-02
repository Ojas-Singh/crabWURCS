//! `crabwurcs` command-line interface.
//!
//! Covers the use cases of GlycanFormatConverter-cli (notation
//! conversion), plus CLI entry points for MolWURCS-equivalent chemistry
//! conversion, PDB glycan extraction, and SNFG rendering — one binary
//! instead of four separate tools.
//!
//! Parsing, notation conversion, de-novo molecular interop, PDB/mmCIF
//! extraction, and SNFG rendering are implemented in the workspace crates.

use clap::{Args, Parser, Subcommand, ValueEnum};
use crabwurcs::{
    Adduct, GridOptions, HighlightSelection, MassOptions, MotifPosition, MotifQuery, Orientation,
    ReducingEnd, RenderOptions,
};
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

/// Use of the bundled table of GlycoShape notations by `convert`.
#[derive(Clone, Copy, Debug, Default, ValueEnum)]
enum TableArg {
    /// Compute every conversion; never consult the table.
    #[default]
    Off,
    /// Answer from the table only when the input cannot be converted.
    Fallback,
    /// Compute, and fail when the input's table record holds another structure.
    Check,
}

impl From<TableArg> for crabwurcs::TableMode {
    fn from(value: TableArg) -> Self {
        match value {
            TableArg::Off => Self::Off,
            TableArg::Fallback => Self::Fallback,
            TableArg::Check => Self::Check,
        }
    }
}

/// What decides when a residue name and the coordinates disagree.
#[derive(Clone, Copy, Debug, Default, ValueEnum)]
enum PolicyArg {
    /// The coordinates decide; names fill in only what the coordinates
    /// cannot show (missing atoms, distorted centres).
    #[default]
    Coordinates,
    /// Names decide wherever they say something.
    Names,
    /// The coordinates decide, and a disagreement is an error.
    Strict,
}

impl From<PolicyArg> for crabwurcs::pdb::NamePolicy {
    fn from(value: PolicyArg) -> Self {
        match value {
            PolicyArg::Coordinates => Self::Coordinates,
            PolicyArg::Names => Self::Names,
            PolicyArg::Strict => Self::Strict,
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

/// How a glycan is drawn, shared by `render` and `grid`.
#[derive(Args, Clone, Debug, Default)]
struct DrawingArgs {
    /// Where the reducing end goes: root-right (the SNFG convention),
    /// root-left, root-bottom or root-top.
    #[arg(long)]
    orientation: Option<Orientation>,
    /// Tighter spacing and no linkage text, for heatmap axes and tables.
    #[arg(long)]
    compact: bool,
    /// Text drawn past the reducing end, such as Asn, Ser/Thr or 2-AB.
    #[arg(long, value_name = "TEXT")]
    reducing_end: Option<String>,
    /// Draw repeat units in square brackets with their repeat count.
    #[arg(long)]
    repeat_brackets: bool,
    /// Draw undefined antennae and substituents above the structure.
    #[arg(long)]
    floating_fragments: bool,
}

impl DrawingArgs {
    fn options(&self) -> RenderOptions {
        RenderOptions {
            orientation: self.orientation.unwrap_or_default(),
            compact: self.compact,
            reducing_end: self.reducing_end.clone(),
            repeat_brackets: self.repeat_brackets,
            floating_fragments: self.floating_fragments,
            ..RenderOptions::default()
        }
    }
}

/// What `render` highlights besides `--highlight-motif`.
#[derive(Args, Clone, Debug, Default)]
struct HighlightArgs {
    /// Comma-separated positions for one --highlight-motif, in the order its
    /// residues are read (terminal, internal, reducing, flexible). Given
    /// once per motif, in the same order; motifs after the last one keep the
    /// positions of their name, or none.
    #[arg(long = "motif-positions", value_name = "POSITIONS")]
    motif_positions: Vec<String>,
    /// Residue to highlight, by its node index. May be repeated.
    #[arg(long = "highlight-residue", value_name = "INDEX")]
    highlight_residues: Vec<usize>,
    /// Fragment ion to mark on its bond, such as B2 or Y3α. May be repeated.
    #[arg(long = "cleavage", value_name = "ION")]
    cleavages: Vec<String>,
}

/// What a glycan is derivatized with, for `mass` and `fragments`.
#[derive(Args, Clone, Debug, Default)]
struct MassArgs {
    /// free, reduced, 2-ab, 2-aa or procainamide.
    #[arg(long, default_value = "free")]
    reducing_end: ReducingEnd,
    /// Methylate every free OH, NH and carboxyl group.
    #[arg(long)]
    permethylated: bool,
}

impl MassArgs {
    fn options(&self) -> MassOptions {
        MassOptions {
            reducing_end: self.reducing_end,
            permethylated: self.permethylated,
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
        /// Use of the bundled table of GlycoShape notations. Conversions
        /// are computed; `fallback` answers inputs that cannot be converted
        /// from the table, `check` fails when the input's table record holds
        /// another structure.
        #[arg(long, value_enum, default_value = "off")]
        table: TableArg,
    },

    /// Normalize a WURCS 2.0 string using the local canonicalizer.
    NormalizeWurcs {
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
    /// through `--to`. Every monosaccharide is read from its coordinates;
    /// where a residue name disagrees, a warning goes to standard error.
    PdbToWurcs {
        input: PathBuf,
        #[arg(long, value_enum, default_value = "wurcs")]
        to: TextFormat,
        /// What decides when a residue name and the coordinates disagree.
        #[arg(long, value_enum, default_value = "coordinates")]
        policy: PolicyArg,
        /// Do not print warnings.
        #[arg(long)]
        quiet: bool,
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
        #[command(flatten)]
        highlights: HighlightArgs,
        #[command(flatten)]
        drawing: DrawingArgs,
    },

    /// Draw many glycans on one SVG or PNG canvas. The input file holds one
    /// notation per line, with an optional tab and caption.
    Grid {
        input: PathBuf,
        #[arg(short, long)]
        output: PathBuf,
        #[arg(long, value_enum, default_value = "auto")]
        from: TextFormat,
        #[arg(long, default_value_t = 4)]
        columns: usize,
        #[command(flatten)]
        drawing: DrawingArgs,
    },

    /// Print the composition, formula and masses of a glycan.
    Mass {
        /// Glycan notation string, or a file path if `--input-file` is set.
        /// Reads stdin if omitted.
        input: Option<String>,
        #[arg(long)]
        input_file: bool,
        #[arg(long, value_enum, default_value = "auto")]
        from: TextFormat,
        #[command(flatten)]
        mass: MassArgs,
        /// Also print the m/z with this adduct: H, Na, K, NH4, Li, -H, Cl,
        /// HCOO or CH3COO. May be repeated.
        #[arg(long = "adduct", allow_hyphen_values = true)]
        adducts: Vec<Adduct>,
        /// Adducts per ion, which is also the charge.
        #[arg(long, default_value_t = 1)]
        charge: u32,
    },

    /// List the B, C, Y and Z ions of every glycosidic cleavage.
    Fragments {
        input: Option<String>,
        #[arg(long)]
        input_file: bool,
        #[arg(long, value_enum, default_value = "auto")]
        from: TextFormat,
        #[command(flatten)]
        mass: MassArgs,
        #[arg(long, default_value = "H", allow_hyphen_values = true)]
        adduct: Adduct,
    },

    /// Count the occurrences of a motif, given as a notation or a name from
    /// the library (see `motifs`). Fails if it does not occur.
    Motif {
        motif: String,
        input: Option<String>,
        #[arg(long)]
        input_file: bool,
        #[arg(long, value_enum, default_value = "auto")]
        from: TextFormat,
        #[arg(long, value_enum, default_value = "auto")]
        motif_from: MotifTextFormat,
        /// Comma-separated positions of the motif's residues.
        #[arg(long)]
        positions: Option<String>,
    },

    /// List the motif library, or the library motifs a glycan holds.
    Motifs {
        input: Option<String>,
        #[arg(long)]
        input_file: bool,
        #[arg(long, value_enum, default_value = "auto")]
        from: TextFormat,
    },

    /// Compare two glycans, in any notations. Fails if they differ.
    Same {
        left: String,
        right: String,
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
            table,
        } => run_convert(from, to, input, input_file, table),
        Command::NormalizeWurcs { input, input_file } => run_normalize_wurcs(input, input_file),
        Command::MolToWurcs { format, input } => run_mol_to_wurcs(format, input),
        Command::WurcsToMol { format, input } => run_wurcs_to_mol(format, input),
        Command::PdbToWurcs {
            input,
            to,
            policy,
            quiet,
        } => run_pdb_convert(input, to, policy, quiet),
        Command::Render {
            input,
            input_file,
            output,
            from,
            highlight_motifs,
            motif_from,
            highlights,
            drawing,
        } => run_render(
            input,
            input_file,
            output,
            from,
            highlight_motifs,
            motif_from,
            &highlights,
            &drawing,
        ),
        Command::Grid {
            input,
            output,
            from,
            columns,
            drawing,
        } => run_grid(&input, &output, from, columns, &drawing),
        Command::Mass {
            input,
            input_file,
            from,
            mass,
            adducts,
            charge,
        } => run_mass(input, input_file, from, &mass, &adducts, charge),
        Command::Fragments {
            input,
            input_file,
            from,
            mass,
            adduct,
        } => run_fragments(input, input_file, from, &mass, adduct),
        Command::Motif {
            motif,
            input,
            input_file,
            from,
            motif_from,
            positions,
        } => run_motif(&motif, input, input_file, from, motif_from, positions),
        Command::Motifs {
            input,
            input_file,
            from,
        } => run_motifs(input, input_file, from),
        Command::Same { left, right, from } => run_same(&left, &right, from),
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
    table: TableArg,
) -> anyhow::Result<()> {
    let text = read_text_input(input, input_file)?;

    let output = crabwurcs::convert_with_table(&text, from.into(), to.into(), table.into())?;

    println!("{output}");
    Ok(())
}

fn run_normalize_wurcs(input: Option<String>, input_file: bool) -> anyhow::Result<()> {
    let text = read_text_input(input, input_file)?;
    println!("{}", crabwurcs::normalize_wurcs(&text)?);
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

fn run_pdb_convert(
    input: PathBuf,
    to: TextFormat,
    policy: PolicyArg,
    quiet: bool,
) -> anyhow::Result<()> {
    if matches!(to, TextFormat::Auto) {
        anyhow::bail!("--to auto is not valid for PDB output");
    }
    let reading = crabwurcs::pdb::read_glycans_from_file(
        &input,
        &crabwurcs::pdb::PdbReadOptions::with_policy(policy.into()),
    )?;
    if !quiet {
        for warning in &reading.warnings {
            eprintln!("warning: {warning}");
        }
    }
    for glycan in reading.glycans {
        let notation = crabwurcs::write_notation(&glycan.graph, to.into())?;
        match glycan.attachment_site {
            Some(site) => println!("{site}\t{notation}"),
            None => println!("{notation}"),
        }
    }
    Ok(())
}

fn parse_positions(text: &str) -> anyhow::Result<Vec<MotifPosition>> {
    text.split(',')
        .map(|position| position.parse().map_err(anyhow::Error::msg))
        .collect()
}

#[allow(clippy::too_many_arguments)]
fn run_render(
    input: Option<String>,
    input_file: bool,
    output: Option<PathBuf>,
    from_format: TextFormat,
    highlight_motifs: Vec<String>,
    motif_from: MotifTextFormat,
    highlights: &HighlightArgs,
    drawing: &DrawingArgs,
) -> anyhow::Result<()> {
    let text = read_text_input(input, input_file)?;
    let requested_format: crabwurcs::Format = from_format.into();
    let actual_format = if requested_format == crabwurcs::Format::Auto {
        crabwurcs::detect_format(&text)
    } else {
        requested_format
    };
    let graph = crabwurcs::parse_notation(&text, actual_format)?;
    let mut render_options = RenderOptions {
        source_notation: Some(crabwurcs::snfg::SourceNotation::new(
            format_name(actual_format),
            text.trim(),
        )),
        ..drawing.options()
    };
    if !highlights.cleavages.is_empty() {
        let fragments =
            crabwurcs::glycan_fragments(&graph, &MassOptions::default(), Adduct::Proton)?;
        for ion in &highlights.cleavages {
            let fragment = fragments
                .iter()
                .find(|fragment| fragment.label == *ion)
                .ok_or_else(|| anyhow::anyhow!("this glycan has no fragment ion {ion}"))?;
            render_options
                .cleavages
                .push(crabwurcs::cleavage_mark(fragment));
        }
    }
    let output_kind = output
        .as_ref()
        .map(|path| render_output_kind(path))
        .transpose()?;
    if highlights.motif_positions.len() > highlight_motifs.len() {
        anyhow::bail!("--motif-positions was given more often than --highlight-motif");
    }
    let mut selection = None::<HighlightSelection>;
    for (index, motif) in highlight_motifs.iter().enumerate() {
        let mut query = crabwurcs::resolve_motif(motif, motif_from.into()).map_err(|error| {
            anyhow::anyhow!(
                "failed to parse highlight motif {} ({motif:?}): {error}",
                index + 1
            )
        })?;
        if let Some(positions) = highlights.motif_positions.get(index) {
            query = MotifQuery::with_positions(query.motif().clone(), &parse_positions(positions)?)
                .map_err(|error| {
                    anyhow::anyhow!("highlight motif {} ({motif:?}): {error}", index + 1)
                })?;
        }
        let matches = query.find(&graph).map_err(|error| {
            anyhow::anyhow!(
                "unsupported highlight motif {} ({motif:?}): {error}",
                index + 1
            )
        })?;
        let found = HighlightSelection::from_matches(&matches);
        let union = selection.get_or_insert_with(HighlightSelection::default);
        union.node_indices.extend(found.node_indices);
        union.edge_indices.extend(found.edge_indices);
    }
    if !highlights.highlight_residues.is_empty() {
        let found =
            HighlightSelection::from_residues(&graph, highlights.highlight_residues.clone())?;
        let union = selection.get_or_insert_with(HighlightSelection::default);
        union.node_indices.extend(found.node_indices);
        union.edge_indices.extend(found.edge_indices);
    }
    let render_svg = || match &selection {
        Some(selection) => {
            crabwurcs::snfg::render_svg_with_selection(&graph, selection, &render_options)
        }
        None => crabwurcs::snfg::render_svg_with_options(&graph, &render_options),
    };
    match (output, output_kind) {
        (None, None) => println!("{}", render_svg()?),
        (Some(path), Some(RenderOutputKind::Svg)) => std::fs::write(path, render_svg()?)?,
        (Some(path), Some(RenderOutputKind::Png)) => {
            let png = match &selection {
                Some(selection) => {
                    crabwurcs::snfg::render_png_with_selection(&graph, selection, &render_options)?
                }
                None => crabwurcs::snfg::render_png_with_options(&graph, &render_options)?,
            };
            std::fs::write(path, png)?;
        }
        _ => unreachable!("output path and inferred kind must be present together"),
    }
    Ok(())
}

fn run_grid(
    input: &std::path::Path,
    output: &std::path::Path,
    from: TextFormat,
    columns: usize,
    drawing: &DrawingArgs,
) -> anyhow::Result<()> {
    let kind = render_output_kind(output)?;
    let text = std::fs::read_to_string(input)?;
    let mut graphs = Vec::new();
    let mut labels = Vec::new();
    for line in text.lines().filter(|line| !line.trim().is_empty()) {
        let (notation, label) = line.split_once('\t').unwrap_or((line, ""));
        graphs.push(crabwurcs::parse_notation(notation.trim(), from.into())?);
        labels.push(label.trim().to_owned());
    }
    let options = GridOptions {
        columns,
        labels,
        render: drawing.options(),
    };
    match kind {
        RenderOutputKind::Svg => {
            std::fs::write(output, crabwurcs::render_grid_svg(&graphs, &options)?)?;
        }
        RenderOutputKind::Png => {
            std::fs::write(output, crabwurcs::render_grid_png(&graphs, &options)?)?;
        }
    }
    Ok(())
}

fn read_graph(
    input: Option<String>,
    input_file: bool,
    from: TextFormat,
) -> anyhow::Result<crabwurcs::ResidueGraph> {
    Ok(crabwurcs::parse_notation(
        read_text_input(input, input_file)?.trim(),
        from.into(),
    )?)
}

fn run_mass(
    input: Option<String>,
    input_file: bool,
    from: TextFormat,
    mass: &MassArgs,
    adducts: &[Adduct],
    charge: u32,
) -> anyhow::Result<()> {
    let graph = read_graph(input, input_file, from)?;
    let found = crabwurcs::glycan_mass(&graph, &mass.options())?;
    println!("composition\t{}", crabwurcs::glycan_composition(&graph));
    println!("formula\t{}", found.formula);
    println!("monoisotopic\t{:.4}", found.monoisotopic);
    println!("average\t{:.4}", found.average);
    for adduct in adducts {
        println!(
            "mz {} x{charge}\t{:.4}",
            adduct.as_str(),
            found.mz(*adduct, charge)?
        );
    }
    Ok(())
}

fn run_fragments(
    input: Option<String>,
    input_file: bool,
    from: TextFormat,
    mass: &MassArgs,
    adduct: Adduct,
) -> anyhow::Result<()> {
    let graph = read_graph(input, input_file, from)?;
    for fragment in crabwurcs::glycan_fragments(&graph, &mass.options(), adduct)? {
        println!(
            "{}\t{:.4}\t{}",
            fragment.label, fragment.mz, fragment.formula
        );
    }
    Ok(())
}

fn run_motif(
    motif: &str,
    input: Option<String>,
    input_file: bool,
    from: TextFormat,
    motif_from: MotifTextFormat,
    positions: Option<String>,
) -> anyhow::Result<()> {
    let graph = read_graph(input, input_file, from)?;
    let mut query = crabwurcs::resolve_motif(motif, motif_from.into())?;
    if let Some(positions) = positions {
        query = MotifQuery::with_positions(query.motif().clone(), &parse_positions(&positions)?)?;
    }
    let count = query.count_in(&graph)?;
    println!("{count}");
    if count == 0 {
        anyhow::bail!("the glycan does not contain {motif}");
    }
    Ok(())
}

fn run_motifs(input: Option<String>, input_file: bool, from: TextFormat) -> anyhow::Result<()> {
    if input.is_none() && !input_file {
        for motif in crabwurcs::named_motifs() {
            println!("{}\t{}", motif.name, motif.iupac);
        }
        return Ok(());
    }
    let graph = read_graph(input, input_file, from)?;
    for name in crabwurcs::annotate_motifs(&graph) {
        println!("{name}");
    }
    Ok(())
}

fn run_same(left: &str, right: &str, from: TextFormat) -> anyhow::Result<()> {
    let left = crabwurcs::parse_notation(left, from.into())?;
    let right = crabwurcs::parse_notation(right, from.into())?;
    let same = crabwurcs::same_structure(&left, &right);
    println!("{same}");
    if !same {
        anyhow::bail!("the glycans differ");
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
    fn pdb_to_wurcs_takes_a_name_policy() {
        let default = Cli::try_parse_from(["crabwurcs", "pdb-to-wurcs", "model.pdb"]).unwrap();
        assert!(matches!(
            default.command,
            Command::PdbToWurcs {
                policy: PolicyArg::Coordinates,
                quiet: false,
                ..
            }
        ));
        for (value, expected) in [("names", "Names"), ("strict", "Strict")] {
            let cli =
                Cli::try_parse_from(["crabwurcs", "pdb-to-wurcs", "--policy", value, "model.pdb"])
                    .unwrap();
            let Command::PdbToWurcs { policy, .. } = cli.command else {
                panic!("expected pdb-to-wurcs");
            };
            assert_eq!(format!("{policy:?}"), expected);
        }
        assert!(
            Cli::try_parse_from(["crabwurcs", "pdb-to-wurcs", "--policy", "x", "a.pdb"]).is_err()
        );
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
    fn normalize_wurcs_command_accepts_literal_or_file_input() {
        let literal = Cli::try_parse_from([
            "crabwurcs",
            "normalize-wurcs",
            "WURCS=2.0/1,1,0/[u2122h]/1/",
        ])
        .unwrap();
        assert!(matches!(literal.command, Command::NormalizeWurcs { .. }));

        let file = Cli::try_parse_from([
            "crabwurcs",
            "normalize-wurcs",
            "--input-file",
            "glycan.wurcs",
        ])
        .unwrap();
        assert!(matches!(
            file.command,
            Command::NormalizeWurcs {
                input_file: true,
                ..
            }
        ));
    }

    #[test]
    fn convert_computes_unless_the_table_is_requested() {
        let default =
            Cli::try_parse_from(["crabwurcs", "convert", "--to", "wurcs", "Glc"]).unwrap();
        assert!(matches!(
            default.command,
            Command::Convert {
                table: TableArg::Off,
                ..
            }
        ));
        for (value, expected) in [("fallback", "Fallback"), ("check", "Check")] {
            let cli = Cli::try_parse_from([
                "crabwurcs",
                "convert",
                "--to",
                "wurcs",
                "--table",
                value,
                "Glc",
            ])
            .unwrap();
            let Command::Convert { table, .. } = cli.command else {
                panic!("expected convert command");
            };
            assert_eq!(format!("{table:?}"), expected);
        }
        assert!(
            run_convert(
                TextFormat::Auto,
                TextFormat::Wurcs,
                Some("G60371D-N".into()),
                false,
                TableArg::Off
            )
            .is_err()
        );
        assert!(
            run_convert(
                TextFormat::Auto,
                TextFormat::Wurcs,
                Some("G60371D-N".into()),
                false,
                TableArg::Fallback
            )
            .is_ok()
        );
    }

    fn scratch(name: &str) -> PathBuf {
        let directory =
            std::env::temp_dir().join(format!("crabwurcs-{name}-test-{}", std::process::id()));
        std::fs::create_dir_all(&directory).unwrap();
        directory
    }

    #[test]
    fn render_takes_the_drawing_flags() {
        let cli = Cli::try_parse_from([
            "crabwurcs",
            "render",
            "--orientation",
            "top-down",
            "--compact",
            "--reducing-end",
            "Asn",
            "--repeat-brackets",
            "--floating-fragments",
            "--highlight-residue",
            "0",
            "--highlight-residue",
            "2",
            "--highlight-motif",
            "LewisX",
            "--motif-positions",
            "terminal,terminal,flexible",
            "--cleavage",
            "B1",
            "Gal(b1-4)GlcNAc",
        ])
        .unwrap();
        let Command::Render {
            highlights,
            drawing,
            ..
        } = cli.command
        else {
            panic!("expected render command");
        };
        assert_eq!(drawing.orientation, Some(Orientation::RootTop));
        let options = drawing.options();
        assert!(options.compact && options.repeat_brackets && options.floating_fragments);
        assert_eq!(options.reducing_end.as_deref(), Some("Asn"));
        assert_eq!(highlights.highlight_residues, [0, 2]);
        assert_eq!(highlights.cleavages, ["B1"]);
        assert!(
            Cli::try_parse_from(["crabwurcs", "render", "--orientation", "sideways", "Glc"])
                .is_err()
        );
    }

    #[test]
    fn default_drawing_flags_are_the_default_options() {
        let options = DrawingArgs::default().options();
        let defaults = RenderOptions::default();
        assert_eq!(options.orientation, defaults.orientation);
        assert_eq!(options.compact, defaults.compact);
        assert_eq!(options.reducing_end, defaults.reducing_end);
        assert_eq!(options.repeat_brackets, defaults.repeat_brackets);
        assert_eq!(options.floating_fragments, defaults.floating_fragments);
    }

    fn render_to_svg(
        highlight_motifs: Vec<String>,
        highlights: HighlightArgs,
        drawing: DrawingArgs,
    ) -> anyhow::Result<String> {
        static NEXT: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);
        let number = NEXT.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        let path = scratch("highlight").join(format!("{number}.svg"));
        run_render(
            Some("Fuc(a1-3)GlcNAc(b1-4)Gal(b1-4)GlcNAc".into()),
            false,
            Some(path.clone()),
            TextFormat::IupacCondensed,
            highlight_motifs,
            MotifTextFormat::Auto,
            &highlights,
            &drawing,
        )?;
        let svg = std::fs::read_to_string(&path).unwrap();
        std::fs::remove_file(path).unwrap();
        Ok(svg)
    }

    #[test]
    fn residues_and_positional_motifs_are_highlighted() {
        let residues = render_to_svg(
            Vec::new(),
            HighlightArgs {
                highlight_residues: vec![0, 1],
                ..HighlightArgs::default()
            },
            DrawingArgs::default(),
        )
        .unwrap();
        assert_eq!(residues.matches("data-node-index").count(), 4);
        assert_eq!(
            residues
                .matches("<g data-node-index=\"0\" class=\"motif-match\"")
                .count(),
            1
        );

        // Gal(b1-4)GlcNAc occurs once, and its Gal carries a GlcNAc.
        let terminal = HighlightArgs {
            motif_positions: vec!["terminal,flexible".into()],
            ..HighlightArgs::default()
        };
        let motif = vec!["Gal(b1-4)GlcNAc".to_owned()];
        let none = render_to_svg(motif.clone(), terminal, DrawingArgs::default()).unwrap();
        assert_eq!(none.matches("class=\"motif-match\"").count(), 0);
        let any = render_to_svg(
            motif.clone(),
            HighlightArgs::default(),
            DrawingArgs::default(),
        )
        .unwrap();
        assert_eq!(any.matches("class=\"motif-match\"").count(), 3);

        let too_many = HighlightArgs {
            motif_positions: vec!["flexible,flexible".into(), "flexible,flexible".into()],
            ..HighlightArgs::default()
        };
        assert!(render_to_svg(motif, too_many, DrawingArgs::default()).is_err());
    }

    #[test]
    fn a_cleavage_names_an_ion_of_the_glycan() {
        let marks = HighlightArgs {
            cleavages: vec!["B1\u{03B1}".into()],
            ..HighlightArgs::default()
        };
        assert!(render_to_svg(Vec::new(), marks, DrawingArgs::default()).is_err());
        let marks = HighlightArgs {
            cleavages: vec!["Y2".into()],
            ..HighlightArgs::default()
        };
        let svg = render_to_svg(Vec::new(), marks, DrawingArgs::default()).unwrap();
        assert!(svg.contains(">Y2</text>"));
    }

    #[test]
    fn the_analysis_commands_parse_and_run() {
        let man5 = "Man(a1-3)[Man(a1-3)[Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc";
        let cli = Cli::try_parse_from([
            "crabwurcs",
            "mass",
            "--reducing-end",
            "2-AB",
            "--adduct",
            "Na",
            "--adduct",
            "-H",
            "--charge",
            "2",
            man5,
        ])
        .unwrap();
        let Command::Mass {
            mass,
            adducts,
            charge,
            ..
        } = cli.command
        else {
            panic!("expected mass command");
        };
        assert_eq!(mass.reducing_end, ReducingEnd::TwoAb);
        assert_eq!(adducts, [Adduct::Sodium, Adduct::Deprotonated]);
        assert_eq!(charge, 2);

        assert!(
            run_mass(
                Some(man5.into()),
                false,
                TextFormat::Auto,
                &MassArgs::default(),
                &[Adduct::Sodium],
                1
            )
            .is_ok()
        );
        assert!(
            run_mass(
                Some("Foo".into()),
                false,
                TextFormat::Auto,
                &MassArgs::default(),
                &[],
                1
            )
            .is_err()
        );
        assert!(
            run_fragments(
                Some(man5.into()),
                false,
                TextFormat::Auto,
                &MassArgs::default(),
                Adduct::Proton
            )
            .is_ok()
        );
        let motif = |motif, positions: Option<&str>| {
            run_motif(
                motif,
                Some("Gal(b1-4)GlcNAc".into()),
                false,
                TextFormat::Auto,
                MotifTextFormat::Auto,
                positions.map(str::to_owned),
            )
        };
        assert!(motif("LacNAc", None).is_ok());
        assert!(motif("Gal(b1-4)GlcNAc", Some("terminal,reducing")).is_ok());
        assert!(motif("Gal(b1-4)GlcNAc", Some("internal,reducing")).is_err());
        assert!(motif("Gal(b1-4)GlcNAc", Some("terminal")).is_err());
        assert!(motif("LewisX", None).is_err());
        assert!(run_motifs(None, false, TextFormat::Auto).is_ok());
        assert!(run_motifs(Some(man5.into()), false, TextFormat::Auto).is_ok());
        assert!(run_same("Gal(b1-4)GlcNAc", "DGalpb1-4DGlcpNAc", TextFormat::Auto).is_ok());
        assert!(run_same("Gal(b1-4)GlcNAc", "Gal(b1-3)GlcNAc", TextFormat::Auto).is_err());
    }

    #[test]
    fn a_grid_reads_notations_and_captions_from_a_file() {
        let directory = scratch("grid");
        let list = directory.join("glycans.tsv");
        std::fs::write(&list, "Gal(b1-4)GlcNAc\tLacNAc\n\nGlc\n").unwrap();
        for name in ["grid.svg", "grid.png"] {
            let output = directory.join(name);
            run_grid(&list, &output, TextFormat::Auto, 2, &DrawingArgs::default()).unwrap();
            assert!(std::fs::metadata(&output).unwrap().len() > 100);
        }
        let svg = std::fs::read_to_string(directory.join("grid.svg")).unwrap();
        assert!(svg.contains(">LacNAc</text>"));
        assert_eq!(svg.matches("<svg ").count(), 3);
        assert!(
            run_grid(
                &list,
                &directory.join("grid.jpg"),
                TextFormat::Auto,
                2,
                &DrawingArgs::default()
            )
            .is_err()
        );
        std::fs::remove_dir_all(directory).unwrap();
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
            &HighlightArgs::default(),
            &DrawingArgs::default(),
        )
        .unwrap();
        run_render(
            Some("Glc".into()),
            false,
            Some(png_path.clone()),
            TextFormat::IupacCondensed,
            Vec::new(),
            MotifTextFormat::Auto,
            &HighlightArgs::default(),
            &DrawingArgs::default(),
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
