use crabwurcs::mol::{ChemFormat, MolError};
use crabwurcs::pdb::{ExtractedGlycanWithProvenance, NamePolicy, PdbReadOptions, PdbReading};
use crabwurcs::{
    Adduct, Format, GridOptions, HighlightSelection, MassOptions, MotifPosition, MotifQuery,
    Orientation, ReducingEnd, RenderOptions, ResidueGraph,
};
use pyo3::conversion::FromPyObjectOwned;
use pyo3::create_exception;
use pyo3::exceptions::PyException;
use pyo3::prelude::*;
use pyo3::types::PyDict;

create_exception!(_native, CrabWurcsError, PyException);
create_exception!(_native, ParseError, CrabWurcsError);
create_exception!(_native, ConversionError, CrabWurcsError);
create_exception!(_native, NormalizationError, ConversionError);
create_exception!(_native, NonConcreteError, ConversionError);
create_exception!(_native, PdbError, CrabWurcsError);
create_exception!(_native, RenderError, CrabWurcsError);
create_exception!(_native, MassError, CrabWurcsError);
create_exception!(_native, MotifError, CrabWurcsError);

fn parse_format(value: &str, allow_auto: bool) -> PyResult<Format> {
    match value.trim().to_ascii_lowercase().replace('_', "-").as_str() {
        "auto" if allow_auto => Ok(Format::Auto),
        "wurcs" => Ok(Format::Wurcs),
        "iupac" | "iupac-condensed" | "condensed" => Ok(Format::IupacCondensed),
        "iupac-extended" | "extended" => Ok(Format::IupacExtended),
        "glycam" => Ok(Format::Glycam),
        "smiles" => Ok(Format::Smiles),
        "mol" | "sdf" => Err(ParseError::new_err(format!(
            "{value} is a molecular format, not a notation format"
        ))),
        _ => Err(ParseError::new_err(format!("unknown format: {value}"))),
    }
}

fn molecular_format(value: &str) -> Option<ChemFormat> {
    match value.trim().to_ascii_lowercase().as_str() {
        "smiles" => Some(ChemFormat::Smiles),
        "mol" => Some(ChemFormat::Mol),
        "sdf" => Some(ChemFormat::Sdf),
        _ => None,
    }
}

fn map_mol_error(error: MolError) -> PyErr {
    if matches!(error, MolError::NonConcrete(_)) {
        NonConcreteError::new_err(error.to_string())
    } else {
        ConversionError::new_err(error.to_string())
    }
}

fn parse_graph(value: &str, format: &str) -> PyResult<ResidueGraph> {
    if matches!(format.trim().to_ascii_lowercase().as_str(), "mol" | "sdf") {
        return crabwurcs::mol::wurcs_from_molecule(
            value,
            molecular_format(format).expect("checked molecular format"),
        )
        .map_err(map_mol_error);
    }
    crabwurcs::parse_notation(value, parse_format(format, true)?)
        .map_err(|error| ParseError::new_err(error.to_string()))
}

fn graph_to(graph: &ResidueGraph, format: &str) -> PyResult<String> {
    if let Some(format) = molecular_format(format) {
        return crabwurcs::mol::molecule_from_wurcs(graph, format).map_err(map_mol_error);
    }
    crabwurcs::write_notation(graph, parse_format(format, false)?)
        .map_err(|error| ConversionError::new_err(error.to_string()))
}

fn map_mass_error(error: crabwurcs::MassError) -> PyErr {
    MassError::new_err(error.to_string())
}

fn map_render_error(error: crabwurcs::SnfgError) -> PyErr {
    RenderError::new_err(error.to_string())
}

fn parse_named<T: std::str::FromStr>(value: &str) -> PyResult<T>
where
    T::Err: std::fmt::Display,
{
    value
        .parse()
        .map_err(|error: T::Err| ParseError::new_err(error.to_string()))
}

/// A motif notation, or with the `auto` format a name of the motif library.
fn parse_motif(value: &str, format: &str, positions: Option<Vec<String>>) -> PyResult<MotifQuery> {
    let named = format
        .trim()
        .eq_ignore_ascii_case("auto")
        .then(|| crabwurcs::named_motif(value))
        .flatten();
    let query = match named {
        Some(motif) => motif.query(),
        None => MotifQuery::new(parse_graph(value, format)?),
    };
    let Some(positions) = positions else {
        return Ok(query);
    };
    let positions = positions
        .iter()
        .map(|position| parse_named::<MotifPosition>(position))
        .collect::<PyResult<Vec<_>>>()?;
    MotifQuery::with_positions(query.motif().clone(), &positions)
        .map_err(|error| MotifError::new_err(error.to_string()))
}

/// An entry of an options dictionary; a missing key and `None` mean unset.
fn entry<'py, T: FromPyObjectOwned<'py>>(
    options: Option<&Bound<'py, PyDict>>,
    key: &str,
) -> PyResult<Option<T>> {
    let Some(options) = options else {
        return Ok(None);
    };
    match options.get_item(key)? {
        Some(value) if !value.is_none() => Ok(Some(value.extract::<T>().map_err(Into::into)?)),
        _ => Ok(None),
    }
}

/// The drawing options of `Glycan.render`, from a dictionary of what was set.
fn drawing_options(drawing: Option<&Bound<'_, PyDict>>) -> PyResult<RenderOptions> {
    let defaults = RenderOptions::default();
    let scale = entry(drawing, "scale")?.unwrap_or(defaults.scale);
    if !(scale.is_finite() && scale > 0.0) {
        return Err(RenderError::new_err("scale must be a positive number"));
    }
    Ok(RenderOptions {
        orientation: match entry::<String>(drawing, "orientation")? {
            Some(orientation) => parse_named::<Orientation>(&orientation)?,
            None => Orientation::default(),
        },
        compact: entry(drawing, "compact")?.unwrap_or(defaults.compact),
        reducing_end: entry(drawing, "reducing_end")?,
        repeat_brackets: entry(drawing, "repeat_brackets")?.unwrap_or(defaults.repeat_brackets),
        floating_fragments: entry(drawing, "floating_fragments")?
            .unwrap_or(defaults.floating_fragments),
        show_labels: entry(drawing, "show_labels")?.unwrap_or(defaults.show_labels),
        show_linkages: entry(drawing, "show_linkages")?.unwrap_or(defaults.show_linkages),
        colour: entry(drawing, "colour")?.unwrap_or(defaults.colour),
        scale,
        ..defaults
    })
}

fn add_selection(union: &mut Option<HighlightSelection>, found: HighlightSelection) {
    let union = union.get_or_insert_with(HighlightSelection::default);
    union.node_indices.extend(found.node_indices);
    union.edge_indices.extend(found.edge_indices);
}

fn mass_options(reducing_end: &str, permethylated: bool) -> PyResult<MassOptions> {
    Ok(MassOptions {
        reducing_end: parse_named::<ReducingEnd>(reducing_end)?,
        permethylated,
    })
}

/// `(node_index, name, parent_index, is_terminal, is_recognized)`.
type ResidueRow = (usize, String, Option<usize>, bool, bool);

#[pyclass(name = "_Glycan", frozen, skip_from_py_object)]
#[derive(Clone)]
struct PyGlycan {
    graph: ResidueGraph,
    source_value: Option<String>,
    source_format: Option<String>,
}

#[pymethods]
impl PyGlycan {
    #[staticmethod]
    fn parse(value: &str, format: &str) -> PyResult<Self> {
        Ok(Self {
            graph: parse_graph(value, format)?,
            source_value: Some(value.trim().to_owned()),
            source_format: Some(format.to_owned()),
        })
    }

    fn to_format(&self, format: &str) -> PyResult<String> {
        graph_to(&self.graph, format)
    }

    #[pyo3(signature = (highlight_motifs=None, motif_format="auto", drawing=None, marks=None))]
    fn render_svg(
        &self,
        highlight_motifs: Option<Vec<String>>,
        motif_format: &str,
        drawing: Option<&Bound<'_, PyDict>>,
        marks: Option<&Bound<'_, PyDict>>,
    ) -> PyResult<String> {
        let (opts, selection) =
            self.render_request(highlight_motifs, motif_format, drawing, marks)?;
        match selection {
            Some(selection) => crabwurcs::render_svg_with_selection(&self.graph, &selection, &opts),
            None => crabwurcs::render_svg_with_options(&self.graph, &opts),
        }
        .map_err(map_render_error)
    }

    #[pyo3(signature = (highlight_motifs=None, motif_format="auto", drawing=None, marks=None))]
    fn render_png(
        &self,
        highlight_motifs: Option<Vec<String>>,
        motif_format: &str,
        drawing: Option<&Bound<'_, PyDict>>,
        marks: Option<&Bound<'_, PyDict>>,
    ) -> PyResult<Vec<u8>> {
        let (opts, selection) =
            self.render_request(highlight_motifs, motif_format, drawing, marks)?;
        match selection {
            Some(selection) => crabwurcs::render_png_with_selection(&self.graph, &selection, &opts),
            None => crabwurcs::render_png_with_options(&self.graph, &opts),
        }
        .map_err(map_render_error)
    }

    fn structure_key(&self) -> Option<String> {
        crabwurcs::structure_key(&self.graph)
    }

    #[pyo3(signature = (motif, format="auto", positions=None))]
    fn find_motif(
        &self,
        motif: &str,
        format: &str,
        positions: Option<Vec<String>>,
    ) -> PyResult<Vec<(Vec<usize>, Vec<usize>)>> {
        let matches = parse_motif(motif, format, positions)?
            .find(&self.graph)
            .map_err(|error| MotifError::new_err(error.to_string()))?;
        Ok(matches
            .into_iter()
            .map(|found| {
                (
                    found.node_indices.into_iter().collect(),
                    found.edge_indices.into_iter().collect(),
                )
            })
            .collect())
    }

    fn annotate_motifs(&self) -> Vec<&'static str> {
        crabwurcs::annotate_motifs(&self.graph)
    }

    fn composition(&self) -> Vec<(String, u32)> {
        crabwurcs::glycan_composition(&self.graph)
            .classes()
            .to_vec()
    }

    /// `(formula, monoisotopic, average)` of the neutral glycan.
    fn mass(&self, reducing_end: &str, permethylated: bool) -> PyResult<(String, f64, f64)> {
        let options = mass_options(reducing_end, permethylated)?;
        let mass = crabwurcs::glycan_mass(&self.graph, &options).map_err(map_mass_error)?;
        Ok((mass.formula.to_string(), mass.monoisotopic, mass.average))
    }

    /// `(monoisotopic, average)` m/z with `count` adducts.
    fn mz(
        &self,
        adduct: &str,
        count: u32,
        reducing_end: &str,
        permethylated: bool,
    ) -> PyResult<(f64, f64)> {
        let options = mass_options(reducing_end, permethylated)?;
        let adduct = parse_named::<Adduct>(adduct)?;
        let mass = crabwurcs::glycan_mass(&self.graph, &options).map_err(map_mass_error)?;
        Ok((
            mass.mz(adduct, count).map_err(map_mass_error)?,
            mass.average_mz(adduct, count).map_err(map_mass_error)?,
        ))
    }

    fn fragments(
        &self,
        adduct: &str,
        reducing_end: &str,
        permethylated: bool,
    ) -> PyResult<Vec<PyFragment>> {
        let options = mass_options(reducing_end, permethylated)?;
        let adduct = parse_named::<Adduct>(adduct)?;
        Ok(crabwurcs::glycan_fragments(&self.graph, &options, adduct)
            .map_err(map_mass_error)?
            .into_iter()
            .map(|fragment| PyFragment {
                kind: fragment.kind.as_str(),
                label: fragment.label,
                edge_index: fragment.edge,
                residues: fragment.residues.into_iter().collect(),
                formula: fragment.formula.to_string(),
                monoisotopic: fragment.monoisotopic,
                mz: fragment.mz,
            })
            .collect())
    }

    /// The residues, by node index.
    fn residues(&self) -> PyResult<Vec<ResidueRow>> {
        let inner = self.graph.inner();
        inner
            .node_indices()
            .map(|node| {
                let residue = &inner[node];
                let symbol = crabwurcs::snfg::symbol_for(residue).map_err(map_render_error)?;
                let recognized = !matches!(
                    crabwurcs::classify_residue(residue),
                    None | Some(crabwurcs::ResidueKind::Assigned)
                );
                Ok((
                    node.index(),
                    residue
                        .display_name
                        .clone()
                        .unwrap_or_else(|| symbol.label.clone()),
                    self.graph.parent(node).map(|parent| parent.index()),
                    self.graph.children(node).is_empty(),
                    recognized,
                ))
            })
            .collect()
    }

    fn __repr__(&self) -> String {
        let wurcs = crabwurcs::core::write_wurcs(&self.graph)
            .unwrap_or_else(|_| "<unrepresentable>".to_owned());
        format!("_Glycan({wurcs:?})")
    }
}

impl PyGlycan {
    /// The options and the selection one render call asks for. The selection
    /// is `None` when nothing was asked to be highlighted.
    fn render_request(
        &self,
        highlight_motifs: Option<Vec<String>>,
        motif_format: &str,
        drawing: Option<&Bound<'_, PyDict>>,
        marks: Option<&Bound<'_, PyDict>>,
    ) -> PyResult<(RenderOptions, Option<HighlightSelection>)> {
        let mut opts = RenderOptions {
            source_notation: self.source_value.as_ref().map(|value| {
                crabwurcs::SourceNotation::new(
                    self.source_format.as_deref().unwrap_or("auto"),
                    value,
                )
            }),
            ..drawing_options(drawing)?
        };

        let motifs = highlight_motifs.unwrap_or_default();
        let positions = entry::<Vec<Option<Vec<String>>>>(marks, "positions")?.unwrap_or_default();
        if positions.len() > motifs.len() {
            return Err(MotifError::new_err(
                "more highlight positions than highlight motifs",
            ));
        }
        let mut selection = None;
        for (index, motif) in motifs.iter().enumerate() {
            let positions = positions.get(index).cloned().flatten();
            let matches = parse_motif(motif, motif_format, positions)?
                .find(&self.graph)
                .map_err(|error| RenderError::new_err(error.to_string()))?;
            add_selection(&mut selection, HighlightSelection::from_matches(&matches));
        }
        if let Some(residues) = entry::<Vec<usize>>(marks, "residues")? {
            let found = HighlightSelection::from_residues(&self.graph, residues)
                .map_err(map_render_error)?;
            add_selection(&mut selection, found);
        }
        if let Some(cleavages) = entry::<Vec<String>>(marks, "cleavages")? {
            let fragments =
                crabwurcs::glycan_fragments(&self.graph, &MassOptions::default(), Adduct::Proton)
                    .map_err(map_mass_error)?;
            for ion in cleavages {
                let fragment = fragments
                    .iter()
                    .find(|fragment| fragment.label == ion)
                    .ok_or_else(|| {
                        RenderError::new_err(format!("this glycan has no fragment ion {ion}"))
                    })?;
                opts.cleavages.push(crabwurcs::cleavage_mark(fragment));
            }
        }
        Ok((opts, selection))
    }
}

#[pyclass(name = "_Fragment", frozen, get_all, skip_from_py_object)]
#[derive(Clone)]
struct PyFragment {
    kind: &'static str,
    label: String,
    edge_index: usize,
    residues: Vec<usize>,
    formula: String,
    monoisotopic: f64,
    mz: f64,
}

#[pyclass(name = "_PdbResidueReference", frozen, get_all, skip_from_py_object)]
#[derive(Clone)]
struct PyPdbResidueReference {
    node_index: usize,
    chain: String,
    sequence_number: isize,
    insertion_code: Option<String>,
}

#[pyclass(name = "_PdbWarning", frozen, get_all, skip_from_py_object)]
#[derive(Clone)]
struct PyPdbWarning {
    node_index: usize,
    chain: String,
    sequence_number: isize,
    insertion_code: Option<String>,
    residue_name: String,
    kind: String,
    named: Option<String>,
    coordinates: Option<String>,
    positions: Vec<u8>,
    aspects: Vec<String>,
    message: String,
}

#[pyclass(name = "_StereoCentre", frozen, get_all, skip_from_py_object)]
#[derive(Clone)]
struct PyStereoCentre {
    position: u8,
    atom: String,
    code: String,
    volume: f64,
}

#[pyclass(name = "_ResidueReading", frozen, get_all, skip_from_py_object)]
#[derive(Clone)]
struct PyResidueReading {
    node_index: usize,
    chain: String,
    sequence_number: isize,
    insertion_code: Option<String>,
    residue_name: String,
    from_coordinates: Option<String>,
    from_name: Option<String>,
    ring_atoms: Vec<String>,
    backbone_atoms: Vec<String>,
    centres: Vec<PyStereoCentre>,
}

#[pyclass(name = "_ExtractedGlycan", frozen, get_all, skip_from_py_object)]
#[derive(Clone)]
struct PyExtractedGlycan {
    glycan: PyGlycan,
    attachment_site: Option<String>,
    residues: Vec<PyPdbResidueReference>,
    warnings: Vec<PyPdbWarning>,
    readings: Vec<PyResidueReading>,
}

fn extracted_to_python(
    index: usize,
    value: ExtractedGlycanWithProvenance,
    reading: &PdbReading,
) -> PyExtractedGlycan {
    PyExtractedGlycan {
        glycan: PyGlycan {
            graph: value.graph,
            source_value: None,
            source_format: None,
        },
        attachment_site: value.attachment_site,
        residues: value
            .residues
            .into_iter()
            .map(|residue| PyPdbResidueReference {
                node_index: residue.node_index,
                chain: residue.chain,
                sequence_number: residue.sequence_number,
                insertion_code: residue.insertion_code,
            })
            .collect(),
        warnings: reading
            .warnings
            .iter()
            .filter(|warning| warning.glycan_index == index)
            .map(|warning| PyPdbWarning {
                node_index: warning.node_index,
                chain: warning.chain.clone(),
                sequence_number: warning.sequence_number,
                insertion_code: warning.insertion_code.clone(),
                residue_name: warning.residue_name.clone(),
                kind: warning.kind.as_str().to_string(),
                named: warning.named.clone(),
                coordinates: warning.coordinates.clone(),
                positions: warning.positions.clone(),
                aspects: warning.aspects.clone(),
                message: warning.to_string(),
            })
            .collect(),
        readings: reading
            .residues
            .iter()
            .filter(|residue| residue.glycan_index == index)
            .map(|residue| PyResidueReading {
                node_index: residue.node_index,
                chain: residue.chain.clone(),
                sequence_number: residue.sequence_number,
                insertion_code: residue.insertion_code.clone(),
                residue_name: residue.residue_name.clone(),
                from_coordinates: residue.from_coordinates.clone(),
                from_name: residue.from_name.clone(),
                ring_atoms: residue.ring_atoms.clone(),
                backbone_atoms: residue.backbone_atoms.clone(),
                centres: residue
                    .centres
                    .iter()
                    .map(|centre| PyStereoCentre {
                        position: centre.position,
                        atom: centre.atom.clone(),
                        code: centre.code.to_string(),
                        volume: centre.volume,
                    })
                    .collect(),
            })
            .collect(),
    }
}

fn reading_to_python(mut reading: PdbReading) -> Vec<PyExtractedGlycan> {
    let glycans = std::mem::take(&mut reading.glycans);
    glycans
        .into_iter()
        .enumerate()
        .map(|(index, glycan)| extracted_to_python(index, glycan, &reading))
        .collect()
}

fn parse_policy(policy: &str) -> PyResult<PdbReadOptions> {
    policy
        .parse::<NamePolicy>()
        .map(PdbReadOptions::with_policy)
        .map_err(ParseError::new_err)
}

#[pyfunction]
#[pyo3(signature = (contents, format="auto", policy="coordinates"))]
fn extract_pdb(contents: &str, format: &str, policy: &str) -> PyResult<Vec<PyExtractedGlycan>> {
    let is_mmcif = match format.trim().to_ascii_lowercase().as_str() {
        "auto" => contents.trim_start().starts_with("data_"),
        "pdb" => false,
        "mmcif" | "cif" => true,
        _ => {
            return Err(ParseError::new_err(format!(
                "unknown structure format: {format}"
            )));
        }
    };
    let options = parse_policy(policy)?;
    crabwurcs::read_glycans_from_str(contents, is_mmcif, &options)
        .map(reading_to_python)
        .map_err(|error| PdbError::new_err(error.to_string()))
}

#[pyfunction]
#[pyo3(signature = (path, policy="coordinates"))]
fn extract_pdb_file(path: &str, policy: &str) -> PyResult<Vec<PyExtractedGlycan>> {
    let options = parse_policy(policy)?;
    crabwurcs::read_glycans_from_file(std::path::Path::new(path), &options)
        .map(reading_to_python)
        .map_err(|error| PdbError::new_err(error.to_string()))
}

#[pyfunction]
fn detect_format(value: &str) -> &'static str {
    match crabwurcs::detect_format(value) {
        Format::Auto => "auto",
        Format::Wurcs => "wurcs",
        Format::IupacCondensed => "iupac-condensed",
        Format::IupacExtended => "iupac-extended",
        Format::Glycam => "glycam",
        Format::Smiles => "smiles",
    }
}

/// Convert between notations with the bundled GlycoShape table used as
/// `table` says (`off`, `fallback` or `check`).
#[pyfunction]
fn convert_with_table(
    value: &str,
    to_format: &str,
    from_format: &str,
    table: &str,
) -> PyResult<String> {
    let mode = match table.trim().to_ascii_lowercase().as_str() {
        "off" => crabwurcs::TableMode::Off,
        "fallback" => crabwurcs::TableMode::Fallback,
        "check" => crabwurcs::TableMode::Check,
        other => {
            return Err(ParseError::new_err(format!(
                "unknown table mode: {other} (expected off, fallback or check)"
            )));
        }
    };
    crabwurcs::convert_with_table(
        value,
        parse_format(from_format, true)?,
        parse_format(to_format, false)?,
        mode,
    )
    .map_err(|error| ConversionError::new_err(error.to_string()))
}

/// `(name, aliases, iupac, positions)` of every library motif.
#[pyfunction]
fn named_motifs() -> Vec<(
    &'static str,
    Vec<&'static str>,
    &'static str,
    Vec<&'static str>,
)> {
    crabwurcs::named_motifs()
        .iter()
        .map(|motif| {
            (
                motif.name,
                motif.aliases.to_vec(),
                motif.iupac,
                motif
                    .positions
                    .iter()
                    .map(|position| position.as_str())
                    .collect(),
            )
        })
        .collect()
}

#[pyfunction]
fn svg_to_png(svg: &str) -> PyResult<Vec<u8>> {
    crabwurcs::svg_to_png(svg).map_err(map_render_error)
}

fn grid_options(
    columns: usize,
    labels: Option<Vec<String>>,
    drawing: Option<&Bound<'_, PyDict>>,
) -> PyResult<GridOptions> {
    Ok(GridOptions {
        columns,
        labels: labels.unwrap_or_default(),
        render: drawing_options(drawing)?,
    })
}

#[pyfunction]
#[pyo3(signature = (glycans, columns=4, labels=None, drawing=None))]
fn render_grid_svg(
    glycans: Vec<PyRef<'_, PyGlycan>>,
    columns: usize,
    labels: Option<Vec<String>>,
    drawing: Option<&Bound<'_, PyDict>>,
) -> PyResult<String> {
    let graphs = glycans
        .iter()
        .map(|glycan| glycan.graph.clone())
        .collect::<Vec<_>>();
    crabwurcs::render_grid_svg(&graphs, &grid_options(columns, labels, drawing)?)
        .map_err(map_render_error)
}

#[pyfunction]
#[pyo3(signature = (glycans, columns=4, labels=None, drawing=None))]
fn render_grid_png(
    glycans: Vec<PyRef<'_, PyGlycan>>,
    columns: usize,
    labels: Option<Vec<String>>,
    drawing: Option<&Bound<'_, PyDict>>,
) -> PyResult<Vec<u8>> {
    let graphs = glycans
        .iter()
        .map(|glycan| glycan.graph.clone())
        .collect::<Vec<_>>();
    crabwurcs::render_grid_png(&graphs, &grid_options(columns, labels, drawing)?)
        .map_err(map_render_error)
}

#[pyfunction]
fn normalize_wurcs(value: &str) -> PyResult<String> {
    crabwurcs::normalize_wurcs(value)
        .map_err(|error| NormalizationError::new_err(error.to_string()))
}

#[pymodule]
#[pyo3(name = "_native")]
fn crabwurcs_native(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyGlycan>()?;
    module.add_class::<PyPdbResidueReference>()?;
    module.add_class::<PyExtractedGlycan>()?;
    module.add_class::<PyPdbWarning>()?;
    module.add_class::<PyResidueReading>()?;
    module.add_class::<PyStereoCentre>()?;
    module.add_class::<PyFragment>()?;
    module.add_function(wrap_pyfunction!(named_motifs, module)?)?;
    module.add_function(wrap_pyfunction!(svg_to_png, module)?)?;
    module.add_function(wrap_pyfunction!(render_grid_svg, module)?)?;
    module.add_function(wrap_pyfunction!(render_grid_png, module)?)?;
    module.add_function(wrap_pyfunction!(extract_pdb, module)?)?;
    module.add_function(wrap_pyfunction!(extract_pdb_file, module)?)?;
    module.add_function(wrap_pyfunction!(detect_format, module)?)?;
    module.add_function(wrap_pyfunction!(normalize_wurcs, module)?)?;
    module.add_function(wrap_pyfunction!(convert_with_table, module)?)?;
    module.add("CrabWurcsError", module.py().get_type::<CrabWurcsError>())?;
    module.add("ParseError", module.py().get_type::<ParseError>())?;
    module.add("ConversionError", module.py().get_type::<ConversionError>())?;
    module.add(
        "NormalizationError",
        module.py().get_type::<NormalizationError>(),
    )?;
    module.add(
        "NonConcreteError",
        module.py().get_type::<NonConcreteError>(),
    )?;
    module.add("PdbError", module.py().get_type::<PdbError>())?;
    module.add("RenderError", module.py().get_type::<RenderError>())?;
    module.add("MassError", module.py().get_type::<MassError>())?;
    module.add("MotifError", module.py().get_type::<MotifError>())?;
    module.add("__version__", env!("CARGO_PKG_VERSION"))?;
    Ok(())
}
