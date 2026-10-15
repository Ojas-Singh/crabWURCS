use crabwurcs::mol::{ChemFormat, MolError};
use crabwurcs::pdb::{ExtractedGlycanWithProvenance, NamePolicy, PdbReadOptions, PdbReading};
use crabwurcs::{Format, ResidueGraph};
use pyo3::create_exception;
use pyo3::exceptions::PyException;
use pyo3::prelude::*;

create_exception!(_native, CrabWurcsError, PyException);
create_exception!(_native, ParseError, CrabWurcsError);
create_exception!(_native, ConversionError, CrabWurcsError);
create_exception!(_native, NormalizationError, ConversionError);
create_exception!(_native, NonConcreteError, ConversionError);
create_exception!(_native, PdbError, CrabWurcsError);
create_exception!(_native, RenderError, CrabWurcsError);

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

    #[pyo3(signature = (highlight_motifs=None, motif_format="auto"))]
    fn render_svg(
        &self,
        highlight_motifs: Option<Vec<String>>,
        motif_format: &str,
    ) -> PyResult<String> {
        let opts = crabwurcs::RenderOptions {
            source_notation: self.source_value.as_ref().map(|value| {
                crabwurcs::SourceNotation::new(
                    self.source_format.as_deref().unwrap_or("auto"),
                    value,
                )
            }),
            ..crabwurcs::RenderOptions::default()
        };
        let motifs = highlight_motifs
            .unwrap_or_default()
            .into_iter()
            .map(|motif| parse_graph(&motif, motif_format))
            .collect::<PyResult<Vec<_>>>()?;
        if motifs.is_empty() {
            crabwurcs::render_svg_with_options(&self.graph, &opts)
        } else {
            crabwurcs::render_svg_with_motifs(&self.graph, &motifs, &opts)
        }
        .map_err(|error| RenderError::new_err(error.to_string()))
    }

    #[pyo3(signature = (highlight_motifs=None, motif_format="auto"))]
    fn render_png(
        &self,
        highlight_motifs: Option<Vec<String>>,
        motif_format: &str,
    ) -> PyResult<Vec<u8>> {
        let opts = crabwurcs::RenderOptions {
            source_notation: self.source_value.as_ref().map(|value| {
                crabwurcs::SourceNotation::new(
                    self.source_format.as_deref().unwrap_or("auto"),
                    value,
                )
            }),
            ..crabwurcs::RenderOptions::default()
        };
        let motifs = highlight_motifs
            .unwrap_or_default()
            .into_iter()
            .map(|motif| parse_graph(&motif, motif_format))
            .collect::<PyResult<Vec<_>>>()?;
        if motifs.is_empty() {
            crabwurcs::render_png_with_options(&self.graph, &opts)
        } else {
            crabwurcs::render_png_with_motifs(&self.graph, &motifs, &opts)
        }
        .map_err(|error| RenderError::new_err(error.to_string()))
    }

    fn __repr__(&self) -> String {
        let wurcs = crabwurcs::core::write_wurcs(&self.graph)
            .unwrap_or_else(|_| "<unrepresentable>".to_owned());
        format!("_Glycan({wurcs:?})")
    }
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
    module.add("__version__", env!("CARGO_PKG_VERSION"))?;
    Ok(())
}
