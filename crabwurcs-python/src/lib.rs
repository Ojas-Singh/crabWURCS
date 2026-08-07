use crabwurcs::mol::{ChemFormat, MolError};
use crabwurcs::pdb::ExtractedGlycanWithProvenance;
use crabwurcs::{Format, ResidueGraph};
use pyo3::create_exception;
use pyo3::exceptions::PyException;
use pyo3::prelude::*;

create_exception!(_native, CrabWurcsError, PyException);
create_exception!(_native, ParseError, CrabWurcsError);
create_exception!(_native, ConversionError, CrabWurcsError);
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

#[pyclass(name = "_ExtractedGlycan", frozen, get_all, skip_from_py_object)]
#[derive(Clone)]
struct PyExtractedGlycan {
    glycan: PyGlycan,
    attachment_site: Option<String>,
    residues: Vec<PyPdbResidueReference>,
}

fn extracted_to_python(value: ExtractedGlycanWithProvenance) -> PyExtractedGlycan {
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
    }
}

#[pyfunction]
#[pyo3(signature = (contents, format="auto"))]
fn extract_pdb(contents: &str, format: &str) -> PyResult<Vec<PyExtractedGlycan>> {
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
    crabwurcs::extract_glycans_with_provenance_from_str(contents, is_mmcif)
        .map(|values| values.into_iter().map(extracted_to_python).collect())
        .map_err(|error| PdbError::new_err(error.to_string()))
}

#[pyfunction]
fn extract_pdb_file(path: &str) -> PyResult<Vec<PyExtractedGlycan>> {
    crabwurcs::extract_glycans_with_provenance_from_file(std::path::Path::new(path))
        .map(|values| values.into_iter().map(extracted_to_python).collect())
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

#[pymodule]
#[pyo3(name = "_native")]
fn crabwurcs_native(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyGlycan>()?;
    module.add_class::<PyPdbResidueReference>()?;
    module.add_class::<PyExtractedGlycan>()?;
    module.add_function(wrap_pyfunction!(extract_pdb, module)?)?;
    module.add_function(wrap_pyfunction!(extract_pdb_file, module)?)?;
    module.add_function(wrap_pyfunction!(detect_format, module)?)?;
    module.add("CrabWurcsError", module.py().get_type::<CrabWurcsError>())?;
    module.add("ParseError", module.py().get_type::<ParseError>())?;
    module.add("ConversionError", module.py().get_type::<ConversionError>())?;
    module.add(
        "NonConcreteError",
        module.py().get_type::<NonConcreteError>(),
    )?;
    module.add("PdbError", module.py().get_type::<PdbError>())?;
    module.add("RenderError", module.py().get_type::<RenderError>())?;
    module.add("__version__", env!("CARGO_PKG_VERSION"))?;
    Ok(())
}
