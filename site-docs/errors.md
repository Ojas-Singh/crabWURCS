# Errors and limitations

Python errors derive from `CrabWurcsError`:

- `ParseError` for invalid or unrecognized input.
- `ConversionError` for unsupported output conversion.
- `NonConcreteError` when a graph does not denote one molecule.
- `PdbError` for structure parsing and extraction failures.
- `RenderError` for SNFG/SVG/PNG failures.

Generic residues, compositions, undefined linkage positions, probability
ensembles, undefined modifications, and variable repeats are not converted to
an arbitrary concrete molecule. Exact-count repeats and finite defined cyclic
graphs may be materialized.

Format conversion cannot restore information absent from the source. PDB
coordinates without adequate component identity, bonds, or stereochemical
geometry may yield conservative unknown values or no glycan rather than a
chemically invented answer.
