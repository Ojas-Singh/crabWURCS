# Python API

## Core types

- `Format`: `AUTO`, `WURCS`, `IUPAC_CONDENSED`, `IUPAC_EXTENDED`, `GLYCAM`,
  `SMILES`, `MOL`, and `SDF`.
- `ImageFormat`: `SVG` and `PNG`.
- `Glycan.parse(value, format=Format.AUTO)` parses a structure.
- `Glycan.to(format)` returns notation or molecular text.
- `Glycan.render(format="svg", highlight_motifs=None, motif_format="auto")`
  returns SVG text or PNG bytes.

## Convenience functions

```python
crabwurcs.detect_format(value)
crabwurcs.convert(value, to_format, from_format="auto")
crabwurcs.render_snfg(value, from_format="auto", image_format="svg")
crabwurcs.extract_pdb(contents, format="auto")
crabwurcs.extract_pdb_file(path)
```

## PDB results

`ExtractedGlycan` is immutable and contains a `Glycan`, optional attachment
site, and an immutable sequence of `PdbResidueReference`. Each reference has
`node_index`, `chain`, `sequence_number`, and `insertion_code`.

The package ships `py.typed` and complete public stubs. Native classes under
`crabwurcs._native` are private implementation details.
