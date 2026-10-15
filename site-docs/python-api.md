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
crabwurcs.convert(value, to_format, from_format="auto", table="off")
crabwurcs.normalize_wurcs(value)
crabwurcs.render_snfg(value, from_format="auto", image_format="svg")
crabwurcs.extract_pdb(contents, format="auto", policy="coordinates")
crabwurcs.extract_pdb_file(path, policy="coordinates")
```

WURCS output is normalized by default. `normalize_wurcs` is strict and raises
`NormalizationError` when the input contains syntax that cannot be normalized
losslessly by the local model.

## PDB results

`ExtractedGlycan` is immutable and contains a `Glycan`, optional attachment
site, and an immutable sequence of `PdbResidueReference`. Each reference has
`node_index`, `chain`, `sequence_number`, and `insertion_code`.

Residues are read from their coordinates. `ExtractedGlycan.warnings` holds a
`PdbWarning` for every disagreement between a residue name and the
coordinates (`kind`, the residue's location and name, `named`, `coordinates`,
`positions`, `aspects`, and the `message` that `str()` returns), and
`ExtractedGlycan.readings` a `ResidueReading` per residue (`from_coordinates`,
`from_name`, `ring_atoms`, `backbone_atoms`, `centres`). `policy="names"` lets
names decide where they say something; `policy="strict"` raises `PdbError`
on a disagreement. See [PDB and mmCIF extraction](pdb.md).

The package ships `py.typed` and complete public stubs. Native classes under
`crabwurcs._native` are private implementation details.
