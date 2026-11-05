# Python API

## Core types

- `Format`: `AUTO`, `WURCS`, `IUPAC_CONDENSED`, `IUPAC_EXTENDED`, `GLYCAM`,
  `SMILES`, `MOL`, and `SDF`.
- `ImageFormat`: `SVG` and `PNG`.
- `Glycan.parse(value, format=Format.AUTO)` parses a structure.
- `Glycan.to(format)` returns notation or molecular text.
- `Glycan.render(format="svg", highlight_motifs=None, motif_format="auto",
  ...)` returns SVG text or PNG bytes. See [SNFG rendering](rendering.md) for
  its other keyword arguments.
- `Glycan.mass()`, `.formula()`, `.composition()`, `.composition_string()` and
  `.fragments()` weigh a glycan; see [Mass and fragments](mass.md).
- `Glycan.contains(motif)`, `.count_motif(motif)` and `.find_motif(motif)` search
  for a motif, given as a notation, a library name or a `Glycan`. `positions`
  holds its residues' positions in reading order.
- `Glycan.motifs()` names the library motifs a glycan holds, and
  `crabwurcs.named_motifs()` lists the library.
- `Glycan.residues()` lists residues with their node index, name and parent.
- `==` and `hash()` compare structures, whatever notation they were read from.
  A stated anomer or position is more than an open one: `GlcNAc` is not
  `GlcNAc(b1-`. A non-reducing sugar such as sucrose is one structure from
  either residue, so two WURCS that differ only in that direction compare equal
  while each still writes as it was given.

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
