# PDB and mmCIF extraction

```python
results = crabwurcs.extract_pdb_file("model.cif")
for result in results:
    print(result.glycan.to("wurcs"))
    for residue in result.residues:
        print(residue.node_index, residue.chain, residue.sequence_number)
```

Recognition uses, in order, the pinned wwPDB CCD component table, concrete
registry names, GLYCAM residue decoding, and atom/coordinate-graph fallback
for renamed components. Explicit connectivity records are authoritative;
coordinate inference is conservative and does not invent inter-residue chains.

PDB text is read from the records the extractor uses: the `ATOM`/`HETATM`
records of the first model and the `CONECT` and `LINK` records. Free-text
records (such as banner `REMARK` lines) and later models are skipped.

Attachment sites and insertion codes are retained. PDB/mmCIF writing,
conformer generation, and coordinate optimization are intentionally outside
the 0.3 scope.
