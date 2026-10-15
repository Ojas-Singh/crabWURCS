# PDB and mmCIF extraction

```python
results = crabwurcs.extract_pdb_file("model.cif")
for result in results:
    print(result.glycan.to("wurcs"))
    for residue in result.residues:
        print(residue.node_index, residue.chain, residue.sequence_number)
    for warning in result.warnings:
        print(warning)  # residue 2uA A:4: named α-L-IdopA, coordinates say α-L-GulpA2S (C2)
```

## Coordinates decide

Every residue is searched for a monosaccharide ring (five or six atoms, one
oxygen) in the bond graph, whatever the residue is called. From the
coordinates alone:

- the carbon numbering comes from the ring graph: the ring carbon next to the
  ring oxygen that carries an exocyclic heteroatom is the anomeric carbon, C1
  precedes it in a ketose, and the numbering walks the ring away from the ring
  oxygen and on along the exocyclic carbon chain;
- every stereocentre's configuration is the sign of the signed volume of its
  three heavy-atom substituents, ordered by the ring graph alone (previous
  carbon, next carbon, exocyclic heteroatom; the ring oxygen takes the place
  of the previous carbon at the anomeric carbon). That gives D/L, the parent
  sugar, deoxy positions and uronic acids;
- the anomer compares the anomeric centre with the anomeric reference atom
  (IUPAC 2-Carb-6.2);
- substituents (N-acetyl, N-glycolyl, N-sulfate, amino, O-sulfate, phosphate,
  phosphocholine, O-methyl, O-acetyl, fluoro) and glycosidic linkages come
  from what the ring atoms are bonded to, in the same residue or in another
  (GLYCAM keeps sulfates in `SO3` residues).

Bonds are perceived from distances within each atom's valence (bonds inside a
residue first, then explicit `CONECT`/`LINK` records, then contacts between
residues), so a clash between residues does not break a ring. Metal ions are
never bonded. wwPDB Chemical Component Dictionary coordinates of 49 components
fix the conventions in the test suite.

## Names are hints

A residue name fills in only what the coordinates cannot show: atoms missing
(a crystal structure without hydrogens: an absent O4 or N-acetyl group), a
centre too close to planar to read (normalised volume below 0.15), or a
residue whose ring atoms are missing (read from its name, as before). Names
come from the pinned wwPDB CCD component table, concrete registry names and
GLYCAM residue codes (upper case D, lower case L).

Every disagreement between a name and the coordinates, and every part the
coordinates could not decide, is a warning. `policy` chooses what decides:

| policy | on a disagreement |
|---|---|
| `coordinates` (default) | the coordinates decide; a warning |
| `names` | the name decides where it says something; a warning |
| `strict` | `PdbError` listing the disagreements |

Warning kinds: `name-contradicts-coordinates`, `missing-atoms`,
`distorted-geometry`, `unreadable-geometry`, `no-ring`, `unrecognized-group`.
`ExtractedGlycan.readings` gives, per residue, the WURCS UniqueRES read from
the coordinates and the one its name describes, its ring and backbone atoms,
and the Fischer side and normalised volume of every centre.

In Rust, `read_glycans_from_file(path, &PdbReadOptions::with_policy(...))`
returns a `PdbReading` with `glycans`, `warnings` and `residues`; the
`extract_glycans_*` functions read with the default policy.

## Records read

PDB text is read from the records the extractor uses: the `ATOM`/`HETATM`
records of the first model and the `CONECT` and `LINK` records. Free-text
records (such as banner `REMARK` lines) and later models are skipped. Of an
mmCIF file, the first model is read. Of alternative locations, the atoms
without one and those of the first are read.

The attachment site is the residue an explicit `CONECT` or `LINK` record bonds
the glycan to (the bond at the root's anomeric atom first; never a substituent
residue such as `SO3`), otherwise the root residue. Insertion codes are
retained. PDB/mmCIF writing, conformer generation, and coordinate optimization
are intentionally outside the 0.3 scope.
