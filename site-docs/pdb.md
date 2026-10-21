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
  (GLYCAM keeps sulfates in `SO3` residues). A glycosidic bond through S or
  N is a WURCS bridge (`a1-b1*S*` for thiodigalactoside). A terminal carbon
  with one N, S or F (6-amino-6-deoxy, 6-deoxy-6-fluoro, a 6-thio link) is a
  hydroxymethyl position carrying that group (`6*F`). A glycan whose residues
  close a ring (a cyclodextrin) keeps its closing linkage.

Bonds are perceived from distances within each atom's valence (bonds inside a
residue first, then explicit records, then contacts between residues), so a
clash between residues does not break a ring. Explicit records are PDB
`CONECT` and `LINK` records and mmCIF `_struct_conn` rows of type `covale*`
(metal coordination, hydrogen bonds and bonds to a symmetry mate are left
out). Metal ions are never bonded. wwPDB Chemical Component Dictionary
coordinates of 49 components fix the conventions in the test suite.

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
mmCIF file, the first model and the `_struct_conn` rows are read. Of
alternative locations, the atoms without one and those of the first are
read. Gzipped files are not read; decompress them first.

Two HETATM residues that put three or more heavy atoms within 0.6 Å of each
other are one residue modelled twice (the α and the β anomer deposited as
two residues, as in 4D52 or 1AXZ); one of them is read: the one in the first
alternative location, then the one of higher occupancy, then the first in
the file.

Crystal items are not used. An mmCIF file whose space-group items disagree
(archive entries such as 31NP state an H-M symbol and a Hall symbol that
pdbtbx reads as different groups) is read again without its space-group and
cell items, as is one whose space group or cell is unknown (`?`); any other
parse error is reported.

Not read from coordinates, and reported as warnings: selenium substituents
(C2 of FSW is left undetermined), C-glycosides whose aglycone carbons belong
to the sugar component (ZDC and similar glycomimetics read as octulose-like
ketoses with the anomeric heteroatom missing; IUPAC condensed has no name
for either), and intramolecular anhydro bridges (AAL,
3,6-anhydro-L-galactose: an unrecognised group at C3 and C6, left out).

## Attachment site

`attachment_site` (`chain/residue/number`) is, in this order:

1. the residue an explicit record bonds to the anomeric atom of the glycan's
   root: the asparagine of an N-glycan, the serine or threonine of an
   O-glycan, the tryptophan of a C-mannoside, an aglycone ligand;
2. any other residue an explicit record bonds to a glycan residue, unless
   the bond graph did not take the record (a metal coordination, atoms more
   than 3 Å apart);
3. when the file records no covalent connection at all (no `LINK` record, no
   `covale` `_struct_conn` row, as in models written by modelling and
   simulation programs), the residue whose N, O, S or C atom is bonded, by
   distance, to the anomeric carbon of the root (within 1.75 Å; 2.0 Å for
   S), or whose carbon is bonded to the root's anomeric heteroatom. Water and
   GLYCAM's `ROH` (a free reducing end) are never an attachment; a file that
   records its connections is never second-guessed from distances;
4. otherwise the root residue itself.

Substituent residues (`SO3`, `MEX`, `ACX`, `PCX`) are never the attachment
site. Insertion codes are retained. PDB/mmCIF writing, conformer generation,
and coordinate optimization are intentionally outside the 0.3 scope.
