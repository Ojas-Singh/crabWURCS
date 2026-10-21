# Formats and residue coverage

| Format | Read | Write | Automatic input detection |
|---|---:|---:|---:|
| WURCS 2.0 | yes | yes | yes |
| IUPAC condensed | yes | yes | yes |
| IUPAC extended | yes | yes | yes |
| GLYCAM | yes | yes | yes |
| stereochemical SMILES | yes | yes | yes for ordinary SMILES |
| V3000 MOL | yes | yes | no; specify `mol` |
| SDF | yes | yes | no; specify `sdf` |
| PDB/mmCIF | extraction | no | by file/text header |

The authoritative registry contains 87 entries. The 74 non-generic entries
have defined backbone chemistry and are release-gated through notation,
molecular, and PDB/mmCIF recognition tests. The complete named list, WURCS
UniqueRES values, aliases, and generic flags is maintained in the
[generated monosaccharide reference](https://github.com/Ojas-Singh/crabWURCS/blob/main/crabwurcs/docs/supported-monosaccharides.md).

The 13 uncertainty-preserving entries are `Hex`, `HexNAc`, `HexN`, `HexA`,
`dHex`, `dHexNAc`, `ddHex`, `Pen`, `NulO`, `Sia`, `ddNulO`, `Unknown`, and
`Assigned`. They can be parsed and rendered, but molecular export fails when
it would require choosing missing stereochemistry.

## IUPAC condensed names

A residue is written with its registry name (`Gal`, `L-Gal`, `Neu5Ac`) and
the substituents the notation spells: `NAc`, `NGc`, `NS`, `Me`, `Ac`, `S`
(O-sulfate), `P` (phosphate, `Man6P`), `PCho`. A residue the registry does
not name gets a name only when that name reads back as the same residue;
otherwise writing IUPAC condensed (and IUPAC extended and GLYCAM, which are
written through it) is a `ConversionError` such as `residue
ax221m-1a_1-5_1*OC (no IUPAC condensed name) cannot be represented in IUPAC
without losing information`. Earlier releases wrote a guess (`Res4`,
`Res41Me`, `GalA` for a heptose, `Fuc` for a 3-deoxy-GlcNAc, `Man` for a
mannose whose C6 was undetermined) that read back as another residue or not
at all.

No IUPAC condensed name, by design:

- a residue with an undetermined centre (`x`, for example read from a
  distorted structure), unless every centre is undetermined (`Hex`, `HexA`);
- skeletons outside the registry: most deoxy sugars (`3dGlc`), heptoses other
  than the L-/D-glycero-D-manno-heptoses, octulosonic acids other than Kdo,
  C-glycosides (ZDC, a C-fucoside, reads from coordinates as an
  octulose-like `ma1121dh-2x_2-6`);
- substituents the notation does not spell: fluoro (`F`), an amino group the
  name does not carry (`4*N`), a thiol (`1*S`), selenium, other WURCS MAPs.

Thioglycosidic and N-glycosidic bonds between residues are bridges:
`Gal(b1-SH-1)Gal(b1-` is WURCS `a1-b1*S*`. A cyclic glycan is written with
its closing position first: `4)Glc(a1-4)Glc(a1-4)…Glc(a1-`.

Conversions are computed from the parsed glycan; the bundled table of
GlycoShape notations is consulted only on request (`convert_with_table`,
`convert(..., table=...)`, `--table`). The registry-derived
molecular index and de-novo recognizer handle concrete residues without a
corpus hit; they do not require network access at runtime.
