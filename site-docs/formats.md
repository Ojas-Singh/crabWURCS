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

Bundled corpus tables accelerate exact known conversions. The registry-derived
molecular index and de-novo recognizer handle concrete residues without a
corpus hit; they do not require network access at runtime.
