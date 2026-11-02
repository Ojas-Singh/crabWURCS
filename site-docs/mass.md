# Mass and fragments

```python
import crabwurcs

man5 = crabwurcs.Glycan.parse(
    "Man(a1-3)[Man(a1-3)[Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc"
)

man5.formula()                               # 'C46H78N2O36'
man5.mass()                                  # 1234.4334
man5.mass(average=True)                      # 1235.1080
man5.mass(adduct="Na")                       # 1257.4226, [M+Na]+
man5.mass(adduct="H", charge=2)              # 618.2240, [M+2H]2+
man5.mass(permethylated=True, adduct="Na")   # 1579.7826
man5.composition()                           # {'HexNAc': 2, 'Hex': 5}
man5.composition_string()                    # 'HexNAc2Hex5'
```

```bash
crabwurcs mass --adduct Na --reducing-end 2-ab 'Gal(b1-4)GlcNAc'
```

## Derivatives

`reducing_end` says what the reducing end carries.

| Value | Adds | Note |
|---|---|---|
| `free` | | the free hemiacetal |
| `reduced` | H2 | the alditol, and one more methylation site |
| `2-ab` | C7H8N2 | 2-aminobenzamide |
| `2-aa` | C7H7NO | 2-aminobenzoic acid |
| `procainamide` | C13H21N3 | |

The label increments are those of reductive amination. A derivative needs a free
reducing end: a cyclic glycan, a glycan linked through its reducing anomeric
carbon such as trehalose, and an alditol raise `MassError`.

`permethylated=True` methylates every free OH, NH and carboxyl group. It works
with `free` and `reduced`. Labelled glycans, and glycans with sulfate or
phosphate, whose acidic groups do not methylate uniformly, raise `MassError`.

## Adducts

`H`, `Na`, `K`, `NH4`, `Li`, `-H`, `Cl`, `HCOO` and `CH3COO`, also written
`[M+Na]+` or `M-H`. `charge` is the number of adducts, each carrying one charge.
Ion masses account for the electron.

## What has a mass

Linkage positions and stereochemistry do not change a mass, so a glycan does not
need to be one concrete molecule to have one. Compositions, `?` linkages,
generic residues such as `Hex`, and undefined fragments all have a mass, and an
exact repeat is multiplied out.

Some things leave the mass open, and raise `MassError`:

- a repeat count that is not a single number,
- a probability on a substituent or a linkage: that is an ensemble of species,
- a residue the parser kept by name without knowing its chemistry,
- residues that are not joined into one glycan.

Residue formulas come from the molecular constructor, and every GlycoShape record
is checked against its constructed molecule and its SMILES. Phosphocholine is
counted as its neutral inner salt, which gives the usual +165.0555 increment.

## Composition

`composition()` counts residue classes in the order glycomics writes them:
`HexNAc`, `Hex`, `Fuc`, `dHex`, `NeuAc`, `NeuGc`, `Kdn`, `HexA`, `HexN`, `Pen`, then
other residues alphabetically, then the substituents a residue carries beyond its
name (`S`, `P`, `Ac`, `Me`, and `Sub` for anything else). The sulfates of GlcNS6S
count; the acetyl of GlcNAc does not.

## Fragment ions

`Glycan.fragments()` lists the Domon–Costello B, C, Y and Z ions of every
glycosidic cleavage.

```python
glycan = crabwurcs.Glycan.parse("Gal(b1-4)GlcNAc(b1-2)Man")

for ion in glycan.fragments(adduct="Na"):
    print(ion.label, round(ion.mz, 4), ion.formula, ion.residues)

svg = glycan.render(cleavages=["B2", "Y2"])
```

```bash
crabwurcs fragments --adduct Na 'Gal(b1-4)GlcNAc(b1-2)Man'
crabwurcs render --cleavage B2 --cleavage Y2 --output marked.svg 'Gal(b1-4)GlcNAc(b1-2)Man'
```

B and C ions keep the non-reducing end, and are numbered from it. Y and Z ions
keep the reducing end, and are numbered from that. A residue past a branch point
adds one Greek letter per branch point on the way, α for the heaviest branch, so
`B1β` is the terminal residue of the lighter arm. A tie goes to the lower node
index.

`residues` holds the node indices an ion keeps, which `highlight_residues`
accepts. `bond` is the edge index of the cleaved bond.

Fragments need one defined tree. Compositions, repeats, cycles and undefined
fragments raise `MassError`. Cross-ring (A and X) ions are not computed.
