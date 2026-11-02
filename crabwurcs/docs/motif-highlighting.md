# Motif highlighting

[`render_svg_with_motifs`] de-emphasises every residue outside a motif of
interest, so a specific substructure stands out against the rest of the glycan.
This is the same GlycoDraw-compatible look used by the `crabwurcs render
--highlight-motif` CLI flag.

[`render_svg_with_motifs`]: https://docs.rs/crabwurcs-snfg/latest/crabwurcs_snfg/fn.render_svg_with_motifs.html

## At a glance

Target: `Neu5Ac(a2-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc`. Highlighted motif:
`Gal(b1-4)[Fuc(a1-3)]GlcNAc`.

| Without highlight | With motif highlighted |
|---|---|
| ![plain](img/examples/motif-before.svg) | ![highlighted](img/examples/motif-after.svg) |

Matched residues and the motif-internal bonds keep their full SNFG colour; the
non-matching `Neu5Ac` branch and its bond drop to the muted palette below.

## How matching works

Matching is structural and wildcard-aware, implemented by
[`find_motif_matches`]. The motif is an ordinary [`ResidueGraph`] parsed from
WURCS, IUPAC (condensed or extended), or GLYCAM.

[`find_motif_matches`]: https://docs.rs/crabwurcs-core/latest/crabwurcs_core/fn.find_motif_matches.html
[`ResidueGraph`]: https://docs.rs/crabwurcs-core/latest/crabwurcs_core/struct.ResidueGraph.html

- **Injective, directed, non-induced.** Each motif node maps to one target
  node; extra branches and extra residue modifications on the target are
  allowed and do not block a match.
- **Unknown anomers and positions are wildcards.** `Fuc(a1-?)GlcNAc` matches
  `Fuc(a1-3)GlcNAc`, `Fuc(a1-4)GlcNAc`, etc.
- **Generic classes match their whole family.** A motif residue `Hex` matches
  `Glc`, `Man`, `Gal`, and every other hexose; `HexNAc` matches every
  N-acetylhexosamine; `Sia` matches `Neu5Ac`, `Neu5Gc`, `Kdn`, `Neu`. This
  follows `ResidueKind::matches_family`, so a single generic motif captures
  every stereochemical variant in its row of the SNFG table.
- **Union over motifs and occurrences.** Every occurrence of every supplied
  motif is highlighted; passing several motifs builds up the union of their
  matches. An empty motif list is identical to a plain render.

Boundary edges — the bonds between a matched region and the rest of the target
— are treated as outside the motif and dimmed with their target-side residue.

## The muted palette

Dimmed residues stay fully opaque (they must occlude the bonds drawn behind
them) but drop to desaturated versions of the SNFG colours:

| SNFG colour | Muted |
|---|---|
| Blue `#0072BC` | `#CDE7EF` |
| Green `#00A651` | `#CDE9DF` |
| Yellow `#FFD400` | `#FFF6DE` |
| Orange `#F47920` | `#FDE7E0` |
| Pink `#F69EA1` | `#FDF0F1` |
| Purple `#A54399` | `#F1E6ED` |
| Light blue `#8FCCE9` | `#EEF8FB` |
| Brown `#A17A4D` | `#F1E9E5` |
| Red `#ED1C24` | `#F7E0E0` |

Bonds and text outside a match draw in light grey `#D9D9D9`.

## Usage

```rust,ignore
use crabwurcs::{render_svg_with_motifs, RenderOptions};

let target = crabwurcs::iupac::parse_iupac_condensed(
    "Neu5Ac(a2-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc",
).unwrap();

// Wildcard linkage position — matches α3 or α4 fucose.
let motif = crabwurcs::iupac::parse_iupac_condensed("Gal(b1-4)[Fuc(a1-?)]GlcNAc").unwrap();

let svg = render_svg_with_motifs(&target, &[motif], &RenderOptions::default()).unwrap();
```

Repeat the call with several motifs to highlight a union of substructures.
PNG output is available through [`render_png_with_motifs`].

[`render_png_with_motifs`]: https://docs.rs/crabwurcs-snfg/latest/crabwurcs_snfg/fn.render_png_with_motifs.html

## Positions

A [`MotifQuery`] says where each motif residue has to sit, in the order the
residues are read in the motif's text, as GlycoDraw's `highlight_termini_list`:

[`MotifQuery`]: https://docs.rs/crabwurcs-core/latest/crabwurcs_core/struct.MotifQuery.html

| [`MotifPosition`] | The target residue |
|---|---|
| `Flexible` | may be anywhere |
| `Terminal` | has nothing on its non-reducing side |
| `Internal` | has a residue on its non-reducing side |
| `Reducing` | is the reducing end of the glycan |

[`MotifPosition`]: https://docs.rs/crabwurcs-core/latest/crabwurcs_core/enum.MotifPosition.html

A candidate parent of an undefined fragment may carry that fragment, so it is
not `Terminal`.

```rust,ignore
use crabwurcs::{HighlightSelection, MotifPosition::*, MotifQuery, RenderOptions};

let target = crabwurcs::iupac::parse_iupac_condensed(
    "Neu5Ac(a2-3)Gal(b1-4)GlcNAc(b1-2)Man(a1-3)[Gal(b1-4)GlcNAc(b1-2)Man(a1-6)]Man",
).unwrap();
let lacnac = crabwurcs::iupac::parse_iupac_condensed("Gal(b1-4)GlcNAc").unwrap();

// Only the LacNAc that is not capped: its Gal must be terminal.
let query = MotifQuery::with_positions(lacnac, &[Terminal, Flexible]).unwrap();
let selection = HighlightSelection::from_matches(&query.find(&target).unwrap());
let svg = crabwurcs::render_svg_with_selection(&target, &selection, &RenderOptions::default());
```

`HighlightSelection::from_residues` selects residues by node index instead, for
residues known from PDB provenance or a fragment ion, with the bonds between
them.

## Named motifs

`crabwurcs::named_motifs` lists a library of 27 motifs, and a motif string in
`auto` format may be one of their names. Names ignore case, spaces and hyphens,
so `sialyl-Lewis X` finds `SialylLewisX`.

| Group | Names |
|---|---|
| Lewis epitopes | `LewisX`, `LewisA`, `LewisY`, `LewisB`, `SialylLewisX`, `SialylLewisA` |
| Blood groups and others | `HType1`, `HType2`, `BloodGroupA`, `BloodGroupB`, `AlphaGal`, `Sda` |
| Antennae | `LacNAc`, `LacNAcType1`, `LacdiNAc`, `PolyLacNAc`, `Sialyl3LacNAc`, `Sialyl6LacNAc` |
| N-glycan cores | `NGlycanCore`, `CoreFucose`, `CoreFucoseAlpha3`, `BisectingGlcNAc` |
| O-glycan cores | `Core1`, `Core2`, `Core3`, `Core4`, `SialylTn` |

An epitope that is a cap by definition requires its non-reducing residues to be
terminal, so `LewisX` is not found inside sialyl-Lewis X. A core requires its
last residue to be the reducing end, so `CoreFucose` is not an antennal fucose.
`LacNAc` is flexible, so it is also found in a capped antenna. A core 2 glycan
holds `Core1` as well as `Core2`, because the first is a substructure of the
second.

`crabwurcs::annotate_motifs` names the library motifs a glycan holds.

## Motif constraints

A motif must describe a connected, directed tree of known residues. The matcher
returns [`MotifError`] otherwise:

[`MotifError`]: https://docs.rs/crabwurcs-core/latest/crabwurcs_core/enum.MotifError.html

| Variant | Meaning |
|---|---|
| `Empty` | The motif graph has no residues. |
| `Composition` | Compositions cannot be used as motifs. |
| `Disconnected` | The motif is not a single connected graph. |
| `NotTree` | The motif is not a directed tree rooted at its reducing end. |
| `CycleOrRepeat` | The motif contains a cycle or a repeat-closing edge. |
| `UndefinedLinkage` | The motif contains an undefined (candidate-parent) linkage. |
| `UndefinedModification` | The motif contains an undefined (candidate-parent) modification. |
| `PositionCount` | The positions given do not match the motif's residues one to one. |
