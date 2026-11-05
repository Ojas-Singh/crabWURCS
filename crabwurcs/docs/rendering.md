# SNFG rendering

[`render_svg`] and friends turn a [`ResidueGraph`] into an SNFG (Symbol
Nomenclature for Glycans) figure. The output follows SNFG 2.0.4: shape encodes
the chemical class, colour encodes the stereochemical family, and the layout
places the reducing end on the right with children extending to the left.

[`render_svg`]: https://docs.rs/crabwurcs-snfg/latest/crabwurcs_snfg/fn.render_svg.html
[`ResidueGraph`]: https://docs.rs/crabwurcs-core/latest/crabwurcs_core/struct.ResidueGraph.html

## Shapes

Each monosaccharide family maps to a distinct geometric primitive. The gallery
below shows one representative per shape; see
[Supported monosaccharides](supported-monosaccharides.md) for the full table.

| Shape | Family | Example |
|---|---|---|
| Circle | Hexoses | ![Glc](img/symbols/Glc.svg) |
| Square | N-Acetylhexosamines | ![GlcNAc](img/symbols/GlcNAc.svg) |
| Notched square | Hexosamines (amino, non-acetylated) | ![GlcN](img/symbols/GlcN.svg) |
| Divided diamond (top) | Hexuronates | ![GlcA](img/symbols/GlcA.svg) |
| Divided diamond (bottom) | Iduronic acid | ![IdoA](img/symbols/IdoA.svg) |
| Triangle | 6-Deoxyhexoses | ![Fuc](img/symbols/Fuc.svg) |
| Divided triangle | 6-Deoxyhexosamines | ![FucNAc](img/symbols/FucNAc.svg) |
| Flat rectangle | 2,6-Dideoxyhexoses | ![Oli](img/symbols/Oli.svg) |
| Star | Pentoses | ![Xyl](img/symbols/Xyl.svg) |
| Diamond | Nonulosonic acids (sialic acids) | ![Neu5Ac](img/symbols/Neu5Ac.svg) |
| Flat diamond | 3,9-Dideoxy-nonulosonic acids | ![Leg](img/symbols/Leg.svg) |
| Flat hexagon | Unknown / bacterial / muramic acid family | ![Kdo](img/symbols/Kdo.svg) |
| Pentagon | Assigned & ketoses | ![Fru](img/symbols/Fru.svg) |

## Colour palette

The official SNFG RGB palette. Generic classes (unspecified stereochemistry)
draw in white while keeping their family shape.

| Name | Hex | Swatch |
|---|---|---|
| White | `#FFFFFF` | ![white](img/symbols/Hex.svg) |
| Blue | `#0072BC` | ![blue](img/symbols/Glc.svg) |
| Green | `#00A651` | ![green](img/symbols/Man.svg) |
| Yellow | `#FFD400` | ![yellow](img/symbols/Gal.svg) |
| Orange | `#F47920` | ![orange](img/symbols/Xyl.svg) |
| Pink | `#F69EA1` | ![pink](img/symbols/Alt.svg) |
| Purple | `#A54399` | ![purple](img/symbols/Neu5Ac.svg) |
| Light blue | `#8FCCE9` | ![light blue](img/symbols/Neu5Gc.svg) |
| Brown | `#A17A4D` | ![brown](img/symbols/IdoA.svg) |
| Red | `#ED1C24` | ![red](img/symbols/Fuc.svg) |

Family → colour assignment follows the standard mapping: Glc blue, Gal yellow,
Man green, Fuc red, sialic acids purple/light-blue, iduronic brown, and so on.
Generic classes (e.g. `Hex`, `HexNAc`) keep the family shape but fill white.

## Layout

The renderer lays the glycan out as a rooted tree:

- **Root on the right.** The reducing-end residue sits at the right edge of the
  canvas; every child extends one step (100 px) to the left.
- **Post-order vertical packing.** A parent is centred vertically between its
  children; leaf residues are stacked at 100 px vertical spacing.
- **Branch order follows the acceptor position.** Higher acceptor positions sit
  above lower ones — for example the N-glycan α1-6 arm is drawn above α1-3, and
  β1-4 above β1-2.
- **Terminal deoxy-sugar lanes.** Terminal fucose and rhamnose residues are
  drawn vertically aligned with their parent. When a parent carries both a
  α1-3 and α1-6 fucose (or α1-2 and α1-4 rhamnose), they are split above and
  below the parent instead of overlapping; a single such branch and the core
  α6-fucose default to one side so no residue is overprinted.
- **Disconnected components** (compositions, undefined antennae) are all laid
  out rather than silently dropped.

Linkage labels (`α3`, `β4`, etc.) are placed beside each bond, rotated to the
bond angle. Repeat/cyclic and undefined bonds draw as dashed grey lines; an
unknown position renders as `?`.

## Composition layout

WURCS compositions (a multiset of residues with no linkage information) render
as a single row of grouped symbols, each with a `×N` count.

![composition](img/examples/composition.svg)

## Render options

[`RenderOptions`] controls the figure:

[`RenderOptions`]: https://docs.rs/crabwurcs-snfg/latest/crabwurcs_snfg/struct.RenderOptions.html

| Field | Default | Effect |
|---|---|---|
| `colour` | `true` | Fill shapes with the SNFG palette; `false` draws outlines only. |
| `show_labels` | `false` | Draw residue abbreviations inside shapes. Off by default — SNFG convention is that shape + colour is the label. |
| `show_linkages` | `true` | Draw anomeric and position labels on each bond. |
| `font_family` | `Arial, Helvetica, sans-serif` | Font stack for all text. |
| `scale` | `1.0` | Uniform scale of the whole figure. |
| `source_notation` | `None` | Records the exact input notation in the SVG metadata. |
| `orientation` | `RootRight` | Where the reducing end is drawn; see [Orientation](#orientation). |
| `compact` | `false` | Tighter spacing and no linkage text, for heatmap axes and tables. |
| `reducing_end` | `None` | Text past the reducing end, such as `Asn`, `Ser/Thr` or `2-AB`. |
| `repeat_brackets` | `false` | Draw repeat units in brackets with their count. |
| `floating_fragments` | `false` | Draw undefined fragments above the structure. |
| `cleavages` | empty | Fragment ions to mark on their bonds. |

Every option is off by default, and a default drawing is the same in every
release that does not say otherwise. Compositions have no tree, so they ignore
`orientation`, `compact`, `reducing_end`, and the bracket options.

`Assigned` pentagons always show their first-letter label regardless of
`show_labels`, because their identity is not recoverable from shape or colour.

## Orientation

[`Orientation`] names where the reducing end sits. The higher acceptor position
of a branch stays above the lower one in the horizontal drawings and to its
right in the vertical ones, so an N-glycan's α1-6 arm is drawn above, or to the
right of, its α1-3 arm.

[`Orientation`]: https://docs.rs/crabwurcs-snfg/latest/crabwurcs_snfg/enum.Orientation.html

| Value | Reducing end | Also accepted |
|---|---|---|
| `RootRight` | right (the SNFG convention) | `right-to-left`, `rl`, `horizontal` |
| `RootLeft` | left | `left-to-right`, `lr` |
| `RootBottom` | bottom | `bottom-up`, `bt`, `vertical` |
| `RootTop` | top | `top-down`, `tb` |

Vertical drawings keep the linkage text of upright bonds horizontal, beside the
bond, and put sulfate labels beside the symbol.

| `RootBottom` with `reducing_end: "Asn"` | `compact` |
|---|---|
| ![vertical](img/examples/vertical.svg) | ![compact](img/examples/compact.svg) |

## Uncertainty

Repeat units, and fragments whose attachment is not known, have two drawings.
The default keeps a dashed bond; `repeat_brackets` and `floating_fragments` use
the SNFG convention.

| Repeat brackets | Floating fragments |
|---|---|
| ![repeat](img/examples/repeat.svg) | ![floating](img/examples/floating.svg) |

A repeat bracket goes around the unit: everything from its reducing-side residue
to its far end, with side branches, but without a cap beyond that far end. The
count reads `n`, `n = 3`, `n = 2–5`, `n ≥ 2` or `n ≤ 5`.

Floating fragments are grouped by their candidate parents. Each group sits over
one bracket that spans its candidates, with a tick above each of them.

## Fragment ions

`cleavages` marks Domon–Costello ions on their bonds, with a stroke across the
bond that hooks towards the fragment that is kept. `crabwurcs::glycan_fragments`
lists the ions and `crabwurcs::cleavage_mark` turns one into a mark.

![cleavage](img/examples/cleavage.svg)

## Grids

[`render_grid_svg`] draws many glycans on one canvas, in rows, with optional
captions. Each cell keeps its own drawing and notations.

[`render_grid_svg`]: https://docs.rs/crabwurcs-snfg/latest/crabwurcs_snfg/fn.render_grid_svg.html

## Output formats

- **SVG** — [`render_svg`] / [`render_svg_with_options`]. Vector output,
  transparent background, embeds accessible `<title>`/`<desc>` and a
  `<metadata>` block carrying the canonical IUPAC condensed and WURCS notation.
- **PNG** — [`render_png`] / [`render_png_with_options`]. Transparent RGBA
  raster at twice the SVG dimensions.

[`render_svg_with_options`]: https://docs.rs/crabwurcs-snfg/latest/crabwurcs_snfg/fn.render_svg_with_options.html
[`render_png`]: https://docs.rs/crabwurcs-snfg/latest/crabwurcs_snfg/fn.render_png.html
[`render_png_with_options`]: https://docs.rs/crabwurcs-snfg/latest/crabwurcs_snfg/fn.render_png_with_options.html

```rust,ignore
use crabwurcs::{render_svg_with_options, render_png, RenderOptions};

let graph = crabwurcs::iupac::parse_iupac_condensed("Gal(b1-4)GlcNAc").unwrap();

let svg = render_svg_with_options(&graph, &RenderOptions::default()).unwrap();
let png = render_png(&graph).unwrap(); // transparent 2× RGBA
```

## Embedded metadata

Every SVG carries an invisible `<metadata>` block (`crabwurcs:notations`) with:

- the canonical IUPAC condensed form,
- the canonical WURCS form, and
- the original source notation and detected format (when available).

Each field is marked `available="true|false"` so a downstream tool knows when a
canonical form could not be produced (for example an `Assigned` residue that
cannot be represented in WURCS) without failing the render.
