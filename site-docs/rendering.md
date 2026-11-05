# SNFG rendering

SVG is returned as text and PNG as bytes:

```python
svg = crabwurcs.render_snfg("Gal(b1-4)GlcNAc")
png = crabwurcs.render_snfg(
    "Fuc(a1-3)GlcNAc(b1-4)Fuc(a1-3)GlcNAc",
    image_format="png",
    highlight_motifs=["Fuc(a1-?)GlcNAc"],
)
```

The renderer implements the SNFG 2.0.4 palette and symbols, accessible SVG
metadata, transparent PNG output, compositions, and wildcard-aware structural
motif highlighting. Motifs accept WURCS, IUPAC, or GLYCAM; SMILES is not a
motif-query notation.

## Layout

Every option is off by default, and a default drawing does not change.

```python
glycan = crabwurcs.Glycan.parse("Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc")

glycan.render(orientation="root-bottom")  # also root-left, root-top
glycan.render(compact=True)               # for heatmap axes and tables
glycan.render(reducing_end="Asn")         # text past the reducing end
```

`orientation` names where the reducing end sits, and also accepts `left-to-right`,
`top-down`, `bottom-up` and `vertical`. The higher acceptor position of a branch
stays above the lower one in horizontal drawings and to its right in vertical
ones, so an N-glycan's α1-6 arm is above, or to the right of, its α1-3 arm.

## Uncertainty

`repeat_brackets=True` draws a repeat unit in square brackets with its count
(`n = 2–5`). `floating_fragments=True` draws undefined antennae and substituents
above the structure, over a bracket across their candidate parents. Without them
both keep a dashed bond.

## Highlighting

```python
target = crabwurcs.Glycan.parse(
    "Neu5Ac(a2-3)Gal(b1-4)GlcNAc(b1-2)Man(a1-3)[Gal(b1-4)GlcNAc(b1-2)Man(a1-6)]Man"
)

# By name, from the motif library.
target.render(highlight_motifs=["LacNAc"])

# Only an uncapped LacNAc: one position per residue, in reading order.
target.render(
    highlight_motifs=["Gal(b1-4)GlcNAc"],
    highlight_positions=[("terminal", "flexible")],
)

# Exact residues, for example from PDB provenance or a fragment ion.
target.render(highlight_residues=[0, 1])

# Fragment ions, by label (see Glycan.fragments()).
target.render(cleavages=["B2\u03b1", "Y3\u03b2"])
```

Positions are `terminal`, `internal`, `reducing` and `flexible`.
`Glycan.residues()` lists the node indices.

## Figures

`crabwurcs.figures` works with many glycans at once.

```python
from crabwurcs import figures

svg = figures.render_grid(glycans, columns=4, labels=names, compact=True)
png = figures.render_grid(glycans, image_format="png")
figures.write_images({"G0F": g0f, "G1F": g1f}, "images/", image_format="png")
```

`annotate_figure` replaces the glycan labels of a matplotlib figure with their
drawings, for heatmaps, volcano plots and bar charts. `annotate_svg` does the
same for any SVG with text labels.

```python
figure, axes = pyplot.subplots()
axes.imshow(abundance)
axes.set_yticks(range(len(glycans)), glycans)
figure.subplots_adjust(left=0.4)

figures.annotate_figure(figure, "heatmap.png", height=70)
figures.annotate_figure(figure, "heatmap.svg", labels={"G0F": g0f, "G1F": g1f})
```

A label counts when every residue in it is recognized, so axis titles and tick
numbers stay as they are; `labels` maps short names to glycans. All drawings
share one scale, so a symbol is the same size in each: `height` is that of the
tallest. Leave room around the axes, because matplotlib lays the figure out for
the text the drawings replace.

matplotlib is optional, and only `annotate_figure` needs it.
