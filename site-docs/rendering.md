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
