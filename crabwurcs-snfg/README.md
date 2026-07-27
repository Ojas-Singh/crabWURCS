# crabwurcs-snfg

[![Crates.io](https://img.shields.io/crates/v/crabwurcs-snfg)](https://crates.io/crates/crabwurcs-snfg)
[![Documentation](https://docs.rs/crabwurcs-snfg/badge.svg)](https://docs.rs/crabwurcs-snfg)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

SNFG (Symbol Nomenclature for Glycans) SVG rendering for the crabWURCS project.

## Overview

This crate implements the complete SNFG 2.0.4 specification for rendering glycan structures:

- **Full SNFG 2.0.4 support**: All shapes, colors, and symbols
- **Vector graphics**: Publication-quality SVG output
- **Raster output**: High-resolution PNG generation
- **Accessible output**: Includes title/description elements and structured metadata
- **Motif highlighting**: Wildcard-aware structural motif matching

## Features

- Collision-free tidy-tree layout
- Complete SNFG 2.0.4 symbol table and RGB palette
- Special handling for terminal fucose branches and fructofuranose
- Transparent linkage labels with bond-aligned rotation
- Support for compositions and disconnected components
- Motif highlighting with GlycoDraw-compatible muted colors

## Usage

```rust
use crabwurcs_snfg::{render_svg, render_png, RenderOptions};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Render to SVG
    let svg = render_svg(&graph)?;
    
    // Render to PNG (transparent RGBA)
    let png = render_png(&graph)?;
    
    // Render with options
    let options = RenderOptions {
        show_labels: true,
        colour: true,
        ..Default::default()
    };
    let svg = render_svg_with_options(&graph, &options)?;
    
    Ok(())
}
```

## Output Formats

- **SVG**: Vector output with embedded metadata and accessibility features
- **PNG**: High-resolution raster output (2x SVG dimensions)

## Embedded Metadata

Every SVG carries invisible `<metadata>` containing:
- Canonical IUPAC condensed form
- Canonical WURCS form  
- Original source notation and detected format (when available)

## Documentation

- **[Full API Documentation](https://docs.rs/crabwurcs-snfg)** - Complete API reference
- **[Rendering Guide](../crabwurcs/docs/rendering.md)** - Detailed rendering documentation
- **[Main Project](../)** - Overall project documentation

## License

MIT - See [LICENSE-MIT](../LICENSE-MIT) for details.
