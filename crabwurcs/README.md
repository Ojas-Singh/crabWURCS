# crabwurcs

[![Crates.io](https://img.shields.io/crates/v/crabwurcs)](https://crates.io/crates/crabwurcs)
[![Documentation](https://docs.rs/crabwurcs/badge.svg)](https://docs.rs/crabwurcs)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A pure-Rust toolkit for glycan notation conversion, chemical-structure interop, and SNFG rendering.

## Overview

This is the **facade crate** for the crabWURCS project, providing a unified API that combines functionality from all specialized sub-crates:

- `crabwurcs-core` - Core WURCS parser/writer + shared ResidueGraph model
- `crabwurcs-iupac` - IUPAC condensed/extended and GLYCAM converters
- `crabwurcs-mol` - MOL/SDF/SMILES molecular structure handling
- `crabwurcs-pdb` - PDB/mmCIF glycan extraction
- `crabwurcs-snfg` - SNFG SVG rendering

## Features

- **Multi-format conversion**: Convert between WURCS, IUPAC (condensed/extended), GLYCAM, SMILES, MOL, and SDF formats
- **Structure extraction**: Extract glycan structures from PDB and GLYCAM coordinate files
- **Molecular chemistry**: Construct concrete WURCS backbones/MAP graphs and round-trip stereochemical SMILES or V3000 MOL/SDF
- **SNFG rendering**: Generate publication-quality SNFG (Symbol Nomenclature for Glycans) SVG or PNG graphics
- **Pure Rust**: No external C/C++ dependencies by default, with optional RDKit backend

## Installation

Add this to your `Cargo.toml`:

```toml
[dependencies]
crabwurcs = "0.3.0"
```

For optional RDKit backend support:

```toml
[dependencies]
crabwurcs = { version = "0.3.0", features = ["rdkit-backend"] }
```

## Quick Start

```rust
use crabwurcs::prelude::*;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Convert between glycan formats
    let iupac = "Gal(b1-4)GlcNAc";
    let wurcs = crabwurcs::convert(iupac, Format::Wurcs)?;
    
    // Render to SNFG SVG
    let svg = crabwurcs::render_snfg(iupac)?;
    
    // Extract from PDB files
    let glycans = crabwurcs::extract_from_pdb("glycan.pdb")?;
    
    Ok(())
}
```

## API Documentation

- **[Full API Documentation](https://docs.rs/crabwurcs)** - Complete API reference
- **[Main Project README](../README.md)** - Project overview and CLI usage
- **[Development Status](../docs/status.md)** - Implementation details and testing coverage

## Crate Structure

This crate re-exports all functionality from the specialized sub-crates:

```rust
use crabwurcs::prelude::*;

// Core types and models
use crabwurcs::{ResidueGraph, Format, RenderOptions};

// Format conversion
use crabwurcs::{convert, convert_to_wurcs, convert_to_iupac};

// SNFG rendering
use crabwurcs::{render_snfg, render_svg, render_png};

// PDB extraction
use crabwurcs::{extract_from_pdb, extract_from_glycam};
```

## License

This project is licensed under the [MIT License](../LICENSE-MIT).

## Repository

[https://github.com/Ojas-Singh/crabWURCS](https://github.com/Ojas-Singh/crabWURCS)
