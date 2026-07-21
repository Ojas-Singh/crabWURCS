# crabWURCS

A pure-Rust toolkit for glycan notation conversion, chemical-structure interop, and SNFG rendering.

## Features

- **Multi-format conversion**: Convert between WURCS, IUPAC (condensed/extended), GLYCAM, SMILES, MOL, and SDF formats
- **Structure extraction**: Extract glycan structures from PDB and GLYCAM coordinate files
- **SNFG rendering**: Generate publication-quality SNFG (Symbol Nomenclature for Glycans) SVG graphics
- **Lossless parsing**: Lossless WURCS 2.0 parsing with full support for ambiguous linkage positions
- **Pure Rust**: No external C/C++ dependencies by default, with optional RDKit backend

## Installation

Add this to your `Cargo.toml`:

```toml
[dependencies]
crabwurcs = "0.1.0"
```

For CLI installation:

```bash
cargo install crabwurcs-cli
```

## Quick Start

### Library Usage

```rust
use crabwurcs::prelude::*;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Convert between glycan formats
    let iupac = "Gal(b1-4)GlcNAc";
    let wurcs = crabwurcs::convert(iupac, Format::Wurcs)?;
    
    // Render to SNFG SVG
    let svg = crabwurcs::render_snfg(iupac)?;
    
    Ok(())
}
```

### CLI Usage

```bash
# Convert formats (auto-detects input)
crabwurcs convert --to wurcs 'Gal(b1-4)GlcNAc'
crabwurcs convert --to glycam 'β-D-Galp-(1→4)-D-GlcNAc'

# Render SNFG graphics
crabwurcs render --output glycan.svg 'Gal(b1-4)GlcNAc'

# Extract from PDB files
crabwurcs pdb-to-wurcs --to iupac-condensed glycan.pdb

# Handle compositions and complex structures
crabwurcs convert --to wurcs '{GlcNAc}2,{Man}3,{Fuc}1'
```

## Architecture

The project is organized as a workspace of specialized crates:

```
crabwurcs-core    — Core WURCS parser/writer + shared ResidueGraph model
crabwurcs-iupac   — IUPAC condensed/extended and GLYCAM converters
crabwurcs-mol     — MOL/SDF/SMILES molecular structure handling
crabwurcs-pdb     — PDB/mmCIF glycan extraction
crabwurcs-snfg    — SNFG SVG rendering
crabwurcs         — Unified facade crate
craburcs-cli      — Command-line interface
```

All formats convert through a shared `ResidueGraph` model, ensuring consistent representation across different notations.

## Documentation

- **[Status & Progress](docs/status.md)**: Detailed development status and testing coverage
- **[API Documentation](https://docs.rs/crabwurcs)**: Full API reference
- Repository: [https://github.com/Ojas-Singh/crabWURCS](https://github.com/Ojas-Singh/crabWURCS)

## License

This project is licensed under the [MIT License](LICENSE-MIT).

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.
