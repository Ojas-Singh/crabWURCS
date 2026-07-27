# crabwurcs-mol

[![Crates.io](https://img.shields.io/crates/v/crabwurcs-mol)](https://crates.io/crates/crabwurcs-mol)
[![Documentation](https://docs.rs/crabwurcs-mol/badge.svg)](https://docs.rs/crabwurcs-mol)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

MOL/SDF/SMILES molecular structure handling for the crabWURCS project.

## Overview

This crate provides molecular structure handling capabilities with both pure-Rust and optional RDKit backends:

- **Pure-Rust backend**: Complete SMILES/MOL/SDF support using hematic
- **Optional RDKit backend**: Industrial-strength chemistry via RDKit
- **MolWURCS extraction**: De novo structure extraction from molecular graphs
- **Stereochemistry**: Full atom and bond stereochemistry support

## Features

- **Pure-Rust default**: No C/C++ dependencies by default
- **938+ molecule support**: All GlycoShape corpus molecules tested
- **V3000 MOL/SDF**: Standard format with CFG stereochemistry attributes
- **De novo construction**: Build WURCS graphs from molecular descriptors
- **Cross-linking support**: Handles phosphate, phosphoethanolamine bridges
- **Ring detection**: Full cyclic structure support

## Backend Options

### Pure-Rust (Default)
```toml
[dependencies]
crabwurcs-mol = "0.2.1"
```

### RDKit Backend
```toml
[dependencies]
crabwurcs-mol = { version = "0.2.1", features = ["rdkit-backend"] }
```

See [README-BACKEND.md](README-BACKEND.md) for RDKit build requirements.

## Usage

```rust
use crabwurcs_mol::{MolParser, MolWurcsExtractor};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Parse SMILES
    let smiles = "OC[C@H]1O[C@H](O)[C@H](O)[C@@H]1O";
    let mol = MolParser::parse_smiles(smiles)?;
    
    // Extract WURCS from molecular graph
    let wurcs = MolWurcsExtractor::extract(&mol)?;
    
    // Read MOL file
    let mol = MolParser::parse_mol_file("structure.mol")?;
    
    Ok(())
}
```

## Testing Coverage

- All 938 isomeric-SMILES pairs from GlycoShape corpus supported
- Stable canonicalization
- Full stereochemistry preservation

## Documentation

- **[Full API Documentation](https://docs.rs/crabwurcs-mol)** - Complete API reference
- **[Backend Documentation](README-BACKEND.md)** - RDKit setup guide
- **[Main Project](../)** - Overall project documentation

## License

MIT - See [LICENSE-MIT](../LICENSE-MIT) for details.
