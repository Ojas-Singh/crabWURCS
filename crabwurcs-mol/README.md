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
crabwurcs-mol = "0.3.0"
```

### RDKit Backend
```toml
[dependencies]
crabwurcs-mol = { version = "0.3.0", features = ["rdkit-backend"] }
```

See [README-BACKEND.md](README-BACKEND.md) for RDKit build requirements.

## Usage

```rust
use crabwurcs_core::{parse_wurcs, write_wurcs};
use crabwurcs_mol::{ChemFormat, molecule_from_wurcs, wurcs_from_molecule};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let smiles = "OC[C@H]1O[C@H](O)[C@H](O)[C@@H]1O";
    let graph = wurcs_from_molecule(smiles, ChemFormat::Smiles)?;
    println!("{}", write_wurcs(&graph)?);

    let graph = parse_wurcs("WURCS=2.0/1,1,0/[a2122h-1a_1-5]/1/")?;
    let mol = molecule_from_wurcs(&graph, ChemFormat::Mol)?;
    println!("{mol}");
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
