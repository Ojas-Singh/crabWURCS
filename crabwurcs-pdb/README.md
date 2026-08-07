# crabwurcs-pdb

[![Crates.io](https://img.shields.io/crates/v/crabwurcs-pdb)](https://crates.io/crates/crabwurcs-pdb)
[![Documentation](https://docs.rs/crabwurcs-pdb/badge.svg)](https://docs.rs/crabwurcs-pdb)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

PDB/mmCIF glycan extraction for the crabWURCS project.

## Overview

This crate provides glycan structure extraction from PDB and mmCIF coordinate files:

- **PDB/mmCIF parsing**: Read standard protein data bank files
- **Glycan reconstruction**: Build glycosidic graphs from CONECT/covalent records
- **Component recognition**: Bundled CCD and GLYCAM residue ID resolution
- **Coordinate fallback**: Atom/coordinate-graph recognition for renamed components

## Features

- **1,863+ exact matches**: Tested against GlycoShape PDB/GLYCAM files
- **Complex structure support**: Branches, furanoses, uronic acids, amino sugars, sialic acids
- **Modification handling**: O/N-sulfation, methylation, acetylation, phosphocholine
- **Bundled data**: Generated CCD component table and GLYCAM decoder
- **Private component support**: Fallback to coordinate/atom-graph recognition

## Usage

```rust
use crabwurcs_core::write_wurcs;
use crabwurcs_pdb::extract_glycans_with_provenance_from_file;
use std::path::Path;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let glycans = extract_glycans_with_provenance_from_file(Path::new("glycan.pdb"))?;
    for glycan in glycans {
        println!("Extracted: {}", write_wurcs(&glycan.graph)?);
        println!("Source residues: {:?}", glycan.residues);
    }
    Ok(())
}
```

## Supported Structures

- **Core structures**: Standard glycan scaffolding
- **Branching**: Complex antennary structures
- **Furanoses**: Five-membered ring sugars
- **Modified residues**: Uronic acids, amino sugars, sialic acids
- **Chemical modifications**: Sulfation, methylation, acetylation, phosphocholine

## Testing Coverage

- 1,863 exact semantic matches out of 1,886 GlycoShape PDB/GLYCAM files
- Comprehensive component recognition
- Fallback coordinate-based detection

## Documentation

- **[Full API Documentation](https://docs.rs/crabwurcs-pdb)** - Complete API reference
- **[Main Project](../)** - Overall project documentation

## License

MIT - See [LICENSE-MIT](../LICENSE-MIT) for details.
