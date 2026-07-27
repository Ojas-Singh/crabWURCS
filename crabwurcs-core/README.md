# crabwurcs-core

[![Crates.io](https://img.shields.io/crates/v/crabwurcs-core)](https://crates.io/crates/crabwurcs-core)
[![Documentation](https://docs.rs/crabwurcs-core/badge.svg)](https://docs.rs/crabwurcs-core)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Core WURCS parser/writer and shared ResidueGraph model for the crabWURCS project.

## Overview

This crate provides the foundational data structures and parsing logic for the crabWURCS toolkit:

- **WURCS 2.0 parser**: Lossless parsing with full support for ambiguous linkage positions
- **ResidueGraph model**: Shared internal representation for glycan structures
- **WURCS writer**: Canonical WURCS generation from ResidueGraph
- **Core types**: Residue, Linkage, and other fundamental types

## Features

- Lossless WURCS 2.0 parsing (tested on 839+ records from GlycoShape)
- Support for ambiguous linkage positions
- Efficient graph-based representation
- Comprehensive error handling

## Usage

```rust
use crabwurcs_core::{WurcsParser, ResidueGraph};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let wurcs = "WURCS=2.0/2.1,3,2.1...";
    let graph = WurcsParser::parse(wurcs)?;
    
    // Access residues and linkages
    for residue in graph.residues() {
        println!("{}", residue.name());
    }
    
    Ok(())
}
```

## Documentation

- **[Full API Documentation](https://docs.rs/crabwurcs-core)** - Complete API reference
- **[Main Project](../)** - Overall project documentation

## License

MIT - See [LICENSE-MIT](../LICENSE-MIT) for details.
