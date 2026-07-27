# crabwurcs-iupac

[![Crates.io](https://img.shields.io/crates/v/crabwurcs-iupac)](https://crates.io/crates/crabwurcs-iupac)
[![Documentation](https://docs.rs/crabwurcs-iupac/badge.svg)](https://docs.rs/crabwurcs-iupac)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

IUPAC condensed/extended and GLYCAM converters for the crabWURCS project.

## Overview

This crate provides conversion between IUPAC, GLYCAM formats and the shared ResidueGraph model:

- **IUPAC Condensed**: Parse and generate condensed IUPAC notation
- **IUPAC Extended**: Full extended IUPAC notation support
- **GLYCAM**: Complete GLYCAM format parsing and generation
- **Cross-conversion**: Bidirectional conversion between all formats

## Features

- Lossless round-trip conversion (839+ IUPAC records tested)
- Full GLYCAM support (943+ records tested)
- Cross-format consistency validation
- Comprehensive error messages

## Usage

```rust
use crabwurcs_iupac::{IupacParser, GlycamParser};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Parse IUPAC condensed
    let iupac = "Gal(b1-4)GlcNAc";
    let graph = IupacParser::parse_condensed(iupac)?;
    
    // Parse GLYCAM
    let glycam = "β-D-Galp-(1→4)-D-GlcNAc";
    let graph = GlycamParser::parse(glycam)?;
    
    // Generate back to IUPAC
    let output = graph.to_iupac_condensed()?;
    
    Ok(())
}
```

## Testing Coverage

- All 839 extended-IUPAC records round-trip losslessly
- All 943 GLYCAM records round-trip losslessly
- Every notation tested against every other output format

## Documentation

- **[Full API Documentation](https://docs.rs/crabwurcs-iupac)** - Complete API reference
- **[Main Project](../)** - Overall project documentation

## License

MIT - See [LICENSE-MIT](../LICENSE-MIT) for details.
