# crabwurcs-cli

[![Crates.io](https://img.shields.io/crates/v/crabwurcs-cli)](https://crates.io/crates/crabwurcs-cli)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Command-line interface for the crabWURCS glycan toolkit.

## Overview

This is the CLI tool for the crabWURCS project, providing command-line access to all glycan conversion, rendering, and extraction capabilities.

## Installation

```bash
cargo install crabwurcs-cli
```

## Features

- **Format conversion**: Auto-detect input format and convert to WURCS, IUPAC, GLYCAM, SMILES, MOL, or SDF
- **SNFG rendering**: Generate publication-quality SVG or PNG graphics
- **PDB extraction**: Extract glycan structures from PDB and GLYCAM files
- **Motif highlighting**: Highlight structural motifs with wildcard support
- **Batch processing**: Handle multiple files and compositions

## Usage

### Format Conversion

```bash
# Convert to WURCS (auto-detects input format)
crabwurcs convert --to wurcs 'Gal(b1-4)GlcNAc'

# Convert to GLYCAM
crabwurcs convert --to glycam 'β-D-Galp-(1→4)-D-GlcNAc'

# Convert to IUPAC extended
crabwurcs convert --to iupac-extended 'WURCS=2.0/...'
```

### SNFG Rendering

```bash
# Render to SVG
crabwurcs render --output glycan.svg 'Gal(b1-4)GlcNAc'

# Render to PNG
crabwurcs render --output glycan.png 'Gal(b1-4)GlcNAc'

# Render with motif highlighting
crabwurcs render \
  --highlight-motif 'Fuc(a1-?)[Gal(b1-?)]GlcNAc' \
  --motif-from iupac-condensed \
  'Neu5Ac(a2-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc'
```

### PDB Extraction

```bash
# Extract glycans from PDB file
crabwurcs pdb-to-wurcs --to iupac-condensed glycan.pdb

# Extract from GLYCAM files
crabwurcs pdb-to-wurcs --to glycam structure.glc
```

### Advanced Features

```bash
# Handle compositions
crabwurcs convert --to wurcs '{GlcNAc}2,{Man}3,{Fuc}1'

# Generic SNFG classes
crabwurcs render 'HexNAc(?1-?)Hex'
crabwurcs render '{Hex}3,{HexNAc}2,{dHex}1'

# Multiple motif highlighting
crabwurcs render \
  --highlight-motif 'Fuc(a1-3)GlcNAc' \
  --highlight-motif 'Gal(b1-4)GlcNAc' \
  'Neu5Ac(a2-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-2)Man'
```

## Command Reference

### `convert`
Convert between glycan notation formats.

```bash
crabwurcs convert --to <FORMAT> <INPUT>
```

**Formats**: `wurcs`, `iupac-condensed`, `iupac-extended`, `glycam`, `smiles`, `mol`, `sdf`

### `render`
Generate SNFG graphics.

```bash
crabwurcs render --output <FILE> <INPUT>
```

**Output formats**: SVG (`.svg`) or PNG (`.png`)

**Options**:
- `--highlight-motif <MOTIF>`: Highlight structural motifs
- `--motif-from <FORMAT>`: Specify motif format
- `--show-labels`: Show residue abbreviations
- `--scale <FACTOR>`: Scale the output

### `pdb-to-wurcs`
Extract glycans from PDB/GLYCAM files.

```bash
crabwurcs pdb-to-wurcs --to <FORMAT> <INPUT_FILE>
```

## Examples

```bash
# Simple conversion
crabwurcs convert --to wurcs 'Gal(b1-4)GlcNAc'

# Complex structure with highlighting
crabwurcs render \
  --highlight-motif 'Fuc(a1-3)GlcNAc' \
  --output sialyl_lewis_x.svg \
  'Neu5Ac(a2-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-3)Gal(b1-4)Glc'

# Batch processing
crabwurcs pdb-to-wurcs --to iupac-condensed --output glycans.txt *.pdb

# Composition rendering
crabwurcs render --output composition.svg '{Hex}3,{HexNAc}2,{dHex}1'
```

## Documentation

- **[Main Project README](../README.md)** - Full project documentation
- **[Development Status](../docs/status.md)** - Implementation details
- **[Rendering Guide](../crabwurcs/docs/rendering.md)** - SNFG rendering details
- **[Motif Highlighting](../crabwurcs/docs/motif-highlighting.md)** - Motif matching

## Repository

[https://github.com/Ojas-Singh/crabWURCS](https://github.com/Ojas-Singh/crabWURCS)

## License

MIT - See [LICENSE-MIT](../LICENSE-MIT) for details.
