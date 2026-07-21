# Development Status & Progress

This document contains detailed development status, testing coverage, and remaining work items.

## Current Implementation Status

The notation and drawing pipeline is working end to end with comprehensive test coverage:

### Format Conversion
- **WURCS 2.0 parsing**: Lossless for all 839 WURCS records in `GLYCOSHAPE.json`, including ambiguous linkage positions
- **IUPAC formats**: Both IUPAC condensed and IUPAC extended parse into a shared residue graph
- **GLYCAM support**: Full GLYCAM parsing and generation
- **Cross-format consistency**: Of the 839 records containing both IUPAC and WURCS, 614 have exact rooted residue/linkage signatures and 225 differ only because the WURCS record omits a reducing-end ring closure

### Testing Coverage
- All 839 extended-IUPAC records round-trip losslessly
- All 943 GLYCAM records round-trip losslessly  
- All 938 isomeric-SMILES pairs from the GlycoShape corpus are supported
- Every included notation is tested against every other output format

### Molecular Structure Support
- **SMILES/MOL/SDF**: Pure-Rust molecular backend parses all 938 corpus molecules and canonicalizes them stably
- **De novo construction**: Finite defined WURCS graphs not in the corpus are constructed as stereochemical molecular graphs
- **Aglycone extraction**: Glycan structures can be extracted from single-bond aglycone attachments (e.g., methyl glycosides)

### Advanced Features
- **MolWURCS extraction**: Previously unseen molecular graphs enter a de-novo pure-Rust MolWURCS extractor with full ring detection, stereochemistry, and MAP substituent support
- **WURCS compositions**: Whole-structure repeats, cyclic closures, and undefined fragment parent candidates represented explicitly
- **Cross-linking**: Common cross-linked WURCS MAP bridges (phosphate, phosphoethanolamine) survive graph edits and convert through all formats
- **Undefined substituents**: Retained as candidate-parent modifications rather than being discarded

### SNFG Rendering
- Collision-free tidy-tree layout with SNFG shapes/colors
- Bond-aligned linkage labels with transparent backgrounds
- Accessible SVG output
- Special handling for terminal fucose branches and fructofuranose
- All 938 molecular-corpus records render without unknown-symbol fallback

### PDB/GLYCAM Extraction  
- Reconstructs real glycosidic graphs from `CONECT`/covalent records
- Handles branches, furanoses, uronic acids, amino sugars, sialic acids
- Supports O/N-sulfation, methylation, acetylation, and phosphocholine
- Audit results: 1,863 exact semantic matches out of 1,886 GlycoShape PDB/GLYCAM files

## Remaining Major Work

1. **Extend MolWURCS coverage**: Beyond the 938 verified GlycoShape molecules with additional MAP templates
2. **IUPAC completeness**: Add representations for undefined substituent fragments and complex nested/partial repeat semantics  
3. **SNFG expansion**: Extend residue dictionary and canonical branch ordering against complete SNFG fixtures
4. **Structure audit**: Expand beyond the current GlycoShape sample with explicit mmCIF fixtures for uncommon residues

## Technical Notes

### Scope Limitations
- Broader MAP chemistry beyond represented substituent families remains MolWURCS porting work
- Polymeric or ambiguous WURCS constructs that don't denote one finite molecule are future work
- Sugar discrimination beyond the 938-structure audit is ongoing
- Non-glycans return `NoGlycanFound`

### Toolchain
- **Minimum Rust version**: 1.88 (required by the pure-Rust chemistry backend's use of let-chain syntax)
- **Development version**: Currently tested on Rust 1.97
- **Dependencies**: Pure-Rust by default with optional RDKit backend via `rdkit` feature

## Design Decisions

- **Pure Rust default**: No C/C++ dependencies by default; optional RDKit backend available
- **No GlycoCT**: Only WURCS and IUPAC formats are supported
- **License**: MIT license for maximum compatibility
- **Testing**: Comprehensive corpus-based testing against GlycoShape reference data
