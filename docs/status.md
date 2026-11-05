# Development status

Version 0.3.3 is the current synchronized Rust and Python release.

## Release-gated capabilities

- Lossless WURCS 2.0 parsing/writing for supported editable and preserved
  constructs, including ambiguous positions, probabilities, repeats, cycles,
  and undefined fragments.
- IUPAC condensed/extended and GLYCAM parsing/writing through the shared graph.
- SMILES and V3000 MOL/SDF conversion with pure-Rust stereochemistry.
- All 74 concrete registry residues construct and reverse-extract through each
  molecular format without relying on a GlycoShape corpus match.
- PDB/mmCIF extraction uses pinned CCD mappings, GLYCAM decoding, all concrete
  canonical registry names, and coordinate-graph fallback with provenance.
- All 87 SNFG 2.0.4 entries render as SVG and PNG, in four orientations, with
  compact layouts, repeat brackets, floating fragments and cleavage marks.
- Formula, monoisotopic and average mass, composition and glycosidic fragment
  ions. Every corpus formula is checked against its molecule and its SMILES.
- Structure comparison across notations, and a library of 27 named motifs.
- CPython 3.9+ bindings expose conversion, extraction, rendering, typed errors,
  immutable result records, typing metadata, and a CLI.

## Deliberate limitations

- Ensemble-valued WURCS constructs never collapse to arbitrary molecules, or
  to one mass.
- PDB/mmCIF output and conformer generation are not implemented.
- Generic and display-only SNFG entries parse/render but cannot always produce
  molecular or GLYCAM output.
- PyPy and free-threaded CPython are not release targets for 0.3.3.

## Verification

The normal suite uses committed fixtures and bundled corpus tables. The full
external GlycoShape JSON audit remains an opt-in maintenance test documented in
`CONTRIBUTING.md`. CI uses Rust 1.97.1 and tests Cargo packages plus CPython 3.9
wheel installation from a clean checkout.
