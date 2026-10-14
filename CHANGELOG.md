# Changelog

## 0.3.2 - 2026-10-14

Fixes for defects found while auditing the GlycoShape database. Several
written strings change; they are listed under "Output changes".

- **Conversion is computed.** `convert` (Rust), `crabwurcs convert` (both
  command lines) and SMILES reading/writing in `parse_notation` /
  `write_notation` no longer answer inputs that are records of the bundled
  GlycoShape table from that table. The table is used only on request:
  `convert_with_table` / `TableMode` (Rust), `convert(..., table=...)`
  (Python) and `--table fallback|check` (command lines). `fallback` answers
  an input that cannot be converted from its record; `check` fails when the
  input's record holds another structure. `table_lookup` returns a record's
  notation without converting. `detect_format` no longer consults the table.
- **GlyTouCan accessions are not notations.** The GLYCAM and IUPAC readers no
  longer resolve `G50508SG`, `G60371DN` / `G60371D-N` and `G77147GI` from
  built-in lists; they return an error. The first two remain available
  through the table fallback. `G77147GI` is no longer resolved at all: the
  structure given for it (four open linkages) is that of G92358QN.
- **GLYCAM substituent lists are read**: `[9A]` (O-acetyl), `[2A]`, `[6PC]`
  (phosphocholine), `[2Ac,4Ac]` (N-acetyl on bacillosamine), `[2Me,3Me]`,
  `[3S,6S]`. Before, `A` and `PC` were dropped without an error; an unknown
  code is now an error.
- **Reducing-end anomer**: IUPAC condensed `GlcNAc(b1-` and IUPAC extended
  `β-D-GlcpNAc-(1→` are read with their anomer and ring; GLYCAM, IUPAC
  condensed and IUPAC extended are written with it.
- **Residue names carry configuration and ring size** where they differ from
  the default of the name: `D-Fuc`, `D-Rha`, `L-Man`, `D-Ara`, `L-Gul`,
  `L-Tal`, `L-Qui`; `Araf`, `Galf`, `Fruf`, `Ribf`, `Psif`, `D-Araf`. Extra
  acetamido groups are written (`dHexNAc4NAc`), and acyl groups on the amino
  groups of a name (`Bac2Ac4Ac`, `Pse5Ac7Ac`) are read as N-acyl groups.
  New: `classify_residue_configuration`, `mirror_skeleton_code`,
  `ResidueKind::default_configuration`, `ResidueKind::default_is_furanose`.
- **SNFG**: a residue with the non-default configuration of its name is drawn
  with that name's symbol and a configuration letter inside it (SNFG: "less
  common configurations need to be stated ... by adding letters inside the
  symbol"). D-Fuc was drawn as Gal, D-Rha as Man, L-Qui as Glc, and the
  enantiomers of the other 6-deoxy and amino sugars likewise.
- **PDB**: files with free-text records (the banner `REMARK` lines of
  GlycoShape's files) are read; the extractor reads the `ATOM`/`HETATM`
  records of the first model and the `CONECT` and `LINK` records
  (`pdb_records_used`). GLYCAM residue codes that begin with a linkage digit
  (`0aU`, `0FA`) are recognized in `ATOM` records too, as GlycoShape writes
  them; amino-acid-like codes still need `HETATM`. mmCIF reading is
  unchanged.
- The wwPDB components AHR (α-L-arabinofuranose), ARA and ARB (α/β-L-
  arabinopyranose) were read as D-arabinose; they are L-arabinose.
- IUPAC condensed kept the phosphocholine of `GlcNAc6PCho` only when the
  name was parsed from IUPAC; written from WURCS it was dropped.

### Output changes

- IUPAC condensed of a glycan whose WURCS declares the reducing-end anomer:
  `Gal(b1-4)GlcNAc(b1-` (was `Gal(b1-4)GlcNAc`); IUPAC extended
  `β-D-Galp-(1→4)-β-D-GlcpNAc-(1→` (was `…-D-GlcpNAc`); GLYCAM
  `DGalpb1-4DGlcpNAcb1-OH` (was `DGalpb1-4DGlcpNAc`).
- Residue names as above: `D-Fuc` (was `Fuc`), `D-Rha` (was `Man`), `L-Man`,
  `L-Gul`, `L-Tal`, `L-Qui` (were `Res4`), `D-Ara` (was `Res3`), `Araf` (was
  `Ara`), `Fruf` (was `Fru`), `GlcNAc6PCho` (was `GlcNAc`); GLYCAM `DFucp`
  (was `LFucp`), `DRhap` (was `DManp`), `LAraf`, `DFruf`, `DGalf` (were
  `LArap`, `DFrup`, `DGalp`); IUPAC extended `α-L-Araf` (was `α-L-Arap`).
- GLYCAM substituents use GLYCAM codes in a comma-separated list:
  `DNeup5Ac[9A]` (was `DNeup5Ac9Ac`), `DGlcpNAc[6PC]` (was
  `DGlcpNAc6PCho`), `LRhap[2Me,3Me]` (was `LRhap[2Me3Me]`); N-sulfated
  hexosamines carry their configuration (`DGlcpNS[6S]`, was `GlcpNS[6S]`).
- IUPAC extended places the ring letter after the stem: `β-L-Rhap2Me3Me` (was
  `β-L-Rha2Me3Mep`), `β-D-Bacp2Ac4Ac` (was `β-D-Bac2Ac4Acp`).
- WURCS read from GLYCAM names with `[9A]`, `[2A]` or `[6PC]` now contains the
  substituent.
- `crabwurcs convert` of table records returns computed strings: for example
  the WURCS of GlycoShape GS00728 converts to an IUPAC sequence (was
  `G60371D-N`), SMILES output is crabWURCS's canonical SMILES rather than the
  string stored in the table, and IUPAC strings are crabWURCS's own spelling.
- SNFG drawings of the residues above change symbol and gain a letter; their
  terminal 6-deoxy sugars take the deoxyhexose lanes of the layout. Other
  drawings are unchanged, but the canonical IUPAC condensed in the SVG title,
  description and metadata carries the reducing-end anomer when the input
  declares it (`…GalNAc(a1-` for a GLYCAM name ending in `a1-OH`).
- `pdb-to-wurcs` of a multi-model PDB file reports the glycans of the first
  model only.

## 0.3.1 - 2026-08-12

- Normalize WURCS output deterministically across Rust, Python, and CLI APIs.
- Add strict `normalize_wurcs` library functions and `normalize-wurcs` CLI command.
- Preserve unsupported parsed WURCS during ordinary serialization while explicit normalization returns a typed error.
- Add an offline GlycoShape conformance audit for all 1,780 populated IUPAC
  condensed/extended notations and an opt-in live GlyCosmos API comparison.

All notable changes to crabWURCS are recorded here.

## 0.3.0 - 2026-08-07

- Add the `crabwurcs` Python package for CPython 3.9+ with typed high-level
  conversion, molecular, PDB/mmCIF, SNFG, and CLI APIs.
- Guarantee and test all 74 concrete registry residues across WURCS, IUPAC,
  SMILES, MOL/SDF, and PDB/mmCIF component recognition.
- Add PDB residue provenance and exact SNFG highlight selections.
- Replace local-only default corpus tests with committed fixtures.
- Add clean-checkout CI, multiplatform wheels, trusted PyPI publishing,
  ordered crates.io publishing, GitHub releases, and GitHub Pages docs.

## 0.2.2 - 2026-07-27

- Expand molecular and PDB sugar recognition and correct uncommon anomer and
  residue handling.
