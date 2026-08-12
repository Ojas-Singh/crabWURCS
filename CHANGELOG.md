# Changelog

## Unreleased

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
