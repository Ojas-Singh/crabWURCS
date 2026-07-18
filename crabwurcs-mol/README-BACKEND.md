# RDKit backend build notes

RDKit is optional and is not required for the bundled GlycoShape molecular
corpus. The default pure-Rust `chematic` backend recognizes equivalent SMILES
graphs, reads/writes corpus-backed MOL/SDF, and performs conservative de-novo
aldose, ketose, nonulosonic, pyranose, and furanose extraction. A lookup-free
audit is semantically exact in both directions for all 938 bundled molecules,
including WURCS → canonical SMILES → WURCS. This document concerns optional
augmentation beyond that audited chemistry.

`crabwurcs-mol` uses the [`rdkit`](https://crates.io/crates/rdkit) crate
(NOT `rdkit-rs`, which is a different, much less mature 0.1.0 crate also on
crates.io — easy to grab the wrong one) as its chemistry backend, gated
behind the `rdkit-backend` feature so the rest of the workspace builds
without it.

## Why RDKit and not OpenBabel

- RDKit is BSD-3-Clause. OpenBabel is GPL-2.0, which would pull a
  statically-linked crabWURCS binary under GPL obligations. If crabWURCS
  is meant to carry a permissive license, RDKit is the only real option
  between the two.
- RDKit's ring perception / aromaticity / valence handling is more
  complete, and MolWURCS's own sugar-vs-non-sugar discrimination algorithm
  (the "modification count" idea from its paper) is built on exactly that
  kind of primitive.

## Build requirements

- Requires RDKit **2023.09.1 or newer**, built with **static libraries**.
  A normal `apt install librdkit-dev` (dynamic libs) is not sufficient on
  Linux per the `rdkit` crate's own docs.
- The crate's upstream CI publishes prebuilt static tarballs for both
  amd64 and arm64 (built on Ubuntu 22.04) — worth checking their releases
  page for a current URL before building RDKit from source yourself:
  <https://github.com/rdkit-rs/rdkit>
- Alternative: a conda-forge RDKit install (Mac Homebrew and conda-forge
  builds are explicitly supported) — the crate's `dynamic-linking-from-conda`
  feature on `rdkit-sys` exists for this path if static linking proves
  impractical in your environment.

## ARM64 / Ampere note

This is a genuinely well-trodden path, not a cross-compilation edge case:
Ampere-class hardware is standard aarch64 Linux server hardware, and both
RDKit upstream and the `rdkit` crate's own CI already build and test
against arm64. Expect this to be no harder than the amd64 build — plan CI
to build/fetch RDKit for both `x86_64-unknown-linux-gnu` and
`aarch64-unknown-linux-gnu` targets from day one rather than treating ARM
as an afterthought.

## TODO for whoever picks this up

- [ ] Confirm current `rdkit` crate version and API surface against
      crates.io (this scaffold pinned `rdkit = "0.4"` as of the RDKit
      crate ecosystem check done while scaffolding — re-verify, this
      ecosystem is young and moves fast).
- [ ] Decide: vendor a static RDKit build per-arch in CI (matching
      upstream's own approach), or depend on conda-forge at build time.
- [ ] Wire up a `xtask` or CI step that builds with `--features
      rdkit-backend` on both x86_64 and aarch64 runners so regressions in
      ARM support are caught immediately, not discovered later.
