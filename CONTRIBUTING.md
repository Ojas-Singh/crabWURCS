# Contributing

crabWURCS requires Rust 1.97 or newer. Python binding work also needs Python
3.9+ and maturin 1.14+.

```bash
cargo fmt --all -- --check
cargo clippy --workspace --all-targets -- -D warnings
cargo test --workspace --all-targets
cargo test --workspace --doc
python -m pip install maturin pytest
maturin develop
python -m pytest
./scripts/check_release.sh
```

Do not add runtime network lookups. Generated chemistry data must include its
source, snapshot, checksum, regeneration command, and a freshness/coverage
test. New concrete registry residues require notation, molecular round-trip,
PDB/mmCIF recognition, and rendering coverage.

The full external GlycoShape audit is optional during normal development:

```bash
CRABWURCS_GLYCOSHAPE_JSON=/path/to/GLYCOSHAPE.json \
  cargo test -p crabwurcs-core --test glycoshape_tests -- --ignored
```
