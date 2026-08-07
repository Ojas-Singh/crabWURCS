#!/usr/bin/env bash
set -euo pipefail

root_version=$(sed -n 's/^version = "\([^"]*\)"/\1/p' Cargo.toml | head -n 1)
python_version=$(sed -n 's/^version = "\([^"]*\)"/\1/p' pyproject.toml | head -n 1)
expected=${1:-$root_version}

test "$root_version" = "$expected"
test "$python_version" = "$expected"
grep -q "version = \"$expected\"" Cargo.lock
grep -q 'rust-version = "1.97"' Cargo.toml

# Cargo cannot fully package a path-dependent crate until the matching new
# dependency version exists on crates.io. Construct the dependency root now,
# and validate the exact publish file list for every dependent crate. The tag
# workflow's ordered `cargo publish` performs full package verification after
# each dependency becomes available.
cargo package -p crabwurcs-core --locked --allow-dirty >/dev/null
for package in crabwurcs-iupac crabwurcs-mol crabwurcs-pdb crabwurcs-snfg crabwurcs crabwurcs-cli; do
  cargo package -p "$package" --locked --allow-dirty --list >/dev/null
done

echo "release metadata and crate packages are consistent at $expected"
