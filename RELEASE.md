# Release process

1. Update `CHANGELOG.md` and remove the release's `Unreleased` marker.
2. Set the same version in workspace metadata and `pyproject.toml`.
3. Run `./scripts/check_release.sh VERSION` and the complete test commands in
   `CONTRIBUTING.md` from a clean checkout.
4. Run the **TestPyPI candidate** workflow and install its wheel in a clean
   environment.
5. Configure repository secrets `PYPI_API_TOKEN` and
   `CARGO_REGISTRY_TOKEN`. The optional `testpypi` environment uses its own
   TestPyPI trusted-publisher relationship.
6. Create and push the signed tag `vVERSION`.

The tag workflow validates all targets, builds and install-tests CPython 3.9+
`abi3` wheels for manylinux x86_64/aarch64, macOS x86_64/arm64, and Windows
x86_64, publishes Rust crates in dependency order, publishes the Python wheel
set and sdist, and creates a GitHub release with SHA-256 checksums.

crates.io and PyPI releases are immutable. Never move or reuse a published
version tag; issue a patch release for corrections.
