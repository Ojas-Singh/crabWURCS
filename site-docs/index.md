# crabWURCS 0.3

crabWURCS converts glycan notations, works with molecular structures, extracts
glycans from coordinate files, and produces SNFG graphics. The same pure-Rust
engine powers its Rust library, Python package, and command-line application.

```bash
pip install crabwurcs
```

```python
import crabwurcs

glycan = crabwurcs.Glycan.parse("Gal(b1-4)GlcNAc")
print(glycan.to("wurcs"))
```

The 0.3 release supports CPython 3.9+, Rust 1.97+, all 87 registry symbols,
and molecular/PDB interoperability for all 74 concrete registry residues.

Use [Getting started](getting-started.md) for practical examples or consult
the [formats and coverage](formats.md) page before processing ambiguous data.
