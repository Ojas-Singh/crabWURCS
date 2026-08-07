# Getting started

## Python

```python
from pathlib import Path
import crabwurcs

glycan = crabwurcs.Glycan.parse(
    "Neu5Ac(a2-3)Gal(b1-4)GlcNAc",
    crabwurcs.Format.IUPAC_CONDENSED,
)

Path("glycan.svg").write_text(glycan.render("svg"))
Path("glycan.mol").write_text(glycan.to("mol"))
```

`Glycan` is immutable and can be serialized repeatedly. `convert()` and
`render_snfg()` are shorter alternatives when an intermediate object is not
needed.

## Rust

```toml
[dependencies]
crabwurcs = "0.3.0"
```

```rust
let graph = crabwurcs::parse_notation(
    "Gal(b1-4)GlcNAc",
    crabwurcs::Format::IupacCondensed,
)?;
let smiles = crabwurcs::write_notation(&graph, crabwurcs::Format::Smiles)?;
```

## Files versus text

Notation and molecule methods accept text. PDB/mmCIF offers separate
`extract_pdb(text)` and `extract_pdb_file(path)` functions so a notation can
never be mistaken for a path.
