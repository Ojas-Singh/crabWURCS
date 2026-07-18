# crabWURCS

A pure-Rust toolkit for glycan notation conversion, chemical-structure
interop, and SNFG rendering — one workspace covering what
[GlycanFormatConverter](https://gitlab.com/glycoinfo/glycanformatconverter),
[GlycanFormatConverter-cli](https://gitlab.com/glycoinfo/GlycanFormatConverter-cli),
[MolWURCS](https://gitlab.com/glycoinfo/molwurcs), and
[seq2snfg](https://gitlab.com/glycoinfo/seq2snfg) do today, all of which are
Java/Maven projects.

## Status

The notation and drawing pipeline is working end to end:

- WURCS 2.0 parsing is lossless for all 839 WURCS records in
  `GLYCOSHAPE.json`, including ambiguous linkage positions.
- IUPAC condensed, IUPAC extended, and GLYCAM parse into one shared residue
  graph. Of the 839 records that contain both IUPAC and WURCS, 614 have
  exact rooted residue/linkage signatures and 225 differ only because the
  WURCS record omits a reducing-end ring closure; there are zero hard
  structural mismatches. Undeclared rings remain `Unknown` rather than
  being guessed.
- All 839 extended-IUPAC records and all 943 GLYCAM records round-trip
  losslessly.
- The four archetypes without source WURCS or SMILES are no longer silently
  guessed: KDO and Bac use the reference WURCS descriptors, the generic
  diamino/trideoxy residue retains unknown stereochemistry explicitly, and
  `G60371DN` resolves through its authoritative GlyTouCan sequence. Thus all
  943 archetype notation records participate in direct conversion.
- The facade and CLI support format autodetection and conversion among
  WURCS, IUPAC condensed, IUPAC extended, GLYCAM, and 938
  isomeric-SMILES pairs from the GlycoShape corpus. Of these, 838 carry source
  WURCS and 100 have reproducibly generated WURCS from their supplied IUPAC.
  A pure-Rust molecular backend parses all 938 and canonicalizes them stably; equivalent
  SMILES atom orderings resolve by graph identity rather than exact text.
  Every included notation is tested against every other output family; where
  the source data has multiple valid branch-order spellings for one WURCS,
  conversion returns one of those verified spellings.
- Corpus structures convert to and from MOL and single-record SDF. Generated
  records include a canonical-SMILES provenance marker so crabWURCS retains
  complete atom-centred stereochemistry even though the current upstream MOL
  writer does not yet encode every SMILES chiral center as CTfile wedges. The
  molecular API and CLI also read multi-record SDF and emit one WURCS per
  record in input order.
- Finite defined WURCS graphs that do not occur in the corpus are constructed
  de novo as stereochemical molecular graphs. A forced audit constructs all
  938 corpus WURCS records without lookup, serializes and reparses their
  canonical SMILES, then recovers semantically identical WURCS for all 938.
  The same fallback powers edited WURCS → SMILES/MOL/SDF conversion.
- A verified corpus glycan can also be extracted from a single-bond aglycone
  attachment (for example, a methyl glycoside). Candidate cuts are
  stereochemically canonicalized, prefer the largest intact glycan, and reject
  internal cuts between two sugar-like rings.
- Previously unseen molecular graphs now enter a de-novo pure-Rust MolWURCS
  extractor. It detects pyranose and furanose rings, aldose, 2-ketose, and
  nonulosonic carbon chains, anomeric bridges, 1→n and 2→n topology, terminal
  hydroxymethyl/methyl/acid chemistry, N-acetyl/N-glycolyl/N-sulfate,
  O-sulfate/O-methyl/O-acetyl/phosphocholine MAP substituents, and WURCS 1/2
  stereodescriptors using MolWURCS's ligand-order parity convention.
  Unspecified stereocentres remain `x`. A forced lookup-free audit recovers
  semantically equivalent editable WURCS graphs for all 938 corpus molecules.
- WURCS compositions, whole-structure repeats, cyclic closures, and
  undefined fragment parent candidates are represented explicitly in the
  editable graph. Condensed/extended IUPAC composition, repeat, cycle, and
  reference-style fragment-anchor notation round-trip through that model.
- Common WURCS linkage and substituent probabilities are editable and
  convert to their distinct IUPAC forms (for example `1-55%4` and
  `6(?%)Me`). Terminal deoxy chemistry is kept separate from ring size.
- Common cross-linked WURCS MAP bridges, including phosphate and
  phosphoethanolamine attachment indices/directions, survive graph edits
  and convert through condensed IUPAC, extended IUPAC, and GLYCAM. SNFG
  bond labels retain the bridge identity.
- Undefined WURCS substituents such as `a?|b?}*OCC/3=O` are retained as
  candidate-parent modifications rather than disappearing after edits.
  SNFG draws converging dashed candidates to an uncertainty label; text
  exporters currently return an explicit error where no faithful IUPAC
  representation has been implemented.
- SNFG SVG rendering accepts any of those inputs and uses a collision-free
  tidy-tree layout, SNFG shapes/colors, linkage labels, and accessible SVG.
  Undefined antennae use dashed candidate bonds; compositions use compact
  counted symbols rather than dropping disconnected residues. All 938
  molecular-corpus records render without an unknown-symbol fallback; rare
  Alt, Gul, All, Tal, and Ido epimers use their standard SNFG families.

Broader MAP chemistry beyond the represented substituent families, polymeric
or ambiguous WURCS constructs that do not denote one finite molecule, and
sugar discrimination beyond the 938-structure audit remain MolWURCS porting
work. Non-glycans still return `NoGlycanFound`; the optional RDKit feature is
not needed for the implemented pure-Rust path.

Quick examples:

```bash
# Input format is detected automatically.
cargo run -p crabwurcs-cli -- convert --to wurcs \
  'Gal(b1-4)GlcNAc'

cargo run -p crabwurcs-cli -- convert --to glycam \
  'β-D-Galp-(1→4)-D-GlcNAc'

cargo run -p crabwurcs-cli -- render --output glycan.svg \
  'Gal(b1-4)GlcNAc'

# Composition and uncertain/cyclic structures are accepted too.
cargo run -p crabwurcs-cli -- convert --to wurcs \
  '{GlcNAc}2,{Man}3,{Fuc}1'
```

## Architecture

```
crabwurcs-core    — lossless WURCS parser/writer + shared ResidueGraph model
crabwurcs-iupac   — IUPAC condensed/extended and GLYCAM ⇄ ResidueGraph
crabwurcs-mol     — MOL/SDF/SMILES ⇄ ResidueGraph, pure Rust + optional RDKit extraction
crabwurcs-pdb     — glycan extraction from PDB/mmCIF structures (new scope, not a port)
crabwurcs-snfg    — ResidueGraph → SNFG SVG (seq2snfg equivalent)
crabwurcs         — facade crate re-exporting the above under one name
crabwurcs-cli     — the `crabwurcs` binary (replaces GlycanFormatConverter-cli)
```

Everything hangs off `crabwurcs_core::ResidueGraph`: residues are nodes,
definite glycosidic bonds are edges, and repeat/cycle/undefined-parent
metadata remains structural rather than being flattened into false bonds.
Every format converts through this model instead of maintaining an N-way
matrix of pairwise converters.

## Deliberate scope decisions (from the design conversation that produced this scaffold)

- **No GlycoCT.** Dropped by request; only WURCS and IUPAC condensed/extended
  are planned notation formats.
- **Pure Rust by default, RDKit over OpenBabel for future general extraction.**
  `chematic` supplies dependency-free SMILES graph canonicalization and
  MOL/SDF I/O. The optional `rdkit-backend` remains off by default. RDKit
  is BSD-3-Clause (OpenBabel is GPL-2.0, which would pull a
  statically-linked crabWURCS binary under GPL), and RDKit's ring
  perception/aromaticity/valence handling is the more complete primitive
  set for porting MolWURCS's sugar-discrimination algorithm onto. See
  `crabwurcs-mol/README-BACKEND.md` for build requirements, including why
  ARM64/Ampere is a well-supported target here, not an edge case.
- **`pdbtbx`** for PDB/mmCIF parsing — pure Rust, no C/C++ dependency.
- **License: dual MIT/Apache-2.0** (the Rust ecosystem convention), which
  is why the RDKit-vs-OpenBabel choice above mattered. If the ported
  algorithms end up transliterated closely from the GPL-2.0-or-later
  Java sources (GlycanFormatConverter, MolWURCS) rather than reimplemented
  from their published algorithm descriptions, that license choice needs
  revisiting — see each crate's module docs for pointers to the papers to
  implement from instead of the Java source directly.

## Toolchain note

The workspace declares **Rust 1.88** as its minimum version. This is required
by the exact-pinned pure-Rust chemistry backend's use of stabilized let-chain
syntax. The full workspace is currently developed and tested on Rust 1.97.

## Remaining major work

1. Extend the bidirectional lookup-free MolWURCS audit beyond the 938 verified
   GlycoShape molecules and add further uncommon MAP templates.
2. Add IUPAC representations for undefined substituent fragments, plus
   endpoint-specific undefined-fragment probabilities, uncommon MAP
   templates, and complex nested/partial repeat semantics.
3. Extend the SNFG residue dictionary and canonical branch ordering against
   the complete SNFG and GlycanFormatConverter fixtures.
4. Implement robust PDB/mmCIF linkage reconstruction.

## License

Licensed under either of [Apache License, Version 2.0](LICENSE-APACHE) or
[MIT license](LICENSE-MIT) at your option.
