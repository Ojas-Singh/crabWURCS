# Changelog

## Unreleased

Drawings can be laid out for a figure, uncertainty is drawn the SNFG way, and
glycans can be weighed, searched and compared. Everything new is opt-in: the
default rendering of every GlycoShape record, and every conversion of one, is
unchanged byte for byte. The molecular fixes under "Output changes" affect
structures outside that corpus.

- **Orientation and compact drawings.** `RenderOptions` gains `orientation`
  (`RootRight`, the SNFG convention and the default, `RootLeft`, `RootBottom`,
  `RootTop`) and `compact` (tighter spacing and no linkage text, for heatmap
  axes and tables). The higher acceptor position of a branch stays above the
  lower one in horizontal drawings and to its right in vertical ones.
- **Reducing-end text.** `reducing_end` draws text such as `Asn`, `Ser/Thr` or
  `2-AB` past the reducing end.
- **Repeat brackets.** `repeat_brackets` draws a repeat unit in square brackets
  with its count (`n`, `n = 3`, `n = 2–5`, `n ≥ 2`) instead of a dashed bond
  closing it. A cap beyond the far end of the unit stays outside.
- **Floating fragments.** `floating_fragments` draws undefined antennae and
  substituents above the structure, over a bracket across their candidate
  parents, instead of a dashed bond to each candidate.
- **Positional and per-residue highlighting.** `MotifQuery` adds a position per
  motif residue (`MotifPosition`: `Terminal`, `Internal`, `Reducing`,
  `Flexible`), listed in the order the residues are read, as GlycoDraw's
  `highlight_termini_list`. `HighlightSelection::from_matches` and
  `from_residues` build a selection from motif occurrences or node indices.
- **A library of named motifs.** `named_motif`, `named_motifs`, `resolve_motif`
  and `annotate_motifs` know 27 motifs: Lewis and blood-group epitopes, LacNAc
  variants, N-glycan core features and O-glycan cores. Their positions keep
  `LewisX` out of sialyl-Lewis X and `CoreFucose` off an antenna. A motif string
  in `auto` format may now be one of these names.
- **Mass and composition.** `glycan_mass`, `glycan_formula` and
  `glycan_composition` give the formula, monoisotopic and average mass, m/z with
  an adduct, and a composition such as `HexNAc4Hex5Fuc1NeuAc2`, with reducing-end
  derivatives (`Reduced`, `TwoAb`, `TwoAa`, `Procainamide`) and permethylation.
  Residue formulas come from the molecular constructor; compositions, `?`
  linkages and undefined fragments have a mass although they have no molecular
  graph. Repeat ranges, probabilities and residues of unknown chemistry are a
  typed `MassError`, and an exact repeat is multiplied out. Phosphocholine is
  counted as its inner salt.
- **Fragment ions.** `glycan_fragments` lists the Domon–Costello B, C, Y and Z
  ions of every glycosidic cleavage, and `cleavage_mark` and
  `RenderOptions::cleavages` mark them on the drawing.
- **Structure comparison.** `structure_key` and `same_structure` compare
  glycans whatever notation or residue numbering they came from. A stated
  anomer or position is more than an open one.
- **Grids and figures.** `render_grid_svg` and `render_grid_png` draw many
  glycans on one canvas, and `svg_to_png` is public. Python's new
  `crabwurcs.figures` has `render_grid`, `write_images`, and `annotate_svg` and
  `annotate_figure`, which replace the glycan labels of any SVG or matplotlib
  figure with their drawings.
- **Command lines.** `render` gains `--orientation`, `--compact`,
  `--reducing-end`, `--repeat-brackets`, `--floating-fragments`,
  `--motif-positions`, `--highlight-residue` and `--cleavage`, and there are new
  commands `grid`, `mass`, `fragments`, `motif`, `motifs` and `same`.
- **Python.** `Glycan` gains `contains`, `count_motif`, `find_motif`, `motifs`,
  `residues`, `composition`, `composition_string`, `formula`, `mass` and
  `fragments`, and `Glycan.render` takes the new options as keyword arguments.
  Glycans now compare and hash by structure, not identity. New exceptions:
  `MassError` and `MotifError`.

### Fixes

- The Python command line writes SVG as UTF-8. `crabwurcs render --output` failed
  on Windows for any glycan with a Greek linkage label.

### Output changes

- **Alditols are open chains.** A backbone with no anomeric prefix (`h2122h`,
  `hxh`) was built with a phantom carbon and a ring: sorbitol came out as a
  seven-carbon ring and glycerol as a furanose. They now write as the polyols
  they are, and SMILES, MOL and SDF read back to the same WURCS.
- **Alditols are read from molecules.** A free polyol, or an alditol at the
  reducing end of a glycan such as `Gal(b1-3)GalNAc-ol`, is now a residue. Before,
  a free polyol was `NoGlycanFound` and a reducing-end alditol was silently dropped
  as an aglycone. A chain counts only if the glycan read accounts for the whole
  formula, so a glycerol lipid keeps the old reading. Both ends of a polyol are
  one molecule, so it is numbered like the sugar it comes from: an amide on the
  lowest carbon, then a D highest stereocentre. A truly symmetric polyol takes the
  lowest position, so `a6-b1` on mannitol reads back as the equivalent `a1-b1`.
- **Open-chain carbonyls keep their carbons.** `hO` and `AOd` backbones lost two
  carbons between them, which `hO122h` wrote as a five-carbon ring. They now write
  in full, and `o` and `O` write as the carbonyl they mean instead of a closed ring.
- **A bonded anomeric carbon has no hydroxyl.** An anomeric carbon kept its hydroxyl
  whatever else was bonded there, which gave it five bonds: trehalose and sucrose
  (bonded through it), a MAP bridge on it, and a substituent such as the `1*OC` of a
  methyl glycoside. Whatever is bonded there now replaces the hydroxyl, and the
  anomer is read through that atom. Trehalose, sucrose and methyl glucoside match
  their PubChem structures.
- **The sugar on an anomeric nitrogen is read.** A ring whose anomeric carbon
  carries nitrogen was skipped without an error, so an N-glycan on an asparagine
  lost its root GlcNAc: LacNAc on asparagine read as Gal. It is now a residue, its
  anomer read through the nitrogen, and the asparagine is an aglycone like any
  other. A bare amine stays a `1*N` substituent. The nitrogen of an aromatic or
  cyclic base does not count, so nucleosides and nucleotide sugars read as before.
  Bridges are still not read from molecules: what is past an anomeric nitrogen is
  an aglycone, so a MAP bridge there reads as the sugar alone.
- **Trehalose and sucrose are read.** Two anomeric carbons sharing one oxygen were
  counted as two bonds, which failed. They are one bond now, from the residue WURCS
  writes first, however the molecule is written. `MolError::AmbiguousGlycanFound`
  no longer says corpus glycans matched, which was only ever one of its causes.

### API notes for Rust callers

- `RenderOptions` has new fields: a struct literal needs `..Default::default()`.
- `MotifError` gains `PositionCount`, and `SnfgError` gains
  `InvalidResidueIndex` and `InvalidBondIndex`: an exhaustive `match` needs arms.

## 0.3.3 - 2026-10-15

PDB/mmCIF reading derives glycan chemistry from the 3D coordinates, so a
glycan reads as the molecule its atoms describe whatever its residue names
say. Several written strings change; they are listed under "Output changes".

- **Coordinates decide.** Every residue is searched for a monosaccharide ring
  (five or six atoms, one oxygen) in the bond graph, whatever it is called.
  The carbon numbering comes from the ring graph; the configuration of every
  stereocentre (so D/L and the parent sugar: Glc, Gal, Man, Ido, Gul, Alt,
  Tal, All, Ara, Lyx, Xyl, Rib, the ketoses, heptoses and nonulosonic acids),
  deoxy positions and uronic acids from signed volumes; the anomer from the
  anomeric centre against the anomeric reference atom (IUPAC 2-Carb-6.2);
  substituents (N-acetyl, N-glycolyl, N-sulfate, amino, O-sulfate, phosphate,
  phosphocholine, O-methyl, O-acetyl, fluoro) and glycosidic linkages from
  what the ring atoms are bonded to. No atom or residue name decides
  anything. A name only fills in what the coordinates cannot show: atoms
  missing (a crystal structure without hydrogens), a centre too flat to read
  (normalised signed volume below 0.15), a residue whose ring atoms are
  missing (read from its name, as before). wwPDB CCD coordinates of 49
  components fix the conventions (`crabwurcs-pdb/tests/data`, regenerated by
  the `regenerate_ccd_stereo_references` example).
- **Every disagreement is reported.** New `read_glycans_from_file` and
  `read_glycans_from_str` take `PdbReadOptions` and return a `PdbReading`:
  the glycans, a `PdbWarning` for every disagreement between a residue name
  and the coordinates (`residue 2uA A:4: named α-L-IdopA, coordinates say
  α-L-GulpA2S (C2)`) and for every part the coordinates could not decide,
  and a `ResidueReading` per residue (the UniqueRES read from the coordinates
  and the one the name describes, ring and backbone atoms, the Fischer side
  and normalised volume of every centre). `NamePolicy` chooses what decides:
  `Coordinates` (default), `Names` (a name decides where it says something),
  `Strict` (a disagreement is `PdbError::NamesContradictCoordinates`).
  Python: `extract_pdb(..., policy=...)`, `extract_pdb_file(..., policy=...)`,
  `ExtractedGlycan.warnings` and `.readings`. Command lines: `crabwurcs
  pdb-to-wurcs --policy coordinates|names|strict` prints the warnings to
  standard error (`--quiet` silences them). The `extract_glycans_*` functions
  read with the default policy and drop the warnings. `PdbError` gains the
  variant `NamesContradictCoordinates` (an exhaustive `match` needs an arm).
- **GLYCAM residue codes**: the case of the sugar letter is the configuration
  (upper case D, lower case L) for every sugar. The registry templates of
  IdoA, Alt and Sor are L, and their lower-case codes were mirrored: `0uA`,
  `4uA`, `YuA` (α-L-IdoA) read as α-D-IdoA, `0eA` as D-Alt, `0bA` as D-Sor.
  The Kdo code (`KO`) carries no anomer (it was read as α).
- **Sulfated glycans stay one glycan.** Sulfate, methyl, acetyl and
  phosphocholine in residues of their own (GLYCAM's `SO3`, `MEX`, `ACX`,
  `PCX`) are substituents of the residue they are bonded to, whatever they
  are called. Residues written with a GLYCAM letter code in `ATOM` records
  (`QYS`, `UYB`, `VMB`) and CCD components outside the bundled table (`RHM`)
  are recognised from their rings, so the glycans they joined are no longer
  split there.
- CCD components carry the substituents of their own atoms: `IDS` reads as
  IdoA2S (was IdoA), `SGN` as GlcNS6S (was GlcN), `NGY` as GlcNAc6S.
- Bonds are perceived within each atom's valence: inside a residue first,
  then the explicit `CONECT`/`LINK` bonds, then contacts between residues.
  A clash between residues (a sulfate packed 1.6 Å from the next ring oxygen
  in a crystal model of heparin, 3IRJ) no longer breaks a ring or a
  substituent. Metal ions are never bonded.
- mmCIF files with several models (NMR entries such as 1HPN) are read (the
  first model); they failed with "Atoms in Models not corresponding".
- Alternative locations: the atoms without one and those of the first
  alternative location are read (all were read together).
- The attachment site keeps its meaning (the residue an explicit `CONECT` or
  `LINK` record bonds the glycan to, otherwise the root residue) but is
  deterministic: the bond at the root's anomeric atom comes first, and a
  substituent residue such as `SO3` is never the site (one bonded residue was
  chosen arbitrarily).

### Output changes

- GLYCAM-named files read the configurations the atoms hold: `LIdopA` (was
  `DIdopA`), `LAltp`, `LSorp` (were `DAltp`, `DSorp`); the reducing-end
  anomer and every configuration are those of the coordinates, not of the
  residue names (GlycoShape Level 1 `beta.pdb` files whose reducing residue
  keeps an α code such as `4YS`, `0KN` or `NGC` read β; `NDG`-named
  Asn-linked GlcNAc reads β).
- Glycans that were split at `SO3`, `QYS`, `RHM` or GLYCAM-Web letter codes
  are one glycan, with its root's site (`A/ROH/1` where it was `A/QYS/2` or
  `A/SO3/7`).
- `IDS`, `SGN` and `NGY` carry their sulfates (`IdoA2S`, `GlcNS6S`,
  `GlcNAc6S`); `IDY` its 1-O-methyl.
- Phosphate monoesters are kept (`6*OPO/3O/3=O` in WURCS; they were dropped,
  and residues outside the CCD table such as ABF lost their stereochemistry,
  `Penf`). IUPAC and GLYCAM writing does not support phosphate yet, so
  `--to iupac-condensed` of such a glycan is an error rather than a name
  without its phosphate.
- A residue named as a sugar that does not match its coordinates reads as
  the coordinates say, with a warning: on the GlycoShape static build, 1,911
  of 1,982 Level 1 files read as their entry names (0.3.2: 1,836, with 71
  files split into several glycans). The other 71 belong to entries the
  entry audit reports: 66 files of 33 entries whose structure contradicts
  the name (2-O-sulfo-α-L-GulA under IdoA2S names, β reducing ends under α
  names, α-D-Galf and Araf anomers, α-GlcNAc on Man O6), three Level 1 files
  holding the other reducing-end anomer (GS00430, GS00954) and the two files
  of GS00039, which hold one sulfate more than its name.

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
