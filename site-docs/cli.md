# Command line

```text
crabwurcs convert [--from FORMAT] --to FORMAT [--table off|fallback|check] [INPUT]
crabwurcs normalize-wurcs [--input-file] [INPUT]
crabwurcs mol-to-wurcs --format {mol,sdf,smiles} [FILE]
crabwurcs wurcs-to-mol --format {mol,sdf,smiles} [FILE]
crabwurcs pdb-to-wurcs [--to FORMAT] [--policy coordinates|names|strict] [--quiet] FILE
crabwurcs render [--from FORMAT] [--output FILE] [INPUT]
crabwurcs grid --output FILE [--columns N] FILE
crabwurcs mass [--reducing-end END] [--permethylated] [--adduct ADDUCT] [--charge N] [INPUT]
crabwurcs fragments [--reducing-end END] [--permethylated] [--adduct ADDUCT] [INPUT]
crabwurcs motif [--positions POSITIONS] MOTIF [INPUT]
crabwurcs motifs [INPUT]
crabwurcs same LEFT RIGHT
```

Omitted notation input is read from standard input. `convert` and `render`
accept `--input-file` when their positional argument is a path. Render output
is SVG by default; `.svg` and `.png` output suffixes select file type.

Both the Cargo and pip distributions provide the same command name and core
behavior, so install only one in a given environment.

`convert` computes every conversion from the parsed glycan. The bundled table
of GlycoShape notations is used only on request: `--table fallback` answers an
input that cannot be converted (for example a GlyTouCan accession stored as a
name) from its table record, and `--table check` fails when the input's
record holds another structure.

Every command that produces WURCS uses the offline normalizer. The dedicated
`normalize-wurcs` command is strict and exits unsuccessfully rather than
claiming an unsupported input was normalized.

`render` and `grid` take `--orientation`, `--compact`, `--reducing-end TEXT`,
`--repeat-brackets` and `--floating-fragments`; see
[SNFG rendering](rendering.md). `render` also takes `--highlight-motif`,
`--motif-positions terminal,flexible` (once per motif, in the same order),
`--highlight-residue INDEX`, and `--cleavage ION`.

`grid` reads one notation per line, with an optional tab and caption.

`mass` prints the composition, formula and masses, and `fragments` one line per
ion; see [Mass and fragments](mass.md). `motif` prints how often a motif, given
as a notation or a library name, occurs, and fails if it does not. `motifs` lists
the library, or the library motifs of a glycan. `same` prints whether two glycans
are the same structure, and fails if they are not.

`pdb-to-wurcs` reads every monosaccharide from its coordinates and prints one
line per glycan (`site<TAB>notation`). Where a residue name disagrees with the
coordinates, or the coordinates cannot decide part of a residue, a
`warning: ...` line goes to standard error (`--quiet` silences them).
`--policy names` lets names decide where they say something; `--policy
strict` fails on a disagreement. See [PDB and mmCIF extraction](pdb.md).
