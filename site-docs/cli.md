# Command line

```text
crabwurcs convert [--from FORMAT] --to FORMAT [--table off|fallback|check] [INPUT]
crabwurcs normalize-wurcs [--input-file] [INPUT]
crabwurcs mol-to-wurcs --format {mol,sdf,smiles} [FILE]
crabwurcs wurcs-to-mol --format {mol,sdf,smiles} [FILE]
crabwurcs pdb-to-wurcs [--to FORMAT] FILE
crabwurcs render [--from FORMAT] [--output FILE] [INPUT]
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
