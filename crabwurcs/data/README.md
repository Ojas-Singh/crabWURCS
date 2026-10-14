# GlycoShape notation pairs

`glycoshape_notations.tsv` is a deterministic projection of the repository's
`GLYCOSHAPE.json`. Its five columns are WURCS, IUPAC condensed, IUPAC extended,
GLYCAM, and isomeric SMILES. Only archetype records containing both WURCS and
SMILES are included.

The facade embeds these tables (with `glycoshape_derived_notations.tsv` and
`glycoshape_notation_only.tsv`). Conversions do not read them: they are
computed from the parsed glycan. The table is used on request only, as a
fallback for inputs that cannot be converted or as a check of a computed
structure (`convert_with_table`, `--table fallback|check`), and as the
offline regression oracle of the test suite. Its IUPAC and GLYCAM columns
keep GlycoShape's strings as stored, including a GlyTouCan accession in
place of a name for two records (G50508SG, G60371DN) and one damaged IUPAC
string (GS00430).
