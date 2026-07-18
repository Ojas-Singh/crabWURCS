# GlycoShape notation pairs

`glycoshape_notations.tsv` is a deterministic projection of the repository's
`GLYCOSHAPE.json`. Its five columns are WURCS, IUPAC condensed, IUPAC extended,
GLYCAM, and isomeric SMILES. Only archetype records containing both WURCS and
SMILES are included.

The facade embeds this compact table so the 838 supplied chemical structures
convert without an RDKit installation. It does not claim to replace general
SMILES parsing or MolWURCS extraction for molecules outside this corpus.
