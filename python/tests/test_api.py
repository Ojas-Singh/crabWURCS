import crabwurcs
from crabwurcs.cli import main


IUPAC = "Gal(b1-4)GlcNAc"
CANONICAL_WURCS = "WURCS=2.0/2,2,1/[u2122h][a2112h-1b_1-5]/1-2/a3-b1"
REORDERED_WURCS = "WURCS=2.0/2,2,1/[a2112h-1b_1-5][u2122h]/1-2/b3-a1"


def test_conversion_and_opaque_glycan_api():
    glycan = crabwurcs.Glycan.parse(IUPAC, crabwurcs.Format.IUPAC_CONDENSED)
    wurcs = glycan.to(crabwurcs.Format.WURCS)
    assert wurcs.startswith("WURCS=2.0/")
    assert crabwurcs.convert(wurcs, "iupac-condensed", "wurcs") == IUPAC
    assert "@" in glycan.to("smiles")


def test_wurcs_is_normalized_by_default_and_explicitly():
    assert crabwurcs.normalize_wurcs(REORDERED_WURCS) == CANONICAL_WURCS
    assert crabwurcs.convert(REORDERED_WURCS, "wurcs", "wurcs") == CANONICAL_WURCS
    assert crabwurcs.Glycan.parse(REORDERED_WURCS, "wurcs").to("wurcs") == CANONICAL_WURCS


def test_explicit_normalization_has_a_typed_unsupported_error():
    unsupported = "WURCS=2.0/1,1,0/[a2122h-1a_1-5_2|3*OC]/1/"
    try:
        crabwurcs.normalize_wurcs(unsupported)
    except crabwurcs.NormalizationError as error:
        assert isinstance(error, crabwurcs.ConversionError)
    else:
        raise AssertionError("expected NormalizationError")

    assert crabwurcs.Glycan.parse(unsupported, "wurcs").to("wurcs") == unsupported


def test_svg_png_and_motif_rendering():
    glycan = crabwurcs.Glycan.parse(IUPAC)
    svg = glycan.render(highlight_motifs=["Gal(b1-?)GlcNAc"])
    png = glycan.render("png")
    assert "<svg" in svg
    assert "motif-match" in svg
    assert png.startswith(b"\x89PNG\r\n\x1a\n")


def test_error_hierarchy_is_typed():
    try:
        crabwurcs.Glycan.parse("C1CC1", "smiles")
    except crabwurcs.ParseError as error:
        assert isinstance(error, crabwurcs.CrabWurcsError)
    else:
        raise AssertionError("expected ParseError")

    generic = crabwurcs.Glycan.parse("{Hex}2")
    try:
        generic.to("smiles")
    except crabwurcs.NonConcreteError as error:
        assert isinstance(error, crabwurcs.ConversionError)
    else:
        raise AssertionError("expected NonConcreteError")


def test_pdb_extraction_exposes_immutable_provenance(tmp_path):
    pdb = "\n".join(
        [
            "HEADER    TEST",
            "HETATM    1 C1   NAG A   1      -1.000   0.000   0.000  1.00  0.00           C",
            "HETATM    2 O4   NAG A   1       0.000   0.000   0.000  1.00  0.00           O",
            "END",
            "",
        ]
    )
    path = tmp_path / "glycan.pdb"
    path.write_text(pdb)
    results = crabwurcs.extract_pdb_file(path)
    assert len(results) == 1
    # The CCD code NAG declares the beta anomer of the reducing end.
    assert results[0].glycan.to("iupac-condensed") == "GlcNAc(b1-"
    assert results[0].residues[0].chain == "A"
    assert results[0].residues[0].sequence_number == 1


def _residue(atoms, name):
    return "\n".join(
        f"HETATM{serial:5d}  {atom:<3} {name:>3} A   1    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {atom[0]:>2}"
        for serial, (atom, x, y, z) in enumerate(atoms, 1)
    ) + "\nEND\n"


# alpha-L-gulopyranuronic acid (wwPDB CCD LGU, ideal coordinates, heavy atoms).
GULA = [
    ("C1", -1.185, 1.216, 0.417),
    ("C2", -1.774, -0.158, 0.746),
    ("O2", -3.192, -0.125, 0.569),
    ("C3", -1.165, -1.204, -0.193),
    ("O3", -1.525, -0.902, -1.543),
    ("C4", 0.359, -1.176, -0.050),
    ("O4", 0.721, -1.542, 1.283),
    ("C5", 0.866, 0.237, -0.348),
    ("O5", 0.238, 1.164, 0.540),
    ("C6", 2.360, 0.287, -0.152),
    ("O6B", 3.166, -0.429, -0.951),
    ("O6A", 2.832, 0.973, 0.724),
    ("O1", -1.536, 1.577, -0.920),
]


def test_pdb_residues_are_read_from_coordinates_and_names_that_disagree_are_reported(tmp_path):
    # GLYCAM names alpha-L-IdoA `2uA`; these coordinates are alpha-L-GulA.
    pdb = _residue(GULA, "2uA")
    [result] = crabwurcs.extract_pdb(pdb)
    assert result.glycan.to("wurcs") == "WURCS=2.0/1,1,0/[a1121A-1a_1-5]/1/"
    [warning] = result.warnings
    assert warning.kind == "name-contradicts-coordinates"
    assert str(warning) == "residue 2uA A:1: named α-L-IdopA, coordinates say α-L-GulpA (C2)"
    assert warning.positions == (2,)
    [reading] = result.readings
    assert reading.from_coordinates == "a1121A-1a_1-5"
    assert reading.from_name == "a2121A-1a_1-5"
    assert [centre.atom for centre in reading.centres] == ["C1", "C2", "C3", "C4", "C5"]

    [named] = crabwurcs.extract_pdb(pdb, policy="names")
    assert named.glycan.to("wurcs") == "WURCS=2.0/1,1,0/[a2121A-1a_1-5]/1/"
    assert len(named.warnings) == 1
    try:
        crabwurcs.extract_pdb(pdb, policy="strict")
    except crabwurcs.PdbError as error:
        assert "coordinates say α-L-GulpA (C2)" in str(error)
    else:
        raise AssertionError("expected PdbError")
    try:
        crabwurcs.extract_pdb(pdb, policy="bogus")
    except crabwurcs.ParseError:
        pass
    else:
        raise AssertionError("expected ParseError")

    path = tmp_path / "gula.pdb"
    path.write_text(pdb)
    assert crabwurcs.extract_pdb_file(path)[0].warnings == result.warnings


def test_python_cli_conversion(capsys):
    assert main(["convert", "--to", "wurcs", IUPAC]) == 0
    assert capsys.readouterr().out.startswith("WURCS=2.0/")


def test_conversion_computes_and_the_table_is_used_on_request(capsys):
    # A GlyTouCan accession stored as a GlycoShape name is not a notation.
    try:
        crabwurcs.convert("G60371D-N", "wurcs")
    except crabwurcs.ParseError:
        pass
    else:
        raise AssertionError("expected ParseError")
    resolved = crabwurcs.convert("G60371D-N", "wurcs", table="fallback")
    assert resolved.startswith("WURCS=2.0/8,23,22/")
    name = "LFucpa1-2DGalpa1-3DGalpb1-3DGalpNAc"
    assert crabwurcs.convert(name, "wurcs", "glycam", table="check") == crabwurcs.convert(
        name, "wurcs", "glycam"
    )
    assert main(["convert", "--to", "wurcs", "--table", "fallback", "G60371D-N"]) == 0
    assert capsys.readouterr().out.startswith("WURCS=2.0/8,23,22/")
    assert main(["convert", "--to", "wurcs", "G60371D-N"]) == 1


def test_written_names_keep_ring_size_configuration_and_reducing_end():
    alpha = crabwurcs.convert("LArafb1-2LArafa1-OH", "wurcs", "glycam")
    assert "[a211h-1a_1-4]" in alpha
    assert crabwurcs.convert(alpha, "iupac-condensed", "wurcs") == "Araf(b1-2)Araf(a1-"
    assert crabwurcs.convert(alpha, "glycam", "wurcs") == "LArafb1-2LArafa1-OH"
    d_fuc = "WURCS=2.0/1,1,0/[u2112m]/1/"
    assert crabwurcs.convert(d_fuc, "iupac-condensed", "wurcs") == "D-Fuc"
    acetylated = crabwurcs.convert("DNeup5Ac[9A]a2-6DGalpb1-OH", "wurcs", "glycam")
    assert "9*OCC/3=O" in acetylated


def test_python_cli_normalization(capsys):
    assert main(["normalize-wurcs", REORDERED_WURCS]) == 0
    assert capsys.readouterr().out.strip() == CANONICAL_WURCS


def test_version_is_synchronized():
    assert crabwurcs.__version__ == "0.3.2"
