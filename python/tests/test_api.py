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
    assert crabwurcs.__version__ == "0.3.1"
