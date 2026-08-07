import crabwurcs
from crabwurcs.cli import main


IUPAC = "Gal(b1-4)GlcNAc"


def test_conversion_and_opaque_glycan_api():
    glycan = crabwurcs.Glycan.parse(IUPAC, crabwurcs.Format.IUPAC_CONDENSED)
    wurcs = glycan.to(crabwurcs.Format.WURCS)
    assert wurcs.startswith("WURCS=2.0/")
    assert crabwurcs.convert(wurcs, "iupac-condensed", "wurcs") == IUPAC
    assert "@" in glycan.to("smiles")


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
    assert results[0].glycan.to("iupac-condensed") == "GlcNAc"
    assert results[0].residues[0].chain == "A"
    assert results[0].residues[0].sequence_number == 1


def test_python_cli_conversion(capsys):
    assert main(["convert", "--to", "wurcs", IUPAC]) == 0
    assert capsys.readouterr().out.startswith("WURCS=2.0/")


def test_version_is_synchronized():
    assert crabwurcs.__version__ == "0.3.0"
