import pytest

import crabwurcs
from crabwurcs import figures
from crabwurcs.cli import main

LACNAC = "Gal(b1-4)GlcNAc"
MAN5 = "Man(a1-3)[Man(a1-3)[Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc"
N_GLYCAN = (
    "Neu5Ac(a2-3)Gal(b1-4)GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-2)Man(a1-6)]"
    "Man(b1-4)GlcNAc(b1-4)[Fuc(a1-6)]GlcNAc"
)
REPEAT = (
    "WURCS=2.0/4,5,5/[a2112h-1b_1-5_2*NCC/3=O][a2122A-1b_1-5][a2112h-1a_1-5]"
    "[a2122h-1b_1-5]/1-2-1-3-4/a3-b1_b4-c1_c4-d1_d3-e1_a1-e4~2-5"
)
FRAGMENT = (
    "WURCS=2.0/6,11,10/[a2122h-1x_1-5_2*NCC/3=O][a2122h-1b_1-5_2*NCC/3=O]"
    "[a1122h-1b_1-5][a1122h-1a_1-5][a2112h-1b_1-5][a1221m-1a_1-5]/1-2-3-4-2-5-4-2-6-2-5/"
    "a4-b1_a6-i1_b4-c1_c3-d1_c6-g1_d2-e1_e4-f1_g2-h1_j4-k1_j1-d4|d6|g4|g6}"
)


def glycan(value=LACNAC):
    return crabwurcs.Glycan.parse(value)


# ── Existing behaviour stays as it was ──────────────────────────────────────


def test_rendering_without_new_arguments_is_the_plain_render():
    target = glycan(N_GLYCAN)
    plain = target.render()
    assert plain == target.render(compact=False, repeat_brackets=False)
    assert plain == target.render(orientation="root-right")
    assert plain == crabwurcs.render_snfg(N_GLYCAN)
    assert target.render(highlight_motifs=[]) == plain


def test_an_extracted_glycan_stays_hashable_and_comparable(tmp_path):
    pdb = "\n".join(
        [
            "HEADER    TEST",
            "HETATM    1 C1   NAG A   1      -1.000   0.000   0.000  1.00  0.00           C",
            "HETATM    2 O4   NAG A   1       0.000   0.000   0.000  1.00  0.00           O",
            "END",
            "",
        ]
    )
    first = crabwurcs.extract_pdb(pdb)[0]
    second = crabwurcs.extract_pdb(pdb)[0]
    assert first == second
    assert len({first, second}) == 1


# ── Drawing ─────────────────────────────────────────────────────────────────


@pytest.mark.parametrize(
    "orientation", ["root-right", "root-left", "root-bottom", "root-top", "top-down", "LR"]
)
def test_orientations_render_as_svg_and_png(orientation):
    target = glycan(MAN5)
    assert "<svg" in target.render(orientation=orientation)
    assert target.render("png", orientation=orientation).startswith(b"\x89PNG")


def test_a_vertical_drawing_is_taller_than_it_is_wide():
    svg = glycan("Gal(b1-4)GlcNAc(b1-2)Man(b1-4)GlcNAc").render(orientation="root-bottom")
    width = float(svg.split(' width="')[1].split('"')[0])
    height = float(svg.split(' height="')[1].split('"')[0])
    assert height > 2 * width


def test_an_unknown_orientation_is_a_parse_error():
    with pytest.raises(crabwurcs.ParseError):
        glycan().render(orientation="sideways")


def test_compact_drops_linkage_text():
    assert 'class="link"' in glycan().render()
    assert 'class="link"' not in glycan().render(compact=True)


def test_the_reducing_end_label_is_drawn_and_escaped():
    svg = glycan().render(reducing_end="<Asn>")
    assert ">&lt;Asn&gt;</text>" in svg


def test_repeats_and_fragments_are_drawn_the_snfg_way_on_request():
    repeat = glycan(REPEAT)
    assert 'class="bracket"' not in repeat.render()
    assert ">n = 2–5</text>" in repeat.render(repeat_brackets=True)
    fragment = glycan(FRAGMENT)
    assert "data-floating-band" not in fragment.render()
    assert "data-floating-band" in fragment.render(floating_fragments=True)


def test_scale_is_checked():
    with pytest.raises(crabwurcs.RenderError):
        glycan().render(scale=0)


# ── Highlighting ────────────────────────────────────────────────────────────


def test_motifs_highlight_by_name_and_by_position():
    capped = glycan("Neu5Ac(a2-3)Gal(b1-4)GlcNAc")
    assert 'class="motif-match"' in capped.render(highlight_motifs=["LacNAc"])
    # The Gal of this LacNAc is not terminal.
    terminal = capped.render(
        highlight_motifs=[LACNAC], highlight_positions=[("terminal", "flexible")]
    )
    assert 'class="motif-match"' not in terminal


def test_positions_can_be_left_open_per_motif():
    capped = glycan("Neu5Ac(a2-3)Gal(b1-4)GlcNAc")
    svg = capped.render(
        highlight_motifs=[LACNAC, "Neu5Ac(a2-3)Gal"],
        highlight_positions=[("terminal", "flexible"), None],
    )
    assert svg.count('class="motif-match"') == 3


def test_more_positions_than_motifs_is_an_error():
    with pytest.raises(crabwurcs.MotifError):
        glycan().render(highlight_motifs=[LACNAC], highlight_positions=[None, None])


def test_residues_highlight_by_node_index():
    svg = glycan().render(highlight_residues=[0])
    assert '<g data-node-index="0" class="motif-match"' in svg
    assert '<g data-node-index="1" class="motif-dimmed"' in svg
    with pytest.raises(crabwurcs.RenderError):
        glycan().render(highlight_residues=[9])


def test_residues_describe_the_tree():
    residues = glycan("Fuc(a1-3)[Gal(b1-4)]GlcNAc").residues()
    assert [residue.parent for residue in residues] == [None, 0, 0]
    assert [residue.terminal for residue in residues] == [False, True, True]
    assert residues[0].name == "GlcNAc"
    assert all(residue.recognized for residue in residues)
    assert not glycan("Foo").residues()[0].recognized


# ── Mass and composition ────────────────────────────────────────────────────


def test_mass_formula_and_composition():
    man5 = glycan(MAN5)
    assert man5.formula() == "C46H78N2O36"
    assert man5.mass() == pytest.approx(1234.4334, abs=0.0001)
    assert man5.mass(average=True) == pytest.approx(1235.08, abs=0.05)
    assert man5.mass(adduct="Na") == pytest.approx(1257.4226, abs=0.0003)
    assert man5.mass(adduct="[M+Na]+") == man5.mass(adduct="Na")
    assert man5.mass(adduct="H", charge=2) == pytest.approx(618.2240, abs=0.0001)
    assert man5.mass(permethylated=True, adduct="Na") == pytest.approx(1579.7826, abs=0.0001)
    assert man5.mass(reducing_end="2-ab") - man5.mass() == pytest.approx(120.0687, abs=0.001)
    assert man5.composition() == {"HexNAc": 2, "Hex": 5}
    assert man5.composition_string() == "HexNAc2Hex5"


def test_a_composition_has_a_mass_but_not_a_molecule():
    composition = glycan("{Hex}3,{HexNAc}2,{dHex}1")
    assert composition.mass() == pytest.approx(1056.3856, abs=0.001)
    with pytest.raises(crabwurcs.NonConcreteError):
        composition.to("smiles")


def test_what_leaves_the_mass_open_is_a_typed_error():
    with pytest.raises(crabwurcs.MassError):
        glycan("Foo").mass()
    with pytest.raises(crabwurcs.MassError):
        glycan(REPEAT).mass()
    with pytest.raises(crabwurcs.ParseError):
        glycan().mass(reducing_end="3-ab")
    with pytest.raises(crabwurcs.ParseError):
        glycan().mass(adduct="xenon")


def test_fragments_can_be_listed_and_marked():
    fragments = glycan("Gal(b1-4)GlcNAc(b1-2)Man").fragments()
    by_label = {fragment.label: fragment for fragment in fragments}
    assert by_label["B2"].mz == pytest.approx(366.1395, abs=0.0002)
    assert by_label["B1"].residues == (2,)
    assert by_label["Y1"].kind == "Y"
    svg = glycan("Gal(b1-4)GlcNAc(b1-2)Man").render(cleavages=["B2", "Y2"])
    assert svg.count('class="cleavage"') == 2
    with pytest.raises(crabwurcs.RenderError):
        glycan().render(cleavages=["B9"])
    with pytest.raises(crabwurcs.MassError):
        glycan(REPEAT).fragments()


# ── Comparison and search ───────────────────────────────────────────────────


def test_glycans_compare_as_structures():
    assert glycan("Fuc(a1-3)[Gal(b1-4)]GlcNAc") == glycan("Gal(b1-4)[Fuc(a1-3)]GlcNAc")
    assert glycan(LACNAC) == crabwurcs.Glycan.parse(glycan(LACNAC).to("wurcs"), "wurcs")
    assert glycan(LACNAC) == crabwurcs.Glycan.parse("DGalpb1-4DGlcpNAc", "glycam")
    assert glycan(LACNAC) != glycan("Gal(b1-3)GlcNAc")
    assert glycan("GlcNAc") != glycan("GlcNAc(b1-")
    assert glycan(LACNAC) != LACNAC
    assert len({glycan(LACNAC), crabwurcs.Glycan.parse("DGalpb1-4DGlcpNAc", "glycam")}) == 1


def test_a_glycan_no_notation_represents_equals_only_itself():
    assert glycan("Foo") == glycan("Foo")  # its IUPAC still identifies it
    assert hash(glycan("Foo")) == hash(glycan("Foo"))


def test_motif_search():
    target = glycan(N_GLYCAN)
    assert target.contains("Fuc(a1-6)GlcNAc")
    assert target.contains("CoreFucose")
    assert target.contains(glycan("Fuc(a1-6)GlcNAc"))
    assert not target.contains("LewisX")
    assert target.count_motif("GlcNAc(b1-2)Man") == 2
    assert target.count_motif("GlcNAc(b1-2)Man", positions=("terminal", "flexible")) == 1
    [occurrence] = target.find_motif("Fuc(a1-6)GlcNAc")
    assert len(occurrence.residues) == 2 and len(occurrence.bonds) == 1
    assert target.count_motif("Gal(b1-4)GlcNAc", "iupac-condensed") == 1
    with pytest.raises(crabwurcs.MotifError):
        target.contains(LACNAC, positions=("terminal",))
    with pytest.raises(crabwurcs.ParseError):
        target.contains(LACNAC, positions=("sideways", "flexible"))


def test_the_library_is_listed_and_annotates():
    library = {motif.name: motif for motif in crabwurcs.named_motifs()}
    assert "LewisX" in library
    assert library["LewisX"].iupac == "Gal(b1-4)[Fuc(a1-3)]GlcNAc"
    assert len(library["CoreFucose"].positions) == 2
    found = glycan(N_GLYCAN).motifs()
    assert {"NGlycanCore", "CoreFucose", "Sialyl3LacNAc"} <= set(found)
    assert "LewisX" not in found


# ── Figures ─────────────────────────────────────────────────────────────────


def test_grids_take_notations_or_glycans():
    svg = figures.render_grid([LACNAC, glycan(MAN5)], columns=1, labels=["LacNAc", "Man5"])
    assert svg.count("<svg ") == 3
    assert ">Man5</text>" in svg
    assert figures.render_grid([LACNAC], image_format="png").startswith(b"\x89PNG")
    with pytest.raises(ValueError):
        figures.render_grid([LACNAC], image_format="gif")
    with pytest.raises(TypeError):
        figures.render_grid([LACNAC], colour_scheme="dark")


def test_write_images_names_files_by_key(tmp_path):
    written = figures.write_images({"a": LACNAC, "b": glycan(MAN5)}, tmp_path, compact=True)
    assert sorted(path.name for path in written) == ["a.svg", "b.svg"]
    assert "<svg" in (tmp_path / "a.svg").read_text(encoding="utf-8")
    pngs = figures.write_images({"a": LACNAC}, tmp_path / "png", image_format="png")
    assert pngs[0].read_bytes().startswith(b"\x89PNG")


SVG_WITH_LABELS = """<svg xmlns="http://www.w3.org/2000/svg" width="400" height="300" viewBox="0 0 400 300">
<text x="100" y="50" style="font: 10px sans-serif; text-anchor: middle">Gal(b1-4)GlcNAc</text>
<text x="100" y="80" style="font: 10px sans-serif; text-anchor: end" transform="rotate(-90 100 80)">Man(a1-3)Man</text>
<text x="10" y="120" font-size="12">0.5</text>
<text x="10" y="150" font-size="12">Time (min)</text>
<text x="10" y="180" font-size="12">G0F</text>
</svg>
"""


def test_annotate_svg_replaces_only_labels_that_name_a_glycan():
    annotated = figures.annotate_svg(SVG_WITH_LABELS)
    assert "Gal(b1-4)GlcNAc</text>" not in annotated
    assert "Man(a1-3)Man</text>" not in annotated
    assert ">0.5</text>" in annotated
    assert ">Time (min)</text>" in annotated
    assert ">G0F</text>" in annotated
    assert annotated.count("<svg ") == 3
    # The rotation moves with the label, and no document metadata comes along.
    assert '<g transform="rotate(-90 100 80)">' in annotated
    assert "<metadata" not in annotated and "snfg-title" not in annotated
    assert figures.svg_to_png(annotated).startswith(b"\x89PNG")


def test_annotate_svg_anchors_a_drawing_like_its_text():
    annotated = figures.annotate_svg(SVG_WITH_LABELS, height=30)
    nested = annotated.split("<svg ")[2].split(">")[0]
    x = float(nested.split('x="')[1].split('"')[0])
    width = float(nested.split('width="')[1].split('"')[0])
    assert width > 0
    # Anchored in the middle at x=100.
    assert x + width / 2 == pytest.approx(100)


def test_annotated_drawings_share_one_scale():
    # A small glycan's symbols must not grow to fill the height of a big one.
    document = """<svg xmlns="http://www.w3.org/2000/svg" width="400" height="300" viewBox="0 0 400 300">
<text x="100" y="50" font-size="10">Gal(b1-4)GlcNAc</text>
<text x="100" y="150" font-size="10">Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc</text>
</svg>
"""
    annotated = figures.annotate_svg(document, height=60)
    small, large = annotated.split("<svg ")[2:4]
    size = lambda header, name: float(header.split(f'{name}="')[1].split('"')[0])
    assert size(large, "height") == pytest.approx(60)
    assert size(small, "height") < 60
    box = lambda header: [float(value) for value in header.split('viewBox="')[1].split('"')[0].split()]
    assert size(small, "height") / box(small)[3] == pytest.approx(size(large, "height") / box(large)[3])
    assert figures.annotate_svg("<svg></svg>") == "<svg></svg>"


def test_annotate_svg_takes_a_table_of_names():
    annotated = figures.annotate_svg(
        SVG_WITH_LABELS, labels={"G0F": MAN5, "Time (min)": glycan(LACNAC)}
    )
    assert ">G0F</text>" not in annotated
    assert ">Time (min)</text>" not in annotated
    # With a table, nothing else is touched.
    assert "Gal(b1-4)GlcNAc</text>" in annotated


def test_annotate_figure_replaces_tick_labels(tmp_path):
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")
    from matplotlib import pyplot

    figure, axes = pyplot.subplots()
    axes.bar([0, 1], [3, 5])
    axes.set_xticks([0, 1], [LACNAC, "Man(a1-3)Man"])
    axes.set_ylabel("Abundance")
    path = tmp_path / "figure.svg"
    svg = figures.annotate_figure(figure, path)
    assert svg.count("<svg ") == 3
    assert path.read_text(encoding="utf-8") == svg
    assert ">Abundance</text>" in svg
    png = tmp_path / "figure.png"
    figures.annotate_figure(figure, png)
    assert png.read_bytes().startswith(b"\x89PNG")
    pyplot.close(figure)


def test_rotated_labels_meet_their_tick():
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")
    from matplotlib import pyplot

    figure, axes = pyplot.subplots()
    axes.bar([0, 1], [3, 5])
    axes.set_xticks([0, 1], [LACNAC, "Man(a1-3)Man"], rotation=90)
    tick = axes.get_xticklabels()[0]
    before = (tick.get_rotation_mode(), tick.get_ha(), tick.get_va())
    svg = figures.annotate_figure(figure)
    # Anchored at the end of the rotated text, with the rotation about that
    # same point: the drawing's right edge is where the text ended.
    [group] = [part for part in svg.split('<g transform="rotate(')[1:] if "<svg" in part][:1]
    pivot_x = float(group.split(" ")[1])
    nested = group.split("<svg ")[1].split(">")[0]
    x = float(nested.split('x="')[1].split('"')[0])
    width = float(nested.split('width="')[1].split('"')[0])
    assert x + width == pytest.approx(pivot_x)
    # The figure is left as it was found.
    assert (tick.get_rotation_mode(), tick.get_ha(), tick.get_va()) == before
    pyplot.close(figure)


def test_a_figure_keeps_its_raster_images_as_png(tmp_path):
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")
    from matplotlib import pyplot
    from PIL import Image

    figure, axes = pyplot.subplots(figsize=(4, 3))
    axes.imshow([[0.0, 1.0], [1.0, 0.0]], cmap="gray", aspect="auto")
    axes.set_yticks([0, 1], [LACNAC, "Man(a1-3)Man"])
    axes.set_xticks([0, 1], ["WT", "KO"])
    figure.subplots_adjust(left=0.4)
    labels_before = [tick.get_visible() for tick in axes.get_yticklabels()]
    png = tmp_path / "heatmap.png"
    figures.annotate_figure(figure, png)
    # The figure is left as it was found.
    assert [tick.get_visible() for tick in axes.get_yticklabels()] == labels_before
    pyplot.close(figure)

    image = Image.open(png).convert("RGB")
    assert image.size == (4 * 192, 3 * 192)
    # The dark quadrant of the heatmap is there: rasterizing the SVG loses it.
    width, height = image.size
    dark = image.getpixel((int(width * 0.55), int(height * 0.3)))
    assert max(dark) < 60, dark
    # So are the drawings: a saturated symbol colour left of the axes.
    colours = image.crop((0, 0, int(width * 0.38), height)).getcolors(maxcolors=1 << 20)
    assert any(max(colour) - min(colour) > 120 for _, colour in colours)


# ── Command line ────────────────────────────────────────────────────────────


def test_render_takes_the_drawing_flags(tmp_path):
    output = tmp_path / "glycan.svg"
    arguments = [
        "render",
        "--orientation",
        "root-top",
        "--compact",
        "--reducing-end",
        "Asn",
        "--highlight-residue",
        "0",
        "--output",
        str(output),
        LACNAC,
    ]
    assert main(arguments) == 0
    svg = output.read_text()
    assert ">Asn</text>" in svg
    assert 'class="link"' not in svg
    assert 'data-node-index="0" class="motif-match"' in svg


def test_svg_files_are_written_as_utf8_on_every_platform(tmp_path):
    # Linkage labels are Greek, which a platform's default encoding may lack.
    output = tmp_path / "glycan.svg"
    assert main(["render", "--output", str(output), LACNAC]) == 0
    assert "\u03b2" in output.read_text(encoding="utf-8")


def test_mass_command(capsys):
    assert main(["mass", "--adduct", "Na", MAN5]) == 0
    out = capsys.readouterr().out
    assert "composition\tHexNAc2Hex5\n" in out
    assert "formula\tC46H78N2O36\n" in out
    assert "monoisotopic\t1234.4334\n" in out
    assert "mz Na x1\t1257.4226\n" in out
    assert main(["mass", "--adduct", "M-H", "--reducing-end", "reduced", MAN5]) == 0
    assert main(["mass", "Foo"]) == 1


def test_motif_commands(capsys):
    assert main(["motif", "LacNAc", N_GLYCAN]) == 0
    assert capsys.readouterr().out.strip() == "1"
    assert main(["motif", "--positions", "terminal,flexible", LACNAC, N_GLYCAN]) == 1
    capsys.readouterr()
    assert main(["motifs"]) == 0
    assert "LewisX\tGal(b1-4)[Fuc(a1-3)]GlcNAc" in capsys.readouterr().out
    assert main(["motifs", N_GLYCAN]) == 0
    assert "CoreFucose" in capsys.readouterr().out


def test_same_and_fragments_commands(capsys):
    assert main(["same", "Fuc(a1-3)[Gal(b1-4)]GlcNAc", "Gal(b1-4)[Fuc(a1-3)]GlcNAc"]) == 0
    assert capsys.readouterr().out.strip() == "true"
    assert main(["same", LACNAC, "Gal(b1-3)GlcNAc"]) == 1
    capsys.readouterr()
    assert main(["fragments", LACNAC]) == 0
    assert capsys.readouterr().out.startswith("B1\t163.0601\t")


def test_grid_command(tmp_path):
    source = tmp_path / "glycans.tsv"
    source.write_text(f"{LACNAC}\tLacNAc\n\n{MAN5}\n")
    output = tmp_path / "grid.svg"
    assert main(["grid", "--output", str(output), "--columns", "1", str(source)]) == 0
    svg = output.read_text(encoding="utf-8")
    assert svg.count("<svg ") == 3 and ">LacNAc</text>" in svg
