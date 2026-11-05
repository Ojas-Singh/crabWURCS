"""Command-line interface installed by the Python distribution."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Optional, Sequence

from . import (
    Format,
    Glycan,
    __version__,
    convert,
    extract_pdb_file,
    named_motifs,
    render_snfg,
)


def _text(value: Optional[str], input_file: bool = False) -> str:
    if value is None:
        return sys.stdin.read()
    return Path(value).read_text(encoding="utf-8") if input_file else value


def _add_drawing(command: argparse.ArgumentParser) -> None:
    command.add_argument("--orientation")
    command.add_argument("--compact", action="store_true")
    command.add_argument("--reducing-end", dest="reducing_end_label")
    command.add_argument("--repeat-brackets", action="store_true")
    command.add_argument("--floating-fragments", action="store_true")


def _drawing(args: argparse.Namespace) -> dict:
    return {
        "orientation": args.orientation,
        "compact": args.compact,
        "reducing_end": args.reducing_end_label,
        "repeat_brackets": args.repeat_brackets,
        "floating_fragments": args.floating_fragments,
    }


def _add_mass(command: argparse.ArgumentParser) -> None:
    command.add_argument("input", nargs="?")
    command.add_argument("--input-file", action="store_true")
    command.add_argument("--from", dest="from_format", default="auto")
    command.add_argument("--reducing-end", default="free")
    command.add_argument("--permethylated", action="store_true")


def _positions(value: Optional[str]) -> Optional[Sequence[str]]:
    return value.split(",") if value else None


def _write_image(output: "str | bytes", path: Optional[str]) -> None:
    if path:
        target = Path(path)
        if isinstance(output, bytes):
            target.write_bytes(output)
        else:
            # Linkage labels are Greek, whatever the platform's default encoding.
            target.write_text(output, encoding="utf-8")
    elif isinstance(output, bytes):
        sys.stdout.buffer.write(output)
    else:
        print(output)


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="crabwurcs", description="Convert, inspect, and render glycans")
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    commands = parser.add_subparsers(dest="command", required=True)

    conversion = commands.add_parser("convert", help="convert glycan formats")
    conversion.add_argument("input", nargs="?")
    conversion.add_argument("--input-file", action="store_true")
    conversion.add_argument("--from", dest="from_format", default="auto")
    conversion.add_argument("--to", dest="to_format", required=True)
    conversion.add_argument("--table", default="off", choices=("off", "fallback", "check"))

    normalization = commands.add_parser("normalize-wurcs", help="normalize WURCS 2.0 text")
    normalization.add_argument("input", nargs="?")
    normalization.add_argument("--input-file", action="store_true")

    mol_to = commands.add_parser("mol-to-wurcs", help="extract WURCS from MOL/SDF/SMILES")
    mol_to.add_argument("input", nargs="?")
    mol_to.add_argument("--format", required=True, choices=("mol", "sdf", "smiles"))

    wurcs_to = commands.add_parser("wurcs-to-mol", help="render WURCS as MOL/SDF/SMILES")
    wurcs_to.add_argument("input", nargs="?")
    wurcs_to.add_argument("--format", required=True, choices=("mol", "sdf", "smiles"))

    pdb = commands.add_parser("pdb-to-wurcs", help="extract glycans from PDB/mmCIF")
    pdb.add_argument("input")
    pdb.add_argument("--to", dest="to_format", default="wurcs")
    pdb.add_argument("--policy", default="coordinates", choices=("coordinates", "names", "strict"))
    pdb.add_argument("--quiet", action="store_true", help="do not print warnings")

    render = commands.add_parser("render", help="render an SNFG SVG or PNG")
    render.add_argument("input", nargs="?")
    render.add_argument("--input-file", action="store_true")
    render.add_argument("--from", dest="from_format", default="auto")
    render.add_argument("--output")
    render.add_argument("--highlight-motif", action="append", default=[])
    render.add_argument("--motif-from", default="auto")
    render.add_argument("--motif-positions", action="append", default=[])
    render.add_argument("--highlight-residue", action="append", default=[], type=int)
    render.add_argument("--cleavage", action="append", default=[])
    _add_drawing(render)

    grid = commands.add_parser("grid", help="draw many glycans on one canvas")
    grid.add_argument("input")
    grid.add_argument("--output", "-o", required=True)
    grid.add_argument("--from", dest="from_format", default="auto")
    grid.add_argument("--columns", type=int, default=4)
    _add_drawing(grid)

    mass = commands.add_parser("mass", help="print composition, formula and masses")
    _add_mass(mass)
    mass.add_argument("--adduct", action="append", default=[])
    mass.add_argument("--charge", type=int, default=1)

    fragments = commands.add_parser("fragments", help="list B, C, Y and Z ions")
    _add_mass(fragments)
    fragments.add_argument("--adduct", default="H")

    motif = commands.add_parser("motif", help="count the occurrences of a motif")
    motif.add_argument("motif")
    motif.add_argument("input", nargs="?")
    motif.add_argument("--input-file", action="store_true")
    motif.add_argument("--from", dest="from_format", default="auto")
    motif.add_argument("--motif-from", default="auto")
    motif.add_argument("--positions")

    motifs = commands.add_parser("motifs", help="list the motif library or a glycan's motifs")
    motifs.add_argument("input", nargs="?")
    motifs.add_argument("--input-file", action="store_true")
    motifs.add_argument("--from", dest="from_format", default="auto")

    same = commands.add_parser("same", help="compare two glycans in any notations")
    same.add_argument("left")
    same.add_argument("right")
    same.add_argument("--from", dest="from_format", default="auto")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = _parser().parse_args(argv)
    try:
        if args.command == "convert":
            print(
                convert(
                    _text(args.input, args.input_file),
                    args.to_format,
                    args.from_format,
                    table=args.table,
                )
            )
        elif args.command == "normalize-wurcs":
            from . import normalize_wurcs

            print(normalize_wurcs(_text(args.input, args.input_file)))
        elif args.command == "mol-to-wurcs":
            print(convert(_text(args.input, bool(args.input)), Format.WURCS, args.format))
        elif args.command == "wurcs-to-mol":
            print(convert(_text(args.input, bool(args.input)), args.format, Format.WURCS))
        elif args.command == "pdb-to-wurcs":
            for extracted in extract_pdb_file(args.input, policy=args.policy):
                if not args.quiet:
                    for warning in extracted.warnings:
                        print(f"warning: {warning}", file=sys.stderr)
                notation = extracted.glycan.to(args.to_format)
                print(f"{extracted.attachment_site}\t{notation}" if extracted.attachment_site else notation)
        elif args.command == "render":
            text = _text(args.input, args.input_file)
            image_format = Path(args.output).suffix.lstrip(".").lower() if args.output else "svg"
            output = render_snfg(
                text,
                from_format=args.from_format,
                image_format=image_format,
                highlight_motifs=args.highlight_motif,
                motif_format=args.motif_from,
                highlight_positions=[_positions(value) for value in args.motif_positions] or None,
                highlight_residues=args.highlight_residue or None,
                cleavages=args.cleavage or None,
                **_drawing(args),
            )
            _write_image(output, args.output)
        elif args.command == "grid":
            from .figures import render_grid

            notations = []
            labels = []
            for line in Path(args.input).read_text(encoding="utf-8").splitlines():
                if line.strip():
                    notation, _, label = line.partition("\t")
                    notations.append(notation.strip())
                    labels.append(label.strip())
            _write_image(
                render_grid(
                    notations,
                    columns=args.columns,
                    labels=labels,
                    image_format=Path(args.output).suffix.lstrip(".").lower(),
                    from_format=args.from_format,
                    **_drawing(args),
                ),
                args.output,
            )
        elif args.command == "mass":
            glycan = Glycan.parse(_text(args.input, args.input_file).strip(), args.from_format)
            derivative = {"reducing_end": args.reducing_end, "permethylated": args.permethylated}
            print(f"composition\t{glycan.composition_string()}")
            print(f"formula\t{glycan.formula(**derivative)}")
            print(f"monoisotopic\t{glycan.mass(**derivative):.4f}")
            print(f"average\t{glycan.mass(average=True, **derivative):.4f}")
            for adduct in args.adduct:
                value = glycan.mass(adduct=adduct, charge=args.charge, **derivative)
                print(f"mz {adduct} x{args.charge}\t{value:.4f}")
        elif args.command == "fragments":
            glycan = Glycan.parse(_text(args.input, args.input_file).strip(), args.from_format)
            for fragment in glycan.fragments(
                adduct=args.adduct,
                reducing_end=args.reducing_end,
                permethylated=args.permethylated,
            ):
                print(f"{fragment.label}\t{fragment.mz:.4f}\t{fragment.formula}")
        elif args.command == "motif":
            glycan = Glycan.parse(_text(args.input, args.input_file).strip(), args.from_format)
            count = glycan.count_motif(
                args.motif, args.motif_from, positions=_positions(args.positions)
            )
            print(count)
            if count == 0:
                print(f"error: the glycan does not contain {args.motif}", file=sys.stderr)
                return 1
        elif args.command == "motifs":
            if args.input is None and not args.input_file:
                for motif in named_motifs():
                    print(f"{motif.name}\t{motif.iupac}")
            else:
                glycan = Glycan.parse(_text(args.input, args.input_file).strip(), args.from_format)
                for name in glycan.motifs():
                    print(name)
        elif args.command == "same":
            same = Glycan.parse(args.left, args.from_format) == Glycan.parse(
                args.right, args.from_format
            )
            print("true" if same else "false")
            if not same:
                print("error: the glycans differ", file=sys.stderr)
                return 1
        return 0
    except Exception as error:
        print(f"error: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
