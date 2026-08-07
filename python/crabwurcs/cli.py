"""Command-line interface installed by the Python distribution."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Optional, Sequence

from . import Format, __version__, convert, extract_pdb_file, render_snfg


def _text(value: Optional[str], input_file: bool = False) -> str:
    if value is None:
        return sys.stdin.read()
    return Path(value).read_text() if input_file else value


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="crabwurcs", description="Convert, inspect, and render glycans")
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    commands = parser.add_subparsers(dest="command", required=True)

    conversion = commands.add_parser("convert", help="convert glycan formats")
    conversion.add_argument("input", nargs="?")
    conversion.add_argument("--input-file", action="store_true")
    conversion.add_argument("--from", dest="from_format", default="auto")
    conversion.add_argument("--to", dest="to_format", required=True)

    mol_to = commands.add_parser("mol-to-wurcs", help="extract WURCS from MOL/SDF/SMILES")
    mol_to.add_argument("input", nargs="?")
    mol_to.add_argument("--format", required=True, choices=("mol", "sdf", "smiles"))

    wurcs_to = commands.add_parser("wurcs-to-mol", help="render WURCS as MOL/SDF/SMILES")
    wurcs_to.add_argument("input", nargs="?")
    wurcs_to.add_argument("--format", required=True, choices=("mol", "sdf", "smiles"))

    pdb = commands.add_parser("pdb-to-wurcs", help="extract glycans from PDB/mmCIF")
    pdb.add_argument("input")
    pdb.add_argument("--to", dest="to_format", default="wurcs")

    render = commands.add_parser("render", help="render an SNFG SVG or PNG")
    render.add_argument("input", nargs="?")
    render.add_argument("--input-file", action="store_true")
    render.add_argument("--from", dest="from_format", default="auto")
    render.add_argument("--output")
    render.add_argument("--highlight-motif", action="append", default=[])
    render.add_argument("--motif-from", default="auto")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = _parser().parse_args(argv)
    try:
        if args.command == "convert":
            print(convert(_text(args.input, args.input_file), args.to_format, args.from_format))
        elif args.command == "mol-to-wurcs":
            print(convert(_text(args.input, bool(args.input)), Format.WURCS, args.format))
        elif args.command == "wurcs-to-mol":
            print(convert(_text(args.input, bool(args.input)), args.format, Format.WURCS))
        elif args.command == "pdb-to-wurcs":
            for extracted in extract_pdb_file(args.input):
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
            )
            if args.output:
                path = Path(args.output)
                path.write_bytes(output) if isinstance(output, bytes) else path.write_text(output)
            elif isinstance(output, bytes):
                sys.stdout.buffer.write(output)
            else:
                print(output)
        return 0
    except Exception as error:
        print(f"error: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
