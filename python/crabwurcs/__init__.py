"""Python interface to crabWURCS.

The public API deliberately wraps the native Rust objects so callers get
stable enums, immutable result records, and useful type annotations.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from os import PathLike
from pathlib import Path
from typing import Iterable, List, Optional, Sequence, Union

from . import _native

__version__ = _native.__version__

CrabWurcsError = _native.CrabWurcsError
ParseError = _native.ParseError
ConversionError = _native.ConversionError
NonConcreteError = _native.NonConcreteError
PdbError = _native.PdbError
RenderError = _native.RenderError


class Format(str, Enum):
    """Input and output formats understood by crabWURCS."""

    AUTO = "auto"
    WURCS = "wurcs"
    IUPAC_CONDENSED = "iupac-condensed"
    IUPAC_EXTENDED = "iupac-extended"
    GLYCAM = "glycam"
    SMILES = "smiles"
    MOL = "mol"
    SDF = "sdf"


class ImageFormat(str, Enum):
    SVG = "svg"
    PNG = "png"


FormatLike = Union[Format, str]
ImageFormatLike = Union[ImageFormat, str]


def _format_value(value: FormatLike) -> str:
    return value.value if isinstance(value, Format) else str(value)


class Glycan:
    """An immutable glycan parsed into crabWURCS's shared graph model."""

    __slots__ = ("_value",)

    def __init__(self, value: object) -> None:
        if not isinstance(value, _native._Glycan):
            raise TypeError("Glycan objects must be created with Glycan.parse()")
        self._value = value

    @classmethod
    def parse(cls, value: str, format: FormatLike = Format.AUTO) -> "Glycan":
        """Parse WURCS, IUPAC, GLYCAM, SMILES, MOL, or SDF text."""

        return cls(_native._Glycan.parse(value, _format_value(format)))

    def to(self, format: FormatLike) -> str:
        """Serialize this glycan in the requested notation or molecular format."""

        return self._value.to_format(_format_value(format))

    def render(
        self,
        format: ImageFormatLike = ImageFormat.SVG,
        *,
        highlight_motifs: Optional[Iterable[str]] = None,
        motif_format: FormatLike = Format.AUTO,
    ) -> Union[str, bytes]:
        """Render an SNFG image, optionally highlighting motif occurrences."""

        image_format = format.value if isinstance(format, ImageFormat) else str(format).lower()
        motifs = list(highlight_motifs) if highlight_motifs is not None else None
        if image_format == ImageFormat.SVG.value:
            return self._value.render_svg(motifs, _format_value(motif_format))
        if image_format == ImageFormat.PNG.value:
            return bytes(self._value.render_png(motifs, _format_value(motif_format)))
        raise ValueError(f"unknown image format: {format}")

    def __repr__(self) -> str:
        try:
            return f"Glycan({self.to(Format.WURCS)!r})"
        except CrabWurcsError:
            return "Glycan(<not representable as WURCS>)"


@dataclass(frozen=True)
class PdbResidueReference:
    """Location of one glycan graph node in the source structure."""

    node_index: int
    chain: str
    sequence_number: int
    insertion_code: Optional[str]


@dataclass(frozen=True)
class ExtractedGlycan:
    """A glycan extracted from PDB/mmCIF with attachment and provenance."""

    glycan: Glycan
    attachment_site: Optional[str]
    residues: Sequence[PdbResidueReference]


def detect_format(value: str) -> Format:
    """Detect a text notation format. MOL/SDF should be supplied explicitly."""

    return Format(_native.detect_format(value))


def convert(
    value: str,
    to_format: FormatLike,
    from_format: FormatLike = Format.AUTO,
) -> str:
    """Convert a glycan between any supported text or molecular formats."""

    return Glycan.parse(value, from_format).to(to_format)


def render_snfg(
    value: str,
    *,
    from_format: FormatLike = Format.AUTO,
    image_format: ImageFormatLike = ImageFormat.SVG,
    highlight_motifs: Optional[Iterable[str]] = None,
    motif_format: FormatLike = Format.AUTO,
) -> Union[str, bytes]:
    """Parse and render a glycan as SNFG SVG text or PNG bytes."""

    return Glycan.parse(value, from_format).render(
        image_format,
        highlight_motifs=highlight_motifs,
        motif_format=motif_format,
    )


def _wrap_extracted(values: Iterable[object]) -> List[ExtractedGlycan]:
    results: List[ExtractedGlycan] = []
    for value in values:
        references = tuple(
            PdbResidueReference(
                node_index=reference.node_index,
                chain=reference.chain,
                sequence_number=reference.sequence_number,
                insertion_code=reference.insertion_code,
            )
            for reference in value.residues
        )
        results.append(
            ExtractedGlycan(
                glycan=Glycan(value.glycan),
                attachment_site=value.attachment_site,
                residues=references,
            )
        )
    return results


def extract_pdb(contents: str, format: str = "auto") -> List[ExtractedGlycan]:
    """Extract glycans and provenance from PDB or mmCIF text."""

    return _wrap_extracted(_native.extract_pdb(contents, format))


def extract_pdb_file(path: Union[str, PathLike[str]]) -> List[ExtractedGlycan]:
    """Extract glycans and provenance from a PDB or mmCIF file."""

    return _wrap_extracted(_native.extract_pdb_file(str(Path(path))))


__all__ = [
    "__version__",
    "Format",
    "ImageFormat",
    "Glycan",
    "PdbResidueReference",
    "ExtractedGlycan",
    "CrabWurcsError",
    "ParseError",
    "ConversionError",
    "NonConcreteError",
    "PdbError",
    "RenderError",
    "detect_format",
    "convert",
    "render_snfg",
    "extract_pdb",
    "extract_pdb_file",
]
