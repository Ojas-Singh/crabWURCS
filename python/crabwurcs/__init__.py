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
NormalizationError = _native.NormalizationError
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
class PdbWarning:
    """Something a residue name and the coordinates disagree on, or that the
    coordinates cannot decide.

    ``kind`` is ``name-contradicts-coordinates``, ``missing-atoms``,
    ``distorted-geometry``, ``unreadable-geometry``, ``no-ring`` or
    ``unrecognized-group``. ``named`` and ``coordinates`` are short IUPAC
    names of what each says (``α-L-IdopA``); ``aspects`` says what differs
    (``C2``, ``anomer``, ``ring``, ``backbone``, ``C2 substituent``).
    """

    node_index: int
    chain: str
    sequence_number: int
    insertion_code: Optional[str]
    residue_name: str
    kind: str
    named: Optional[str]
    coordinates: Optional[str]
    positions: Sequence[int]
    aspects: Sequence[str]
    message: str

    def __str__(self) -> str:
        return self.message


@dataclass(frozen=True)
class StereoCentre:
    """One configuration read from a signed volume: the Fischer side ``1``
    (left) or ``2`` (right) and the volume normalised by the bond lengths."""

    position: int
    atom: str
    code: str
    volume: float


@dataclass(frozen=True)
class ResidueReading:
    """What the coordinates and the residue name say about one residue, as
    WURCS UniqueRES strings, with the ring, backbone and centres read."""

    node_index: int
    chain: str
    sequence_number: int
    insertion_code: Optional[str]
    residue_name: str
    from_coordinates: Optional[str]
    from_name: Optional[str]
    ring_atoms: Sequence[str]
    backbone_atoms: Sequence[str]
    centres: Sequence[StereoCentre]


@dataclass(frozen=True)
class ExtractedGlycan:
    """A glycan extracted from PDB/mmCIF with attachment and provenance.

    Every residue is read from its coordinates; ``warnings`` lists where a
    residue name disagrees with them, and ``readings`` what each residue's
    coordinates and name say.
    """

    glycan: Glycan
    attachment_site: Optional[str]
    residues: Sequence[PdbResidueReference]
    warnings: Sequence[PdbWarning] = ()
    readings: Sequence[ResidueReading] = ()


def detect_format(value: str) -> Format:
    """Detect a text notation format. MOL/SDF should be supplied explicitly."""

    return Format(_native.detect_format(value))


def convert(
    value: str,
    to_format: FormatLike,
    from_format: FormatLike = Format.AUTO,
    *,
    table: str = "off",
) -> str:
    """Convert a glycan between any supported text or molecular formats.

    Conversions are computed. ``table="fallback"`` answers notation inputs
    that cannot be converted from the bundled table of GlycoShape records;
    ``table="check"`` raises ``ConversionError`` when the input's record in
    that table holds another structure.
    """

    if table != "off":
        return _native.convert_with_table(
            value, _format_value(to_format), _format_value(from_format), table
        )
    return Glycan.parse(value, from_format).to(to_format)


def normalize_wurcs(value: str) -> str:
    """Return the normalized WURCS 2.0 representation of a WURCS string."""

    return _native.normalize_wurcs(value)


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
        warnings = tuple(
            PdbWarning(
                node_index=warning.node_index,
                chain=warning.chain,
                sequence_number=warning.sequence_number,
                insertion_code=warning.insertion_code,
                residue_name=warning.residue_name,
                kind=warning.kind,
                named=warning.named,
                coordinates=warning.coordinates,
                positions=tuple(warning.positions),
                aspects=tuple(warning.aspects),
                message=warning.message,
            )
            for warning in value.warnings
        )
        readings = tuple(
            ResidueReading(
                node_index=reading.node_index,
                chain=reading.chain,
                sequence_number=reading.sequence_number,
                insertion_code=reading.insertion_code,
                residue_name=reading.residue_name,
                from_coordinates=reading.from_coordinates,
                from_name=reading.from_name,
                ring_atoms=tuple(reading.ring_atoms),
                backbone_atoms=tuple(reading.backbone_atoms),
                centres=tuple(
                    StereoCentre(
                        position=centre.position,
                        atom=centre.atom,
                        code=centre.code,
                        volume=centre.volume,
                    )
                    for centre in reading.centres
                ),
            )
            for reading in value.readings
        )
        results.append(
            ExtractedGlycan(
                glycan=Glycan(value.glycan),
                attachment_site=value.attachment_site,
                residues=references,
                warnings=warnings,
                readings=readings,
            )
        )
    return results


def extract_pdb(
    contents: str, format: str = "auto", *, policy: str = "coordinates"
) -> List[ExtractedGlycan]:
    """Extract glycans and provenance from PDB or mmCIF text.

    Every monosaccharide is read from its coordinates. ``policy`` decides
    what happens when a residue name disagrees: ``coordinates`` (the default:
    the coordinates decide and the disagreement is a warning), ``names`` (the
    name decides where it says something; still a warning) or ``strict``
    (``PdbError``).
    """

    return _wrap_extracted(_native.extract_pdb(contents, format, policy))


def extract_pdb_file(
    path: Union[str, PathLike[str]], *, policy: str = "coordinates"
) -> List[ExtractedGlycan]:
    """Extract glycans and provenance from a PDB or mmCIF file (see
    ``extract_pdb`` for ``policy``)."""

    return _wrap_extracted(_native.extract_pdb_file(str(Path(path)), policy))


__all__ = [
    "__version__",
    "Format",
    "ImageFormat",
    "Glycan",
    "PdbResidueReference",
    "PdbWarning",
    "StereoCentre",
    "ResidueReading",
    "ExtractedGlycan",
    "CrabWurcsError",
    "ParseError",
    "ConversionError",
    "NormalizationError",
    "NonConcreteError",
    "PdbError",
    "RenderError",
    "detect_format",
    "convert",
    "normalize_wurcs",
    "render_snfg",
    "extract_pdb",
    "extract_pdb_file",
]
