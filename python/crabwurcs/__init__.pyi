from enum import Enum
from os import PathLike
from typing import Iterable, List, Optional, Sequence, Union

__version__: str

class Format(str, Enum):
    AUTO: Format
    WURCS: Format
    IUPAC_CONDENSED: Format
    IUPAC_EXTENDED: Format
    GLYCAM: Format
    SMILES: Format
    MOL: Format
    SDF: Format

class ImageFormat(str, Enum):
    SVG: ImageFormat
    PNG: ImageFormat

FormatLike = Union[Format, str]
ImageFormatLike = Union[ImageFormat, str]

class CrabWurcsError(Exception): ...
class ParseError(CrabWurcsError): ...
class ConversionError(CrabWurcsError): ...
class NormalizationError(ConversionError): ...
class NonConcreteError(ConversionError): ...
class PdbError(CrabWurcsError): ...
class RenderError(CrabWurcsError): ...

class Glycan:
    @classmethod
    def parse(cls, value: str, format: FormatLike = ...) -> Glycan: ...
    def to(self, format: FormatLike) -> str: ...
    def render(self, format: ImageFormatLike = ..., *, highlight_motifs: Optional[Iterable[str]] = ..., motif_format: FormatLike = ...) -> Union[str, bytes]: ...

class PdbResidueReference:
    node_index: int
    chain: str
    sequence_number: int
    insertion_code: Optional[str]

class PdbWarning:
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

class StereoCentre:
    position: int
    atom: str
    code: str
    volume: float

class ResidueReading:
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

class ExtractedGlycan:
    glycan: Glycan
    attachment_site: Optional[str]
    residues: Sequence[PdbResidueReference]
    warnings: Sequence[PdbWarning]
    readings: Sequence[ResidueReading]

def detect_format(value: str) -> Format: ...
def convert(value: str, to_format: FormatLike, from_format: FormatLike = ..., *, table: str = ...) -> str: ...
def normalize_wurcs(value: str) -> str: ...
def render_snfg(value: str, *, from_format: FormatLike = ..., image_format: ImageFormatLike = ..., highlight_motifs: Optional[Iterable[str]] = ..., motif_format: FormatLike = ...) -> Union[str, bytes]: ...
def extract_pdb(contents: str, format: str = ..., *, policy: str = ...) -> List[ExtractedGlycan]: ...
def extract_pdb_file(path: Union[str, PathLike[str]], *, policy: str = ...) -> List[ExtractedGlycan]: ...
