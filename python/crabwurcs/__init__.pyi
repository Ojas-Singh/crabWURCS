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

class ExtractedGlycan:
    glycan: Glycan
    attachment_site: Optional[str]
    residues: Sequence[PdbResidueReference]

def detect_format(value: str) -> Format: ...
def convert(value: str, to_format: FormatLike, from_format: FormatLike = ...) -> str: ...
def render_snfg(value: str, *, from_format: FormatLike = ..., image_format: ImageFormatLike = ..., highlight_motifs: Optional[Iterable[str]] = ..., motif_format: FormatLike = ...) -> Union[str, bytes]: ...
def extract_pdb(contents: str, format: str = ...) -> List[ExtractedGlycan]: ...
def extract_pdb_file(path: Union[str, PathLike[str]]) -> List[ExtractedGlycan]: ...
