"""Glycan figures: grids of glycans, and glycans in place of text labels.

``annotate_svg`` replaces the text labels of any SVG that name a glycan with
its SNFG drawing, and ``annotate_figure`` does so for a matplotlib figure, for
heatmaps, volcano plots and bar charts with glycans on their axes. Nothing here
needs matplotlib until ``annotate_figure`` is called.
"""

from __future__ import annotations

import html
import io
import re
from contextlib import contextmanager
from dataclasses import dataclass
from os import PathLike
from pathlib import Path
from typing import Any, Dict, Iterator, List, Mapping, Optional, Sequence, Tuple, Union

from . import CrabWurcsError, Format, FormatLike, Glycan, ImageFormat, ImageFormatLike, _native

_TEXT = re.compile(r"<text\b([^>]*)>(.*?)</text>", re.DOTALL)
_ATTRIBUTE = re.compile(r"""([\w:-]+)\s*=\s*(?:"([^"]*)"|'([^']*)')""")
_TAG = re.compile(r"<[^>]+>")
_FONT_SIZE = re.compile(r"font(?:-size)?:[^;]*?([\d.]+)px")
_ANCHOR = re.compile(r"text-anchor:\s*(start|middle|end)")
_SVG_START = re.compile(r"<svg\b[^>]*>")
_DOCUMENT_ONLY = re.compile(
    r"<title\b.*?</title>\s*|<desc\b.*?</desc>\s*|<metadata\b.*?</metadata>\s*", re.DOTALL
)
_VIEWBOX = re.compile(r'viewBox="([^"]*)"')

#: How many times the font size the tallest drawing in a document is.
LABEL_HEIGHT = 3.0
#: ``svg_to_png`` draws 2 pixels per CSS pixel, of which an inch has 96.
_PNG_DPI = 192


def _drawing(options: Mapping[str, Any]) -> Dict[str, Any]:
    known = {
        "orientation",
        "compact",
        "reducing_end",
        "repeat_brackets",
        "floating_fragments",
        "show_labels",
        "show_linkages",
        "colour",
        "scale",
    }
    unknown = set(options) - known
    if unknown:
        raise TypeError(f"unknown drawing options: {', '.join(sorted(unknown))}")
    return {key: value for key, value in options.items() if value is not None}


def render_grid(
    glycans: Sequence[Union[Glycan, str]],
    *,
    columns: int = 4,
    labels: Optional[Sequence[str]] = None,
    image_format: ImageFormatLike = ImageFormat.SVG,
    from_format: FormatLike = Format.AUTO,
    **drawing: Any,
) -> Union[str, bytes]:
    """Draw several glycans on one canvas, as SVG text or PNG bytes.

    ``labels`` are captions, one per glycan. ``drawing`` takes the drawing
    options of ``Glycan.render``, such as ``compact`` and ``orientation``.
    """

    parsed = [
        glycan if isinstance(glycan, Glycan) else Glycan.parse(glycan, from_format)
        for glycan in glycans
    ]
    values = [glycan._value for glycan in parsed]
    captions = list(labels) if labels is not None else None
    options = _drawing(drawing)
    image = image_format.value if isinstance(image_format, ImageFormat) else str(image_format).lower()
    if image == ImageFormat.SVG.value:
        return _native.render_grid_svg(values, columns, captions, options)
    if image == ImageFormat.PNG.value:
        return bytes(_native.render_grid_png(values, columns, captions, options))
    raise ValueError(f"unknown image format: {image_format}")


def write_images(
    glycans: Mapping[str, Union[Glycan, str]],
    directory: Union[str, PathLike],
    *,
    image_format: ImageFormatLike = ImageFormat.SVG,
    from_format: FormatLike = Format.AUTO,
    **drawing: Any,
) -> Sequence[Path]:
    """Write one image per glycan, named by the keys of ``glycans``."""

    image = image_format.value if isinstance(image_format, ImageFormat) else str(image_format).lower()
    target = Path(directory)
    target.mkdir(parents=True, exist_ok=True)
    written = []
    for name, glycan in glycans.items():
        parsed = glycan if isinstance(glycan, Glycan) else Glycan.parse(glycan, from_format)
        output = parsed.render(image, **drawing)
        path = target / f"{name}.{image}"
        if isinstance(output, bytes):
            path.write_bytes(output)
        else:
            path.write_text(output, encoding="utf-8")
        written.append(path)
    return written


def svg_to_png(svg: str) -> bytes:
    """Rasterize any SVG as a transparent PNG at twice its dimensions."""

    return bytes(_native.svg_to_png(svg))


def _glycan_in(text: str, labels: Optional[Mapping[str, Any]]) -> Optional[Glycan]:
    """The glycan a text label names. Without ``labels`` a text counts only if
    every residue in it is recognized: the parser reads any word as a residue
    of unknown chemistry, and an axis title is not a glycan."""

    if labels is not None:
        notation = labels.get(text)
        if notation is None:
            return None
        return notation if isinstance(notation, Glycan) else Glycan.parse(notation)
    if not any(character.isalpha() for character in text):
        return None
    try:
        glycan = Glycan.parse(text)
        residues = glycan.residues()
    except CrabWurcsError:
        return None
    return glycan if residues and all(residue.recognized for residue in residues) else None


def _size_of(svg: str) -> Tuple[float, float]:
    box = _VIEWBOX.search(svg)
    if box is None:
        raise ValueError("a generated drawing has no viewBox")
    _, _, width, height = box.group(1).split()
    return float(width), float(height)


@dataclass
class _Label:
    """A text element of a document that names a glycan."""

    start: int
    end: int
    attributes: Dict[str, str]
    font_size: float
    anchor: str
    body: str
    width: float
    height: float

    def nested(self, factor: float) -> str:
        """The drawing as a nested ``<svg>`` where the text was, ``factor``
        times its natural size."""

        width = self.width * factor
        height = self.height * factor
        x = float(self.attributes.get("x", "0"))
        y = float(self.attributes.get("y", "0"))
        left = {"start": x, "middle": x - width / 2.0, "end": x - width}[self.anchor]
        # A text's y is its baseline, about a third of the font below its middle.
        top = y - self.font_size * 0.35 - height / 2.0
        box = _VIEWBOX.search(self.body).group(1)  # type: ignore[union-attr]
        start = (
            f'<svg x="{left}" y="{top}" width="{width}" height="{height}" viewBox="{box}"'
            ' xmlns="http://www.w3.org/2000/svg">'
        )
        nested = _SVG_START.sub(start, self.body, count=1)
        transform = self.attributes.get("transform")
        return f'<g transform="{transform}">{nested}</g>' if transform else nested


def _find_labels(
    svg: str, labels: Optional[Mapping[str, Any]], drawing: Mapping[str, Any]
) -> List[_Label]:
    found = []
    for match in _TEXT.finditer(svg):
        text = html.unescape(_TAG.sub("", match.group(2))).strip()
        glycan = _glycan_in(text, labels) if text else None
        if glycan is None:
            continue
        attributes = {
            name: double if double else single
            for name, double, single in _ATTRIBUTE.findall(match.group(1))
        }
        style = attributes.get("style", "")
        font = _FONT_SIZE.search(style)
        size = font.group(1) if font else attributes.get("font-size", "10").rstrip("px")
        anchor = _ANCHOR.search(style)
        rendered = glycan.render(**drawing)
        assert isinstance(rendered, str)
        width, height = _size_of(rendered)
        found.append(
            _Label(
                start=match.start(),
                end=match.end(),
                attributes=attributes,
                font_size=float(size),
                anchor=anchor.group(1) if anchor else attributes.get("text-anchor", "start"),
                body=_DOCUMENT_ONLY.sub("", rendered),
                width=width,
                height=height,
            )
        )
    return found


def _scale_of(found: Sequence[_Label], height: Optional[float]) -> float:
    """One scale for every drawing, so a symbol is the same size in all of
    them: the tallest drawing gets the requested height."""

    tallest = max(label.height for label in found)
    target = height or max(label.font_size for label in found) * LABEL_HEIGHT
    return target / tallest


def annotate_svg(
    svg: str,
    *,
    labels: Optional[Mapping[str, Any]] = None,
    height: Optional[float] = None,
    **drawing: Any,
) -> str:
    """Replace the text labels of an SVG that name a glycan with its drawing.

    The rest of the document stays as it is, and a drawing takes the position,
    anchor and rotation of the text it replaces. ``labels`` maps label texts to
    the glycans they stand for, for short names such as ``G0F``; without it,
    any label that reads as a fully recognized glycan is replaced. All drawings
    share one scale, so their symbols match: ``height`` is the height of the
    tallest one in the document's units, by default three times the font size.
    ``drawing`` takes the options of ``Glycan.render``, compact by default.
    """

    options = dict(drawing)
    options.setdefault("compact", True)
    found = _find_labels(svg, labels, options)
    if not found:
        return svg
    factor = _scale_of(found, height)
    parts = []
    position = 0
    for label in found:
        parts.append(svg[position : label.start])
        parts.append(label.nested(factor))
        position = label.end
    parts.append(svg[position:])
    return "".join(parts)


def _glycan_texts(figure: Any, labels: Optional[Mapping[str, Any]]) -> List[Any]:
    """The visible matplotlib text artists that name a glycan."""

    from matplotlib.text import Text

    return [
        text
        for text in figure.findobj(Text)
        if text.get_visible()
        and text.get_text().strip()
        and _glycan_in(text.get_text().strip(), labels) is not None
    ]


@contextmanager
def _rotated_labels_anchored(figure: Any, labels: Optional[Mapping[str, Any]]) -> Iterator[None]:
    """Anchor rotated glycan labels at the tick they belong to.

    A rotated text is placed by its own measured box, which says nothing about
    a drawing of another length. Anchored at its end, the reducing end of the
    drawing meets the tick whatever its length.
    """

    rotated = [text for text in _glycan_texts(figure, labels) if text.get_rotation() % 360]
    saved = [(text.get_rotation_mode(), text.get_ha(), text.get_va()) for text in rotated]
    for text in rotated:
        text.set_rotation_mode("anchor")
        text.set_ha("right")
        text.set_va("center")
    try:
        yield
    finally:
        for text, (mode, horizontal, vertical) in zip(rotated, saved):
            text.set_rotation_mode(mode)
            text.set_ha(horizontal)
            text.set_va(vertical)


def _figure_png(
    figure: Any,
    svg: str,
    labels: Optional[Mapping[str, Any]],
    height: Optional[float],
    drawing: Mapping[str, Any],
) -> bytes:
    """The figure rasterized by matplotlib with the drawings laid over it.

    ``svg_to_png`` cannot be used on the whole document: it has no raster
    decoders, so it would drop what ``imshow`` embeds.
    """

    from PIL import Image

    options = dict(drawing)
    options.setdefault("compact", True)
    found = _find_labels(svg, labels, options)
    hidden = _glycan_texts(figure, labels)
    for text in hidden:
        text.set_visible(False)
    try:
        buffer = io.BytesIO()
        figure.savefig(buffer, format="png", dpi=_PNG_DPI)
    finally:
        for text in hidden:
            text.set_visible(True)
    base = Image.open(buffer).convert("RGBA")
    if found:
        factor = _scale_of(found, height)
        header = _SVG_START.search(svg).group(0)  # type: ignore[union-attr]
        overlay = header + "".join(label.nested(factor) for label in found) + "</svg>\n"
        drawn = Image.open(io.BytesIO(svg_to_png(overlay))).convert("RGBA")
        base.alpha_composite(drawn.crop((0, 0, base.width, base.height)))
    output = io.BytesIO()
    base.save(output, format="PNG")
    return output.getvalue()


def annotate_figure(
    figure: Any,
    path: Optional[Union[str, PathLike]] = None,
    *,
    labels: Optional[Mapping[str, Any]] = None,
    height: Optional[float] = None,
    **drawing: Any,
) -> str:
    """Save a matplotlib figure with glycan labels drawn as SNFG symbols.

    Returns the SVG text, and writes it to ``path`` if given: a ``.png`` path
    is rasterized, anything else is written as SVG. See ``annotate_svg`` for
    ``labels``, ``height`` and ``drawing``. Leave room around the axes (with
    ``tick_params(pad=...)`` or ``subplots_adjust``): matplotlib lays the
    figure out for the text the drawings replace.
    """

    import matplotlib

    buffer = io.StringIO()
    # Text has to stay text, or there is nothing to find; and both outputs
    # have to share one layout, which a tight bounding box would change.
    settings = {"svg.fonttype": "none", "savefig.bbox": "standard"}
    with matplotlib.rc_context(settings), _rotated_labels_anchored(figure, labels):
        figure.savefig(buffer, format="svg")
        source = buffer.getvalue()
        annotated = annotate_svg(source, labels=labels, height=height, **drawing)
        if path is not None:
            target = Path(path)
            if target.suffix.lower() == ".png":
                target.write_bytes(_figure_png(figure, source, labels, height, drawing))
            else:
                target.write_text(annotated, encoding="utf-8")
    return annotated


__all__ = [
    "LABEL_HEIGHT",
    "annotate_figure",
    "annotate_svg",
    "render_grid",
    "svg_to_png",
    "write_images",
]
