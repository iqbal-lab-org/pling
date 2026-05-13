import logging
import os
import re
from pathlib import Path
from typing import Optional, Union

CYTOSCAPE_CDN_URL = "https://unpkg.com/cytoscape@3.33.1/dist/cytoscape.min.js"

CYTOSCAPE_SOURCE_RE = re.compile(
    r"<script\b[^>]*\bsrc=[\"'][^\"']*cytoscape(?:\.min)?\.js[\"'][^>]*>\s*</script>",
    re.IGNORECASE,
)
LIBS_REFERENCE_RE = re.compile(
    r"<(?:script|link)\b[^>]*(?:src|href)=[\"'](?P<src>[^\"']*libs/[^\"']+)[\"']",
    re.IGNORECASE,
)


def repair_cytoscape_html_outputs(root: Union[Path, str]) -> int:
    """Ensure generated Cytoscape HTML files load the core Cytoscape.js library."""
    root = Path(root)
    if not root.exists():
        return 0

    repaired = 0
    for html_path in root.rglob("*.html"):
        if _repair_html_file(html_path):
            repaired += 1

    if repaired:
        logging.info("Added missing Cytoscape.js source tags to %d HTML file(s).", repaired)
    return repaired


def _repair_html_file(html_path: Path) -> bool:
    html = html_path.read_text(encoding="utf-8")
    if not _html_uses_cytoscape(html) or CYTOSCAPE_SOURCE_RE.search(html):
        return False

    cytoscape_src = _get_cytoscape_source(html_path, html)
    repaired_html = _insert_script_tag(html, cytoscape_src)
    if repaired_html == html:
        return False

    html_path.write_text(repaired_html, encoding="utf-8")
    return True


def _html_uses_cytoscape(html: str) -> bool:
    lowered = html.lower()
    return "cytoscape(" in lowered or "cytoscape." in lowered or "cytoscape-" in lowered


def _get_cytoscape_source(html_path: Path, html: str) -> str:
    libs_reference = LIBS_REFERENCE_RE.search(html)
    if libs_reference:
        src = libs_reference.group("src").replace("\\", "/")
        libs_prefix = src.split("libs/", 1)[0]
        return f"{libs_prefix}libs/cytoscape.min.js"

    local_cytoscape = _find_nearest_local_cytoscape(html_path)
    if local_cytoscape:
        return local_cytoscape

    return CYTOSCAPE_CDN_URL


def _find_nearest_local_cytoscape(html_path: Path) -> Optional[str]:
    for directory in (html_path.parent, *html_path.parents):
        cytoscape_path = directory / "libs" / "cytoscape.min.js"
        if cytoscape_path.exists():
            relative_path = os.path.relpath(cytoscape_path, html_path.parent)
            return relative_path.replace(os.sep, "/")
    return None


def _insert_script_tag(html: str, cytoscape_src: str) -> str:
    newline = "\r\n" if "\r\n" in html else "\n"
    lines = html.splitlines()

    for index, line in enumerate(lines):
        if "cytoscape-" in line.lower() and "<script" in line.lower():
            indent = line[: len(line) - len(line.lstrip())]
            script_tag = f'{indent}<script src="{cytoscape_src}"></script>'
            lines.insert(index, script_tag)
            return _join_lines(lines, newline, html)

    for index, line in enumerate(lines):
        if "<script" in line.lower():
            indent = line[: len(line) - len(line.lstrip())]
            script_tag = f'{indent}<script src="{cytoscape_src}"></script>'
            lines.insert(index, script_tag)
            return _join_lines(lines, newline, html)

    for index, line in enumerate(lines):
        if "</head>" in line.lower():
            indent = line[: len(line) - len(line.lstrip())]
            script_tag = f'{indent}<script src="{cytoscape_src}"></script>'
            lines.insert(index, script_tag)
            return _join_lines(lines, newline, html)

    script_tag = f'<script src="{cytoscape_src}"></script>'
    return f"{script_tag}{newline}{html}"


def _join_lines(lines: list[str], newline: str, original_html: str) -> str:
    repaired = newline.join(lines)
    if original_html.endswith(("\n", "\r\n")):
        repaired += newline
    return repaired
