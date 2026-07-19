from __future__ import annotations

import json
import re
from pathlib import Path

NOTEBOOK_PATH = Path("sshmodel.ipynb")
INLINE_MATH_RE = re.compile(r"(?<!\$)\$[^\n$]+\$(?!\$)")


def source_text(cell: dict) -> str:
    source = cell.get("source", [])
    return "".join(source) if isinstance(source, list) else str(source)


def set_source(cell: dict, text: str) -> None:
    cell["source"] = text.splitlines(keepends=True)


def normalize_plain_text_inline_math(text: str) -> str:
    """Add spaces outside inline math without changing its contents."""
    output: list[str] = []
    last = 0

    for match in INLINE_MATH_RE.finditer(text):
        prefix = text[last:match.start()]
        if match.start() > 0 and not text[match.start() - 1].isspace():
            prefix += " "
        output.append(prefix)
        output.append(match.group(0))
        if match.end() < len(text) and not text[match.end()].isspace():
            output.append(" ")
        last = match.end()

    output.append(text[last:])
    return "".join(output)


def normalize_inline_math_spacing(text: str) -> str:
    """Normalize inline $...$ outside fenced code and display math."""
    fenced_parts = re.split(r"(```.*?```)", text, flags=re.DOTALL)

    for fenced_index in range(0, len(fenced_parts), 2):
        display_parts = re.split(
            r"(\$\$.*?\$\$)",
            fenced_parts[fenced_index],
            flags=re.DOTALL,
        )
        for display_index in range(0, len(display_parts), 2):
            display_parts[display_index] = normalize_plain_text_inline_math(
                display_parts[display_index]
            )
        fenced_parts[fenced_index] = "".join(display_parts)

    return "".join(fenced_parts)


def inline_spacing_errors(text: str) -> list[str]:
    """Return malformed inline-math snippets outside code/display math."""
    errors: list[str] = []
    fenced_parts = re.split(r"(```.*?```)", text, flags=re.DOTALL)

    for fenced_index in range(0, len(fenced_parts), 2):
        display_parts = re.split(
            r"(\$\$.*?\$\$)",
            fenced_parts[fenced_index],
            flags=re.DOTALL,
        )
        for display_index in range(0, len(display_parts), 2):
            plain = display_parts[display_index]
            for match in INLINE_MATH_RE.finditer(plain):
                content = match.group(0)[1:-1]
                malformed = (
                    (match.start() > 0 and not plain[match.start() - 1].isspace())
                    or (match.end() < len(plain) and not plain[match.end()].isspace())
                    or content.startswith(" ")
                    or content.endswith(" ")
                )
                if malformed:
                    start = max(0, match.start() - 20)
                    end = min(len(plain), match.end() + 20)
                    errors.append(plain[start:end].replace("\n", "\\n"))

    return errors


def repair_github_math(text: str) -> str:
    """Convert notebook Markdown to GitHub-supported math syntax."""
    text = text.replace(r"\(", "$").replace(r"\)", "$")
    text = re.sub(
        r"\\operatorname\{([^{}]+)\}",
        lambda match: rf"\mathrm{{{match.group(1)}}}",
        text,
    )
    text = text.replace(r"\!", "")

    phase_classification = (
        r"\nu=0 \quad (v>w,\ \mathrm{trivial}),"
        "\n"
        r"\qquad"
        "\n"
        r"\nu=1 \quad (v<w,\ \mathrm{topological})."
    )
    text = re.sub(
        r"\\nu\s*=\s*\\begin\{cases\}.*?\\end\{cases\}",
        lambda _match: phase_classification,
        text,
        flags=re.DOTALL,
    )
    text = re.sub(
        r"\\begin\{aligned\}\s*\\nu.*?\\mathrm\{trivial\}.*?"
        r"\\mathrm\{topological\}.*?\\end\{aligned\}",
        lambda _match: phase_classification,
        text,
        flags=re.DOTALL,
    )
    text = re.sub(
        r"(\\mathrm\{(?:Arg|atan2|diag)\})(?=[A-Za-z\\])",
        r"\1\,",
        text,
    )
    return normalize_inline_math_spacing(text)


def validate(markdown: str) -> None:
    forbidden = {
        r"\(": "Jupyter inline-math opener remains",
        r"\)": "Jupyter inline-math closer remains",
        r"\operatorname": "unsupported operatorname remains",
        r"\begin{cases}": "unsupported cases environment remains",
        r"\end{cases}": "unsupported cases environment remains",
        r"\!": "negative thin-space command remains",
        r"\nu&=0": "fragile aligned winding block remains",
    }
    for token, message in forbidden.items():
        if token in markdown:
            raise ValueError(f"{message}: {token}")

    if markdown.count("$$") % 2:
        raise ValueError("Unbalanced display-math delimiters ($$).")

    spacing_errors = inline_spacing_errors(markdown)
    if spacing_errors:
        preview = "\n".join(spacing_errors[:10])
        raise ValueError(
            "Inline math must have external spaces and no internal edge spaces:\n"
            + preview
        )


if not NOTEBOOK_PATH.exists():
    raise FileNotFoundError(f"{NOTEBOOK_PATH} was not found.")

notebook = json.loads(NOTEBOOK_PATH.read_text(encoding="utf-8"))

for cell in notebook.get("cells", []):
    if cell.get("cell_type") != "markdown":
        continue
    set_source(cell, repair_github_math(source_text(cell)))

markdown = "\n".join(
    source_text(cell)
    for cell in notebook.get("cells", [])
    if cell.get("cell_type") == "markdown"
)
validate(markdown)

NOTEBOOK_PATH.write_text(
    json.dumps(notebook, ensure_ascii=False, indent=1) + "\n",
    encoding="utf-8",
)
