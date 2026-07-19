from __future__ import annotations

import json
import re
from pathlib import Path

NOTEBOOK_PATH = Path("sshmodel.ipynb")



def source_text(cell: dict) -> str:
    source = cell.get("source", [])
    return "".join(source) if isinstance(source, list) else str(source)



def set_source(cell: dict, text: str) -> None:
    cell["source"] = text.splitlines(keepends=True)



def normalize_inline_math_spacing(text: str) -> str:
    """Add ASCII spaces around inline $...$ outside code and display math."""
    fenced_parts = re.split(r"(```.*?```)", text, flags=re.DOTALL)

    for fenced_index in range(0, len(fenced_parts), 2):
        display_parts = re.split(r"(\$\$.*?\$\$)", fenced_parts[fenced_index], flags=re.DOTALL)

        for display_index in range(0, len(display_parts), 2):
            plain = display_parts[display_index]

            # Insert one ASCII space when an inline-math delimiter directly
            # touches Japanese text, punctuation, or another non-space token.
            plain = re.sub(
                r"([^\s\n$])(\$[^\n$]+?\$)",
                r"\1 \2",
                plain,
            )
            plain = re.sub(
                r"(\$[^\n$]+?\$)([^\s\n$])",
                r"\1 \2",
                plain,
            )
            display_parts[display_index] = plain

        fenced_parts[fenced_index] = "".join(display_parts)

    return "".join(fenced_parts)



def inline_spacing_errors(text: str) -> list[str]:
    """Return snippets containing inline math without surrounding spaces."""
    errors: list[str] = []
    fenced_parts = re.split(r"(```.*?```)", text, flags=re.DOTALL)

    for fenced_index in range(0, len(fenced_parts), 2):
        display_parts = re.split(r"(\$\$.*?\$\$)", fenced_parts[fenced_index], flags=re.DOTALL)

        for display_index in range(0, len(display_parts), 2):
            plain = display_parts[display_index]
            patterns = (
                r"[^\s\n$]\$[^\n$]+?\$",
                r"\$[^\n$]+?\$[^\s\n$]",
            )
            for pattern in patterns:
                for match in re.finditer(pattern, plain):
                    start = max(0, match.start() - 20)
                    end = min(len(plain), match.end() + 20)
                    errors.append(plain[start:end].replace("\n", "\\n"))

    return errors



def repair_github_math(text: str) -> str:
    """Convert notebook Markdown to GitHub-supported math syntax."""
    # GitHub Markdown renders $...$ reliably, but not Jupyter-style \(...\).
    text = text.replace(r"\(", "$").replace(r"\)", "$")

    # GitHub rejects \operatorname in repository README math.
    text = re.sub(
        r"\\operatorname\{([^{}]+)\}",
        lambda match: rf"\mathrm{{{match.group(1)}}}",
        text,
    )

    # Negative thin space is exposed as a literal exclamation mark by GitHub's
    # renderer in this README, so omit it entirely.
    text = text.replace(r"\!", "")

    # Avoid cases/aligned for the two-line phase classification. GitHub has
    # reported false missing-end errors for both variants in this block.
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

    # Keep function names visually separated from their arguments.
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
        raise ValueError(f"Inline math without surrounding ASCII spaces:\n{preview}")


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
