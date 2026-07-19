from pathlib import Path

path = Path("README.md")
text = path.read_text(encoding="utf-8")

# Repair both a literal LaTeX command and the accidental newline produced when
# a backslash-n sequence was interpreted as a Python/JSON newline.
replacements = {
    "\\(\\nu=1\\)": "\\(ν=1\\)",
    "\\(\\nu=0\\)": "\\(ν=0\\)",
    "\\(\\nu\\)": "\\(ν\\)",
    "\\(\nu=1\\)": "\\(ν=1\\)",
    "\\(\nu=0\\)": "\\(ν=0\\)",
    "\\(\nu\\)": "\\(ν\\)",
    "|\\Delta\\nu|": "|\\Delta ν|",
    "|\\Delta\nu|": "|\\Delta ν|",
}

for old, new in replacements.items():
    text = text.replace(old, new)

path.write_text(text, encoding="utf-8")
