#!/usr/bin/env python3
# ==============================================================================
# Project     : Chromatin-state context of SREBF2 knockout-responsive genes in HeLa cells: an exploratory integrative epigenomic analysis
# File        : add_headers.py
# Purpose     : Preview and add consistent attribution headers to tracked code.
# Author      : Md Tariqul Islam
# Institution : Northeastern University | Bioinformatics Department
# License     : MIT (see LICENSE)
# ==============================================================================
"""Add attribution headers to repository source files without changing code logic.

Default behavior is a read-only preview. Use --apply to write changes.
Only tracked files are selected unless --include-untracked is specified.
"""

from __future__ import annotations

import argparse
import ast
import re
import subprocess
import sys
from pathlib import Path

PROJECT = (
    "Chromatin-state context of SREBF2 knockout-responsive genes in HeLa "
    "cells: an exploratory integrative epigenomic analysis"
)
AUTHOR = "Md Tariqul Islam"
INSTITUTION = "Northeastern University | Bioinformatics Department"
EXCLUDED_DIRS = {
    ".git", ".venv", ".venv-cuttag", ".r_libs_phase15",
    ".r_libs_phase15_docker", "env", "venv", "node_modules",
    "__pycache__", "site-packages", "results", "metadata", "data",
    "raw", "reports", "results_public",
}
HEADER_MARKER = "# Project     : Chromatin-state context of SREBF2"
CODING_RE = re.compile(r"^#.*coding[:=]\s*[-\w.]+")
PY_EXTENSIONS = {".py"}
R_EXTENSIONS = {".r"}
SNAKEMAKE_EXTENSIONS = {".smk", ".snakefile"}
MAX_PURPOSE = 72


def source_files(root: Path, include_untracked: bool) -> list[Path]:
    command = ["git", "-C", str(root), "ls-files", "-z", "--cached"]
    if include_untracked:
        command.extend(["--others", "--exclude-standard"])
    output = subprocess.check_output(command)
    paths = []
    for raw_path in output.split(b"\0"):
        if not raw_path:
            continue
        relative = Path(raw_path.decode("utf-8", "surrogateescape"))
        if any(part in EXCLUDED_DIRS or part.startswith(".venv") for part in relative.parts[:-1]):
            continue
        name = relative.name.lower()
        if name != "snakefile" and relative.suffix.lower() not in (
            PY_EXTENSIONS | R_EXTENSIONS | SNAKEMAKE_EXTENSIONS
        ):
            continue
        path = root / relative
        if path.is_file() and not path.is_symlink():
            paths.append(relative)
    return sorted(set(paths), key=str)


def first_sentence(value: str) -> str:
    value = re.sub(r"\s+", " ", value.strip().strip("#* \t"))
    return re.split(r"(?<=[.!?])\s+", value, maxsplit=1)[0].rstrip(".")


def purpose_for(path: Path, source: str) -> str:
    # Prefer the author's actual module documentation to guessed descriptions.
    if path.suffix.lower() == ".py":
        try:
            doc = ast.get_docstring(ast.parse(source))
        except SyntaxError:
            doc = None
        if doc:
            return first_sentence(doc)
    else:
        lines = source.splitlines()
        for line in lines[:18]:
            stripped = line.strip()
            if stripped.startswith("#") and not stripped.startswith("#!"):
                candidate = first_sentence(stripped.lstrip("# "))
                if candidate and not candidate.lower().startswith(
                    ("project", "file", "author", "license", "institution", "====")
                ):
                    return candidate
    return "Run " + path.stem.replace("_", " ").replace("-", " ") + " analysis"


def wrap_purpose(text: str) -> list[str]:
    import textwrap

    chunks = textwrap.wrap(text, width=MAX_PURPOSE, break_long_words=False)
    return chunks[:5] or ["Project analysis script"]


def header(relative: Path, purpose: str, newline: str) -> str:
    lines = [
        "# " + "=" * 78,
        f"# Project     : {PROJECT}",
        f"# File        : {relative.as_posix()}",
    ]
    purpose_lines = wrap_purpose(purpose)
    lines.append(f"# Purpose     : {purpose_lines[0]}")
    lines.extend(f"#               {line}" for line in purpose_lines[1:])
    lines.extend([
        f"# Author      : {AUTHOR}",
        f"# Institution : {INSTITUTION}",
        "# License     : MIT (see LICENSE)",
        "# " + "=" * 78,
        "",
    ])
    return newline.join(lines) + newline


def insert_header(data: bytes, relative: Path) -> tuple[bytes, str]:
    if b"\0" in data:
        return data, "binary-skip"
    try:
        decoded = data.decode("utf-8-sig")
    except UnicodeDecodeError:
        return data, "non-utf8-skip"
    if HEADER_MARKER in decoded[:4000]:
        return data, "already-present"
    newline = "\r\n" if b"\r\n" in data else "\n"
    # Retain shebang and Python encoding declaration within the first two lines.
    lines = decoded.splitlines(keepends=True)
    insertion = 0
    if lines and lines[0].startswith("#!"):
        insertion = 1
    if insertion < min(2, len(lines)) and CODING_RE.match(lines[insertion].rstrip("\r\n")):
        insertion += 1
    prefix = "".join(lines[:insertion])
    suffix = "".join(lines[insertion:])
    if prefix and not prefix.endswith(("\n", "\r")):
        prefix += newline
    modified = prefix + header(relative, purpose_for(relative, decoded), newline) + suffix
    return (b"\xef\xbb\xbf" if data.startswith(b"\xef\xbb\xbf") else b"") + modified.encode("utf-8"), "new"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--apply", action="store_true", help="Write changes (default: preview only)")
    parser.add_argument("--include-untracked", action="store_true", help="Include untracked, nonignored source files")
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parent, help="Repository root")
    args = parser.parse_args()
    root = args.root.resolve()
    if not (root / ".git").exists():
        parser.error(f"Not a Git working-tree root: {root}")
    if args.apply and not (root / "LICENSE").is_file():
        parser.error("LICENSE missing. Add and verify MIT LICENSE before applying MIT headers.")
    paths = source_files(root, args.include_untracked)
    counts = {"new": 0, "already-present": 0, "binary-skip": 0, "non-utf8-skip": 0}
    for relative in paths:
        path = root / relative
        data = path.read_bytes()
        updated, status = insert_header(data, relative)
        counts[status] += 1
        if status == "new":
            print(f"{'UPDATE' if args.apply else 'WOULD UPDATE'}  {relative.as_posix()}")
            if args.apply:
                path.write_bytes(updated)
        else:
            print(f"SKIP ({status})  {relative.as_posix()}")
    print(f"\nFiles inspected: {len(paths)}")
    print("Counts:", counts)
    print("Mode:", "APPLIED" if args.apply else "DRY RUN — no changes made")
    return 0


if __name__ == "__main__":
    sys.exit(main())
