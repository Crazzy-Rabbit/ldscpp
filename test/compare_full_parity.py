#!/usr/bin/env python3
"""Compare broader Python2 LDSC and C++ ldsc parity outputs."""

from __future__ import annotations

import argparse
import difflib
import gzip
import shlex
from pathlib import Path


def read_text(path: Path) -> str:
    if path.suffix == ".gz":
        with gzip.open(path, "rt", encoding="utf-8") as handle:
            return handle.read()
    return path.read_text(encoding="utf-8")


def normalize_log(text: str) -> str:
    lines = []
    skip_call = False
    skip_ldscore_diagnostics = False
    for line in text.splitlines():
        stripped = line.rstrip()
        if skip_call:
            if stripped == "":
                skip_call = False
            continue
        if skip_ldscore_diagnostics:
            if stripped.startswith("Analysis finished at "):
                skip_ldscore_diagnostics = False
                lines.append("Analysis finished at <TIME>")
            continue
        if stripped == "Call:":
            lines.append("Call: <PARAMETERS_CHECKED_SEPARATELY>")
            skip_call = True
        elif stripped.startswith("Summary of LD Scores in "):
            lines.append("Summary of LD Scores in <LDSCORE_OUTPUT>")
            lines.append("<LDSCORE_DIAGNOSTICS>")
            skip_ldscore_diagnostics = True
        elif line.startswith("* (C) "):
            lines.append("* (C) <AUTHOR>")
        elif line.startswith("* Broad Institute") or line.startswith("* West China Hospital"):
            lines.append("* <INSTITUTION>")
        elif line.startswith("Beginning analysis at "):
            lines.append("Beginning analysis at <TIME>")
        elif line.startswith("Analysis finished at "):
            lines.append("Analysis finished at <TIME>")
        elif line.startswith("Conversion finished at "):
            lines.append("Conversion finished at <TIME>")
        elif line.startswith("Total time elapsed: "):
            lines.append("Total time elapsed: <ELAPSED>")
        else:
            lines.append(stripped)
    return "\n".join(lines) + "\n"


def compare(label: str, py_text: str, cpp_text: str) -> bool:
    if py_text == cpp_text:
        print(f"OK {label}")
        return True
    print(f"DIFF {label}")
    for i, line in enumerate(
        difflib.unified_diff(
            py_text.splitlines(),
            cpp_text.splitlines(),
            fromfile=f"py/{label}",
            tofile=f"cpp/{label}",
            lineterm="",
        )
    ):
        if i >= 260:
            print("... diff truncated ...")
            break
        print(line)
    return False


def compare_numeric(label: str, py_text: str, cpp_text: str, tol: float = 1e-10) -> bool:
    try:
        py_vals = [float(x) for x in py_text.split()]
        cpp_vals = [float(x) for x in cpp_text.split()]
    except ValueError:
        return compare(label, py_text, cpp_text)
    if len(py_vals) != len(cpp_vals):
        return compare(label, py_text, cpp_text)
    max_diff = 0.0
    for a, b in zip(py_vals, cpp_vals):
        diff = abs(a - b)
        scale = max(1.0, abs(a), abs(b))
        max_diff = max(max_diff, diff)
        if diff > tol * scale:
            return compare(label, py_text, cpp_text)
    print(f"OK {label} (numeric max diff {max_diff:.3g})")
    return True


def compare_command_params(flow_dir: Path) -> bool:
    py_commands = (flow_dir / "commands_py.sh").read_text(encoding="utf-8").splitlines()
    cpp_commands = (flow_dir / "commands_cpp.sh").read_text(encoding="utf-8").splitlines()
    ok = True
    if len(py_commands) != len(cpp_commands):
        print("DIFF command parameter count")
        return False
    for idx, (py_cmd, cpp_cmd) in enumerate(zip(py_commands, cpp_commands), start=1):
        py_args = shlex.split(py_cmd)[1:]
        cpp_args = shlex.split(cpp_cmd)[1:]
        if py_args != cpp_args:
            ok = False
            print(f"DIFF command {idx} parameters")
            print(f"py : {' '.join(py_args)}")
            print(f"cpp: {' '.join(cpp_args)}")
    if ok:
        print("OK command parameters")
    return ok


def collect_outputs(root: Path) -> set[str]:
    out = root / "out"
    return {str(p.relative_to(root)).replace("\\", "/") for p in out.rglob("*") if p.is_file()}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-dir", required=True, type=Path)
    args = parser.parse_args()
    py = args.run_dir / "py"
    cpp = args.run_dir / "cpp"
    py_files = collect_outputs(py)
    cpp_files = collect_outputs(cpp)
    ok = True
    ok &= compare_command_params(args.run_dir.parent)
    if py_files != cpp_files:
        ok = False
        for rel in sorted(py_files - cpp_files):
            print(f"MISSING cpp/{rel}")
        for rel in sorted(cpp_files - py_files):
            print(f"EXTRA cpp/{rel}")

    for rel in sorted(py_files & cpp_files):
        py_text = read_text(py / rel)
        cpp_text = read_text(cpp / rel)
        if rel.endswith(".log"):
            py_text = normalize_log(py_text)
            cpp_text = normalize_log(cpp_text)
        if rel.endswith((".cov", ".delete", ".part_delete")):
            ok &= compare_numeric(rel, py_text, cpp_text)
        else:
            ok &= compare(rel, py_text, cpp_text)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
