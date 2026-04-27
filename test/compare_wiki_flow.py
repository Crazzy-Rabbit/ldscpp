#!/usr/bin/env python3
"""Compare Python2 LDSC and C++ ldsc outputs for the wiki-style flow."""

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
    for line in text.splitlines():
        stripped = line.rstrip()
        if skip_call:
            if stripped == "":
                skip_call = False
            continue
        if stripped == "Call:":
            lines.append("Call: <PARAMETERS_CHECKED_SEPARATELY>")
            skip_call = True
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
    for i, line in enumerate(difflib.unified_diff(
        py_text.splitlines(),
        cpp_text.splitlines(),
        fromfile=f"py/{label}",
        tofile=f"cpp/{label}",
        lineterm="",
    )):
        if i >= 220:
            print("... diff truncated ...")
            break
        print(line)
    return False


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


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-dir", required=True, type=Path)
    args = parser.parse_args()
    py = args.run_dir / "py"
    cpp = args.run_dir / "cpp"
    ok = True
    ok &= compare_command_params(args.run_dir.parent)

    exact = [
        "out/scz.sumstats.gz",
        "out/bip.sumstats.gz",
        "out/22.l2.ldscore.gz",
        "out/22.l2.M",
        "out/22.l2.M_5_50",
    ]
    for rel in exact:
        ok &= compare(rel, read_text(py / rel), read_text(cpp / rel))

    logs = [
        "out/scz.log",
        "out/bip.log",
        "out/rg.log",
        "out/h2.log",
        "out/22.log",
    ]
    for rel in logs:
        ok &= compare(rel, normalize_log(read_text(py / rel)), normalize_log(read_text(cpp / rel)))

    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
