"""Generate split f2py interface files for regular and reduced spharm builds."""

from __future__ import annotations

import re
import sys
from pathlib import Path

BLOCK_START = re.compile(r"^\s*subroutine\s+([a-zA-Z0-9_]+)\b")
BLOCK_END = re.compile(r"^\s*end\s+subroutine\b")
MODULE_END = re.compile(
    r"^(\s*end\s+python\s+module\s+)([a-zA-Z0-9_]+)(.*)$", re.IGNORECASE
)

REDUCED_SUBROUTINES = {
    "invlap",
    "lap",
    "multsmoothfact",
    "gaqd",
    "reduced_gaussian_legendre_basis",
    "reduced_gaussian_legendre_derivative_from_basis",
    "reduced_gaussian_legendre_derivative_basis",
    "reduced_gaussian_grdtospec",
    "reduced_gaussian_spectogrd",
    "reduced_gaussian_spectogrd_pair",
    "reduced_gaussian_getgrad",
    "reduced_gaussian_getgrad_pair",
    "reduced_gaussian_getvrtdivspec",
    "reduced_gaussian_getvrtspec",
    "reduced_gaussian_getdivspec",
}

DIRECT_REGULAR_SUBROUTINES = {
    "invlap",
    "lap",
    "multsmoothfact",
    "gaqd",
    "vhaes",
    "vhaesdiv",
    "vhaesvrt",
    "vhaesi",
    "vhags",
    "vhagsdiv",
    "vhagsvrt",
    "vhagsi",
    "vhses",
    "vhsesdiv",
    "vhsesvrt",
    "vhsesi",
    "vhsgs",
    "vhsgsdiv",
    "vhsgsvrt",
    "vhsgsi",
    "shaes",
    "shaesi",
    "shags",
    "shagsi",
    "shses",
    "shsesi",
    "shsgs",
    "shsgsi",
    "shaec",
    "shaeci",
    "shagc",
    "shagci",
    "shsec",
    "shseci",
    "shsgc",
    "shsgci",
    "vhaec",
    "vhaecdiv",
    "vhaecvrt",
    "vhaeci",
    "vhagc",
    "vhagcdiv",
    "vhagcvrt",
    "vhagci",
    "vhsec",
    "vhsecdiv",
    "vhsecvrt",
    "vhseci",
    "vhsgc",
    "vhsgcdiv",
    "vhsgcvrt",
    "vhsgci",
    "onedtotwod",
    "onedtotwod_vrtdiv",
    "onedtotwod_vrt",
    "onedtotwod_div",
    "twodtooned",
    "twodtooned_vrtdiv",
    "twodtooned_vrt",
    "twodtooned_div",
    "ihgeod",
    "getlegfunc",
    "specintrp",
}


def parse_blocks(lines: list[str]) -> tuple[list[str], dict[str, list[str]], list[str]]:
    header: list[str] = []
    blocks: dict[str, list[str]] = {}
    tail: list[str] = []
    i = 0

    while i < len(lines):
        match = BLOCK_START.match(lines[i])
        if match:
            break
        header.append(lines[i])
        i += 1

    while i < len(lines):
        if lines[i].lstrip().startswith("end interface"):
            tail = lines[i:]
            break
        match = BLOCK_START.match(lines[i])
        if not match:
            i += 1
            continue

        name = match.group(1)
        block: list[str] = [lines[i]]
        i += 1
        while i < len(lines):
            block.append(lines[i])
            if BLOCK_END.match(lines[i]):
                i += 1
                break
            i += 1
        blocks[name] = block

    return header, blocks, tail


def build_header(module_name: str) -> list[str]:
    return [
        "!    -*- f90 -*-\n",
        f"python module {module_name} ! in\n",
        f"    interface  ! in :{module_name}\n",
        "\n",
    ]


def main() -> int:
    if len(sys.argv) != 5:
        raise SystemExit(
            "usage: build_split_pyf.py <regular|reduced> <input.pyf> <output.pyf> <module_name>"
        )

    mode = sys.argv[1].strip().lower()
    input_path = Path(sys.argv[2])
    output_path = Path(sys.argv[3])
    module_name = sys.argv[4].strip()

    if mode not in {"regular", "reduced"}:
        raise SystemExit(f"unsupported mode: {mode}")

    header, blocks, tail = parse_blocks(
        input_path.read_text(encoding="utf-8").splitlines(True)
    )

    if mode == "regular":
        selected_names = [
            name
            for name in blocks
            if name not in REDUCED_SUBROUTINES
            and name not in DIRECT_REGULAR_SUBROUTINES
        ]
    else:
        selected_names = [name for name in blocks if name in REDUCED_SUBROUTINES]

    output_lines = build_header(module_name)
    for name in selected_names:
        output_lines.extend(blocks[name])
    if tail:
        output_lines.append("\n")
        for line in tail:
            match = MODULE_END.match(line)
            if match:
                output_lines.append(f"{match.group(1)}{module_name}{match.group(3)}\n")
            else:
                output_lines.append(line)
    else:
        output_lines.extend(
            ["    end interface\n", f"end python module {module_name}\n"]
        )

    output_path.write_text("".join(output_lines), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
