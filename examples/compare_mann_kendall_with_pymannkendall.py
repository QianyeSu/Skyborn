"""
Compare the current Skyborn Mann-Kendall implementation with pymannkendall.

This script intentionally loads the local development copy of
``src/skyborn/calc/mann_kendall.py`` so the comparison reflects the code in the
current workspace rather than any previously installed Skyborn package.
"""

from __future__ import annotations

import argparse
import importlib.util
import pathlib
from typing import Iterable

import numpy as np
import pymannkendall as pmk

ROOT = pathlib.Path(__file__).resolve().parents[1]
MANN_KENDALL_PATH = ROOT / "src" / "skyborn" / "calc" / "mann_kendall.py"


def load_local_mann_kendall():
    """Load the local development mann_kendall module directly from source."""
    spec = importlib.util.spec_from_file_location(
        "skyborn_dev_mann_kendall", MANN_KENDALL_PATH
    )
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Unable to load module from {MANN_KENDALL_PATH}")

    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def equal_or_both_nan(left, right, atol: float = 1e-12) -> bool:
    """Treat NaN/NaN as equal while keeping exact-value comparisons strict."""
    left_is_nan = np.isscalar(left) and np.isnan(left)
    right_is_nan = np.isscalar(right) and np.isnan(right)
    if left_is_nan and right_is_nan:
        return True
    return bool(np.allclose(left, right, rtol=0.0, atol=atol))


def build_cases(include_strict_linear: bool) -> dict[str, np.ndarray]:
    """Return deterministic sample cases for comparison."""
    cases = {
        "ties_up": np.array([1, 1, 2, 2, 3, 3, 4, 4, 5, 5], dtype=float),
        "ties_mixed": np.array([1, 2, 2, 3, 3, 3, 4, 4, 5, 5], dtype=float),
        "oscillating": np.array([1, 3, 2, 4, 3, 5, 4, 6, 5, 7], dtype=float),
        "random_seeded": np.array(
            [
                0.12573022,
                -0.13210486,
                0.64042265,
                0.10490012,
                -0.53566937,
                0.36159505,
                1.30400005,
                0.94708096,
                -0.70373524,
                -1.26542147,
            ],
            dtype=float,
        ),
    }
    if include_strict_linear:
        cases = {
            "strict_increase": np.arange(1, 11, dtype=float),
            **cases,
        }
    return cases


def iter_modified_flags(mode: str) -> Iterable[bool]:
    """Map CLI mode to modified flags."""
    if mode == "standard":
        return (False,)
    if mode == "modified":
        return (True,)
    return (False, True)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Compare the current local Skyborn Mann-Kendall implementation "
            "with pymannkendall on a fixed set of deterministic sample cases."
        )
    )
    parser.add_argument(
        "--mode",
        choices=("standard", "modified", "both"),
        default="both",
        help="Which Mann-Kendall mode to compare.",
    )
    parser.add_argument(
        "--method",
        choices=("theilslopes", "linregress"),
        default="theilslopes",
        help=(
            "Slope method passed to Skyborn. "
            "Direct slope parity with pymannkendall is meaningful for 'theilslopes'."
        ),
    )
    parser.add_argument(
        "--include-strict-linear",
        action="store_true",
        help=(
            "Also include a perfectly linear increasing series. "
            "In modified mode this case can legitimately yield NaN z/p in both libraries."
        ),
    )
    args = parser.parse_args()

    mk = load_local_mann_kendall()
    cases = build_cases(include_strict_linear=args.include_strict_linear)

    print(f"Loaded local Skyborn development module from: {MANN_KENDALL_PATH}")
    print(f"pymannkendall version: {getattr(pmk, '__version__', 'unknown')}")
    print(f"Skyborn slope method: {args.method}")

    if args.method != "theilslopes":
        print(
            "Note: pymannkendall always reports Sen/Theil-Sen slope. "
            "With method='linregress', slope equality is not expected."
        )

    for modified in iter_modified_flags(args.mode):
        print("\n" + "=" * 88)
        print(f"modified = {modified}")
        print("=" * 88)

        for name, arr in cases.items():
            sky = mk.mann_kendall_test(arr, method=args.method, modified=modified)
            pmk_result = (
                pmk.yue_wang_modification_test(arr)
                if modified
                else pmk.original_test(arr)
            )

            trend_equal = equal_or_both_nan(sky["trend"], pmk_result.slope, atol=0.0)
            h_equal = sky["h"] == pmk_result.h
            p_equal = equal_or_both_nan(sky["p"], pmk_result.p)
            z_equal = equal_or_both_nan(sky["z"], pmk_result.z)
            tau_equal = equal_or_both_nan(sky["tau"], pmk_result.Tau)
            all_equal = trend_equal and h_equal and p_equal and z_equal and tau_equal

            print(f"\nCASE: {name}")
            print(
                f"  skyborn trend / pmk slope : {sky['trend']} / {pmk_result.slope} "
                f"equal={trend_equal}"
            )
            print(
                f"  skyborn h     / pmk h     : {sky['h']} / {pmk_result.h} "
                f"equal={h_equal}"
            )
            print(
                f"  skyborn p     / pmk p     : {sky['p']} / {pmk_result.p} "
                f"equal={p_equal}"
            )
            print(
                f"  skyborn z     / pmk z     : {sky['z']} / {pmk_result.z} "
                f"equal={z_equal}"
            )
            print(
                f"  skyborn tau   / pmk Tau   : {sky['tau']} / {pmk_result.Tau} "
                f"equal={tau_equal}"
            )
            print(f"  all_equal = {all_equal}")


if __name__ == "__main__":
    main()
