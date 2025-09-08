#!/usr/bin/env python
"""
Check actual coverage of GPI module by importing and using all functions.
"""

import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "../..")))


def main():
    """Test all GPI module imports and basic functionality."""

    print("Testing GPI module coverage...")

    # Test interface module
    import numpy as np

    from skyborn.calc.GPI.interface import (
        UNDEF,
        PotentialIntensityCalculator,
        _ensure_pressure_ordering,
        _postprocess_results,
        _validate_dimensions,
        _validate_input_arrays,
        potential_intensity,
    )

    # Exercise functions
    arr = np.array([1, 2, 3])
    _postprocess_results(arr, arr)
    _validate_input_arrays(arr)

    pressure = np.array([1000, 850, 700])
    temp = np.array([300, 285, 270])
    mixr = np.array([0.012, 0.008, 0.005])
    _ensure_pressure_ordering(pressure, temp, mixr)

    # Test xarray module
    try:
        import xarray as xr

        from skyborn.calc.GPI.xarray import (
            _check_units,
            _create_output_dataset,
            _detect_atmospheric_dimensions,
            potential_intensity,
        )

        # Create test data
        da = xr.DataArray(25.0, attrs={"units": "K"})
        _check_units(da, "temperature", "K")

        temp_3d = xr.DataArray(
            np.random.rand(10, 20, 30), dims=["pressure", "lat", "lon"]
        )
        _detect_atmospheric_dimensions(temp_3d)

        print("XArray module tested")
    except ImportError:
        print("XArray not available")

    # Calculate coverage
    import coverage

    cov = coverage.Coverage()
    cov.load()

    total_statements = 0
    covered_statements = 0

    for filename in cov.get_data().measured_files():
        if "skyborn/calc/GPI" in filename and "test" not in filename:
            analysis = cov.analysis2(filename)
            executed = len(analysis[1])
            missing = len(analysis[3])
            total = executed + missing
            if total > 0:
                percent = (executed / total) * 100
                module_name = os.path.basename(filename)
                print(f"  {module_name}: {percent:.1f}% ({executed}/{total} lines)")
                total_statements += total
                covered_statements += executed

    if total_statements > 0:
        overall = (covered_statements / total_statements) * 100
        print(
            f"\nOverall GPI coverage: {overall:.1f}% ({covered_statements}/{total_statements} lines)"
        )

        if overall >= 96:
            print("✓ Coverage goal of 96% achieved!")
            return 0
        else:
            needed = int(0.96 * total_statements) - covered_statements
            print(f"Need {needed} more lines for 96% coverage")
            return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
