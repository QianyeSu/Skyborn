"""
Accuracy comparison between new Fortran 90 implementation and original Fortran 77 code.

Compares:
1. Weight computation
2. Filter response
3. Filtering operations

Against the reference implementation from:
https://github.com/QianyeSu/Meteorology-Fortran-Functions/blob/main/fortran/filtrx.f
"""

import numpy as np

from skyborn.calc.filter import lanczos_highpass, lanczos_lowpass, lanczos_weights


def reference_fortran77_weights(nwt, fca, nsigma=1.0):
    """
    Reference implementation matching original DFILWTQ Fortran 77 code.

    This is a direct translation of the original algorithm.
    """
    pi = np.pi
    twopi = 2.0 * pi
    nw = (nwt - 1) // 2
    fnw = float(nw)

    arg = twopi * fca

    wt = np.zeros(nwt)

    # Central weight
    wt[0] = 2.0 * fca

    # Compute weights with sigma factor
    for i in range(1, nw + 1):
        sinx = np.sin(arg * float(i)) / (pi * float(i))
        siny = fnw * np.sin(float(i) * pi / fnw) / (float(i) * pi)
        wt[i] = sinx * (siny**nsigma)

    # Normalize weights
    sum_wt = wt[0]
    for i in range(1, nw + 1):
        sum_wt = sum_wt + 2.0 * wt[i]

    for i in range(nw + 1):
        wt[i] = wt[i] / sum_wt

    # Make symmetric (expand to full nwt length)
    work = wt[1 : nw + 1].copy()
    wt_full = np.zeros(nwt)
    wt_full[nw] = wt[0]  # Center
    for n in range(nw):
        wt_full[n] = work[nw - 1 - n]
        wt_full[nw + 1 + n] = work[n]

    return wt_full


def reference_fortran77_highpass(nwt, fca, nsigma=1.0):
    """
    Reference high-pass filter matching original DFILTRQ with IHP=1.
    """
    # Get low-pass weights (before symmetrization)
    pi = np.pi
    twopi = 2.0 * pi
    nw = (nwt - 1) // 2
    fnw = float(nw)

    arg = twopi * fca

    wt = np.zeros(nw + 1)
    wt[0] = 2.0 * fca

    for i in range(1, nw + 1):
        sinx = np.sin(arg * float(i)) / (pi * float(i))
        siny = fnw * np.sin(float(i) * pi / fnw) / (float(i) * pi)
        wt[i] = sinx * (siny**nsigma)

    # Normalize
    sum_wt = wt[0]
    for i in range(1, nw + 1):
        sum_wt = sum_wt + 2.0 * wt[i]

    for i in range(nw + 1):
        wt[i] = wt[i] / sum_wt

    # High-pass transformation
    wt[0] = 1.0 - wt[0]
    for i in range(1, nw + 1):
        wt[i] = -wt[i]

    # Make symmetric
    work = wt[1 : nw + 1].copy()
    wt_full = np.zeros(nwt)
    wt_full[nw] = wt[0]
    for n in range(nw):
        wt_full[n] = work[nw - 1 - n]
        wt_full[nw + 1 + n] = work[n]

    return wt_full


def compute_response_function(weights):
    """
    Compute frequency response function for given weights.

    Matches the response computation from DFILTRQ.
    """
    nwt = len(weights)
    nw = (nwt - 1) // 2
    nf = 2 * nwt - 1

    freq = np.zeros(nf)
    resp = np.zeros(nf)

    twopi = 2.0 * np.pi
    frqint = 0.5 / float(nf - 1)

    # DC response
    freq[0] = 0.0
    sum_resp = 0.0
    for j in range(1, nw + 1):
        sum_resp = sum_resp + 2.0 * weights[nw + j]
    resp[0] = sum_resp + weights[nw]

    # Frequency response
    for i in range(1, nf):
        sum_resp = 0.0
        freq[i] = float(i) * frqint
        for j in range(1, nw + 1):
            sum_resp = sum_resp + weights[nw + j] * np.cos(twopi * freq[i] * float(j))
        resp[i] = weights[nw] + 2.0 * sum_resp

    return freq, resp


def test_weight_comparison():
    """Compare weight computation with original Fortran 77."""
    print("=" * 80)
    print("WEIGHT COMPUTATION COMPARISON")
    print("=" * 80)

    test_cases = [
        # (nwt, fca, description)
        (61, 1 / 30.0, "30-day low-pass, 61 weights"),
        (81, 0.05, "0.05 cutoff, 81 weights"),
        (121, 1 / 365.0, "365-day low-pass, 121 weights"),
        (41, 0.2, "0.2 cutoff, 41 weights"),
    ]

    max_error = 0.0

    for nwt, fca, desc in test_cases:
        print(f"\nTest: {desc}")

        # Our Fortran 90 implementation
        wt_f90 = lanczos_weights(fca, nwt, "low")

        # Reference Fortran 77 implementation
        wt_f77 = reference_fortran77_weights(nwt, fca, nsigma=1.0)

        # Compare
        abs_error = np.abs(wt_f90 - wt_f77)
        max_abs_err = abs_error.max()
        mean_abs_err = abs_error.mean()
        rel_error = abs_error / (np.abs(wt_f77) + 1e-15)
        max_rel_err = rel_error.max()

        max_error = max(max_error, max_abs_err)

        print(f"  Max absolute error: {max_abs_err:.3e}")
        print(f"  Mean absolute error: {mean_abs_err:.3e}")
        print(f"  Max relative error: {max_rel_err:.3e}")
        print(f"  F90 sum: {wt_f90.sum():.15f}")
        print(f"  F77 sum: {wt_f77.sum():.15f}")

        # Compute response functions
        freq_f90, resp_f90 = compute_response_function(wt_f90)
        freq_f77, resp_f77 = compute_response_function(wt_f77)

        resp_error = np.abs(resp_f90 - resp_f77).max()
        print(f"  Response function max error: {resp_error:.3e}")

        # Check key properties
        assert abs(wt_f90.sum() - 1.0) < 1e-12, "F90 weights not normalized"
        assert abs(wt_f77.sum() - 1.0) < 1e-12, "F77 weights not normalized"
        assert max_abs_err < 1e-13, f"Weight error too large: {max_abs_err}"

    print("\n" + "=" * 80)
    print(f"OVERALL MAX ERROR: {max_error:.3e}")
    print("[PASS] Low-pass weight computation matches original Fortran 77")
    print("=" * 80)

    return True


def test_highpass_comparison():
    """Compare high-pass filter with original Fortran 77."""
    print("\n" + "=" * 80)
    print("HIGH-PASS FILTER COMPARISON")
    print("=" * 80)

    test_cases = [
        (61, 1 / 30.0, "30-day high-pass"),
        (81, 0.1, "0.1 cutoff high-pass"),
    ]

    max_error = 0.0

    for nwt, fca, desc in test_cases:
        print(f"\nTest: {desc}")

        # Our Fortran 90 implementation
        wt_f90 = lanczos_weights(fca, nwt, "high")

        # Reference Fortran 77 implementation
        wt_f77 = reference_fortran77_highpass(nwt, fca, nsigma=1.0)

        # Compare
        abs_error = np.abs(wt_f90 - wt_f77)
        max_abs_err = abs_error.max()
        mean_abs_err = abs_error.mean()

        max_error = max(max_error, max_abs_err)

        print(f"  Max absolute error: {max_abs_err:.3e}")
        print(f"  Mean absolute error: {mean_abs_err:.3e}")
        print(f"  F90 sum: {wt_f90.sum():.15f}")
        print(f"  F77 sum: {wt_f77.sum():.15f}")

        # High-pass should sum to zero
        assert abs(wt_f90.sum()) < 1e-12, "F90 high-pass not zero-sum"
        assert abs(wt_f77.sum()) < 1e-12, "F77 high-pass not zero-sum"
        assert max_abs_err < 1e-13, f"High-pass error too large: {max_abs_err}"

    print("\n" + "=" * 80)
    print(f"OVERALL MAX ERROR: {max_error:.3e}")
    print("[PASS] High-pass filter computation matches original Fortran 77")
    print("=" * 80)

    return True


def test_nsigma_parameter():
    """
    Test different nsigma values (Lanczos window power).

    Note: Our implementation uses nsigma=1 (standard Lanczos).
    The original code allows nsigma as parameter.
    """
    print("\n" + "=" * 80)
    print("NSIGMA PARAMETER TEST")
    print("=" * 80)

    nwt = 61
    fca = 0.1

    # Our implementation (nsigma=1 hardcoded)
    wt_f90 = lanczos_weights(fca, nwt, "low")

    # Reference with nsigma=1
    wt_f77_ns1 = reference_fortran77_weights(nwt, fca, nsigma=1.0)

    # Reference with nsigma=2 (sharper cutoff)
    wt_f77_ns2 = reference_fortran77_weights(nwt, fca, nsigma=2.0)

    error_ns1 = np.abs(wt_f90 - wt_f77_ns1).max()
    error_ns2 = np.abs(wt_f90 - wt_f77_ns2).max()

    print(f"\nOur F90 vs F77 nsigma=1: {error_ns1:.3e}")
    print(f"Our F90 vs F77 nsigma=2: {error_ns2:.3e}")

    assert error_ns1 < 1e-13, "Should match nsigma=1"
    assert error_ns2 > 1e-3, "Should NOT match nsigma=2"

    print("\n[PASS] Our implementation correctly uses nsigma=1 (standard Lanczos)")
    print("=" * 80)

    return True


def test_filtering_operation():
    """Test actual filtering operation accuracy."""
    print("\n" + "=" * 80)
    print("FILTERING OPERATION TEST")
    print("=" * 80)

    # Create test signal
    np.random.seed(42)
    n = 1000
    t = np.arange(n)

    # Signal = low freq + high freq noise
    signal = np.sin(2 * np.pi * t / 100) + 0.3 * np.random.randn(n)

    # Filter with our implementation
    nwt = 61
    fca = 0.05
    filtered_f90 = lanczos_lowpass(signal, cutoff_freq=fca, window=nwt)

    # Manual convolution with reference weights
    wt_f77 = reference_fortran77_weights(nwt, fca)
    nw = (nwt - 1) // 2

    # Apply convolution manually (reflect boundary)
    filtered_manual = np.zeros(n)
    for i in range(n):
        sum_val = 0.0
        for j in range(nwt):
            idx = i + j - nw
            # Reflect boundary
            if idx < 0:
                idx = -idx
            elif idx >= n:
                idx = 2 * n - idx - 2
            sum_val += signal[idx] * wt_f77[j]
        filtered_manual[i] = sum_val

    # Compare (avoid edges due to boundary effects)
    center = slice(nw + 10, n - nw - 10)
    error = np.abs(filtered_f90[center] - filtered_manual[center])
    max_err = error.max()
    mean_err = error.mean()

    print(f"\nMax filtering error (center region): {max_err:.3e}")
    print(f"Mean filtering error (center region): {mean_err:.3e}")

    assert max_err < 1e-10, f"Filtering error too large: {max_err}"

    print("\n[PASS] Filtering operation matches reference implementation")
    print("=" * 80)

    return True


if __name__ == "__main__":
    print("\n")
    print("╔" + "═" * 78 + "╗")
    print("║" + " " * 78 + "║")
    print(
        "║"
        + " " * 10
        + "LANCZOS FILTER: FORTRAN 90 vs FORTRAN 77 COMPARISON"
        + " " * 16
        + "║"
    )
    print("║" + " " * 78 + "║")
    print(
        "║"
        + " " * 5
        + "Comparing against original reference implementation"
        + " " * 21
        + "║"
    )
    print(
        "║" + " " * 5 + "from Meteorology-Fortran-Functions repository" + " " * 26 + "║"
    )
    print("║" + " " * 78 + "║")
    print("╚" + "═" * 78 + "╝")
    print("\n")

    try:
        test_weight_comparison()
        test_highpass_comparison()
        test_nsigma_parameter()
        test_filtering_operation()

        print("\n")
        print("╔" + "═" * 78 + "╗")
        print("║" + " " * 78 + "║")
        print(
            "║"
            + " " * 15
            + "ALL COMPARISON TESTS PASSED SUCCESSFULLY!"
            + " " * 21
            + "║"
        )
        print("║" + " " * 78 + "║")
        print(
            "║"
            + " " * 5
            + "Fortran 90 implementation matches Fortran 77 reference"
            + " " * 18
            + "║"
        )
        print(
            "║"
            + " " * 5
            + "to machine precision (< 1e-13 absolute error)"
            + " " * 28
            + "║"
        )
        print("║" + " " * 78 + "║")
        print("╚" + "═" * 78 + "╝")
        print("\n")

    except AssertionError as e:
        print("\n")
        print("╔" + "═" * 78 + "╗")
        print("║" + " " * 78 + "║")
        print("║" + " " * 28 + "COMPARISON TEST FAILED" + " " * 28 + "║")
        print("║" + " " * 78 + "║")
        print("╚" + "═" * 78 + "╝")
        print(f"\nError: {e}\n")
        raise
