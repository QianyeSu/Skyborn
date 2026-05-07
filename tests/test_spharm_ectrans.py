from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

REPO_SRC = Path(__file__).resolve().parents[1] / "src"
if str(REPO_SRC) not in sys.path:
    sys.path.insert(0, str(REPO_SRC))

from skyborn.spharm.ectrans_backend_api import scalar_synthesis_stub  # noqa: E402
from skyborn.spharm.reduced_gaussian import (  # noqa: E402
    ReducedGaussianGrid,
    ReducedGaussianSpharmt,
)


def test_reduced_gaussian_grid_basic_metadata():
    grid = ReducedGaussianGrid(np.array([20, 24, 28, 24, 20], dtype=np.int32))
    assert grid.ndgl == 5
    assert grid.ngptot == 116
    assert grid.max_nlon == 28


def test_reduced_gaussian_spharmt_repr_contains_key_sizes():
    backend = ReducedGaussianSpharmt(np.array([20, 24, 28, 24, 20], dtype=np.int32))
    text = repr(backend)
    assert "ndgl=5" in text
    assert "ngptot=116" in text


def test_reduced_gaussian_spharmt_scalar_prototype_returns_result():
    backend = ReducedGaussianSpharmt(np.array([20, 24, 28, 24, 20], dtype=np.int32))
    field = np.zeros((116,), dtype=np.float64)
    spec = backend.grdtospec(field)
    assert spec.shape == (15,)
    np.testing.assert_allclose(spec, 0.0, rtol=0.0, atol=1e-12)
    backend.close()


def test_reduced_gaussian_vector_prototype_matches_bridge_reference():
    backend = ReducedGaussianSpharmt(np.array([20, 24, 28, 24, 20], dtype=np.int32))
    bridge = backend._get_bridge_spharmt()

    spec = np.zeros((6,), dtype=np.complex128)
    spec[1] = 0.25 - 0.1j
    spec[2] = -0.2 + 0.05j
    spec[4] = 0.15 + 0.12j

    full_scalar = bridge.spectogrd(spec)
    packed_scalar = backend._bridge_to_packed_grid(full_scalar, "test packed scalar")
    rebuilt_full_scalar = backend._packed_to_bridge_grid(
        backend._validate_packed_grid_data(packed_scalar, "test packed scalar")[1],
        (),
    )
    np.testing.assert_allclose(rebuilt_full_scalar, full_scalar, rtol=0.0, atol=5e-8)

    u_expected, v_expected = bridge.getgrad(spec)
    packed_u = backend._bridge_to_packed_grid(u_expected, "test packed u")
    packed_v = backend._bridge_to_packed_grid(v_expected, "test packed v")

    vrt_expected, div_expected = bridge.getvrtdivspec(u_expected, v_expected, ntrunc=2)
    vrt_actual, div_actual = backend.getvrtdivspec(packed_u, packed_v, ntrunc=2)
    np.testing.assert_allclose(vrt_actual, vrt_expected, rtol=0.0, atol=1e-10)
    np.testing.assert_allclose(div_actual, div_expected, rtol=0.0, atol=1e-10)

    u_actual, v_actual = backend.getuv(vrt_expected, div_expected)
    np.testing.assert_allclose(u_actual, packed_u, rtol=0.0, atol=1e-10)
    np.testing.assert_allclose(v_actual, packed_v, rtol=0.0, atol=1e-10)

    grad_u_actual, grad_v_actual = backend.getgrad(spec)
    np.testing.assert_allclose(grad_u_actual, packed_u, rtol=0.0, atol=1e-10)
    np.testing.assert_allclose(grad_v_actual, packed_v, rtol=0.0, atol=1e-10)

    backend.close()


def test_reduced_gaussian_spharmt_vector_shape_validation():
    backend = ReducedGaussianSpharmt(np.array([20, 24, 28, 24, 20], dtype=np.int32))
    u = np.zeros((116,), dtype=np.float64)
    v = np.zeros((116, 2), dtype=np.float64)
    spec = np.zeros((15,), dtype=np.complex128)
    spec_2d = np.zeros((15, 2), dtype=np.complex128)

    try:
        backend.getvrtdivspec(u, v)
    except Exception as exc:
        assert "same shape" in str(exc)
    else:
        raise AssertionError("Expected same-shape validation error for vector analysis")

    try:
        backend.getuv(spec, spec_2d)
    except Exception as exc:
        assert "same shape" in str(exc)
    else:
        raise AssertionError("Expected same-shape validation error for wind synthesis")


def test_reduced_gaussian_scalar_ectrans_pack_roundtrip():
    backend = ReducedGaussianSpharmt(np.array([20, 24, 28, 24, 20], dtype=np.int32))
    spec = np.array(
        [
            1.0 + 2.0j,
            3.0 + 4.0j,
            5.0 + 6.0j,
            7.0 + 8.0j,
            9.0 + 10.0j,
            11.0 + 12.0j,
        ],
        dtype=np.complex128,
    )

    packed, ntrunc, extra_shape = backend._pack_scalar_spectrum_for_ectrans(
        spec, "test_pack"
    )
    restored = backend._unpack_scalar_spectrum_from_ectrans(packed, "test_unpack")

    assert ntrunc == 2
    assert extra_shape == ()
    assert packed.shape == (12,)
    np.testing.assert_allclose(packed[0::2], spec.real)
    np.testing.assert_allclose(packed[1::2], spec.imag)
    np.testing.assert_allclose(restored, spec)


def test_reduced_gaussian_validate_ectrans_scalar_packed_layout():
    backend = ReducedGaussianSpharmt(np.array([20, 24, 28, 24, 20], dtype=np.int32))
    packed = np.zeros((12, 3), dtype=np.float64)

    nt, ntrunc, normalized, extra_shape = backend._validate_ectrans_scalar_packed(
        packed, "packed_layout"
    )

    assert nt == 3
    assert ntrunc == 2
    assert normalized.shape == (12, 3)
    assert extra_shape == (3,)


def test_reduced_gaussian_native_setup_handle_lifecycle():
    backend = ReducedGaussianSpharmt(np.array([20, 24, 28, 24, 20], dtype=np.int32))

    handle = backend._require_setup_handle()
    assert handle is backend._require_setup_handle()

    backend.close()
    assert backend._setup_handle is None

    handle = backend._require_setup_handle()
    assert handle is not None
    backend.close()


def test_reduced_gaussian_native_setup_layout_metadata():
    backend = ReducedGaussianSpharmt(np.array([20, 24, 28, 24, 20], dtype=np.int32))

    description = backend._describe_setup()

    assert description["ndgl"] == 5
    assert description["ngptot"] == 116
    assert description["max_nlon"] == 28
    np.testing.assert_array_equal(
        description["nloen"],
        np.array([20, 24, 28, 24, 20], dtype=np.int32),
    )
    np.testing.assert_array_equal(
        description["lat_offsets"],
        np.array([0, 20, 44, 72, 96], dtype=np.int32),
    )
    assert description["mu"].shape == (5,)
    assert description["weights"].shape == (5,)
    assert np.all(np.diff(description["mu"]) < 0.0)
    assert np.all(description["weights"] > 0.0)
    np.testing.assert_allclose(description["weights"].sum(), 2.0, rtol=0.0, atol=1e-14)

    backend.close()


def test_reduced_gaussian_close_is_idempotent():
    backend = ReducedGaussianSpharmt(np.array([20, 24, 28, 24, 20], dtype=np.int32))
    backend.close()
    backend.close()


def test_reduced_gaussian_spectogrd_constant_field_prototype():
    backend = ReducedGaussianSpharmt(np.array([20, 24, 28, 24, 20], dtype=np.int32))
    spec = np.zeros((6,), dtype=np.complex128)
    spec[0] = np.sqrt(2.0)

    grid = backend.spectogrd(spec)

    assert grid.shape == (116,)
    np.testing.assert_allclose(grid, 1.0, rtol=0.0, atol=1e-6)
    backend.close()


def test_reduced_gaussian_native_scalar_synthesis_stub_returns_success():
    backend = ReducedGaussianSpharmt(np.array([20, 24, 28, 24, 20], dtype=np.int32))
    spec = np.zeros((6,), dtype=np.complex128)
    spec[0] = np.sqrt(2.0)
    packed_spec, _, _ = backend._pack_scalar_spectrum_for_ectrans(
        spec, "native stub test"
    )

    datagrid, ierror = scalar_synthesis_stub(
        backend._require_setup_handle(),
        np.asarray(packed_spec, dtype=np.float64),
    )

    assert ierror == 0
    assert datagrid.shape == (116,)
    np.testing.assert_allclose(datagrid, 1.0, rtol=0.0, atol=1e-6)
    backend.close()


def test_reduced_gaussian_grdtospec_constant_field_prototype():
    backend = ReducedGaussianSpharmt(np.array([20, 24, 28, 24, 20], dtype=np.int32))
    grid = np.ones((116,), dtype=np.float64)

    spec = backend.grdtospec(grid, ntrunc=2)

    assert spec.shape == (6,)
    np.testing.assert_allclose(spec[0], np.sqrt(2.0), rtol=0.0, atol=1e-6)
    np.testing.assert_allclose(spec[1:], 0.0, rtol=0.0, atol=1e-6)
    backend.close()


def test_reduced_gaussian_grdtospec_recovers_native_basis_vector():
    backend = ReducedGaussianSpharmt(np.array([20, 24, 28, 24, 20], dtype=np.int32))
    spec_in = np.zeros((6,), dtype=np.complex128)
    spec_in[4] = 1.0 + 0.0j

    grid = backend.spectogrd(spec_in)
    spec_out = backend.grdtospec(grid, ntrunc=2)

    np.testing.assert_allclose(spec_out, spec_in, rtol=0.0, atol=1e-5)
    backend.close()


def test_reduced_gaussian_scalar_constant_roundtrip_subset():
    backend = ReducedGaussianSpharmt(np.array([20, 24, 28, 24, 20], dtype=np.int32))
    spec_in = np.zeros((6,), dtype=np.complex128)
    spec_in[0] = np.sqrt(2.0)

    grid = backend.spectogrd(spec_in)
    spec_out = backend.grdtospec(grid, ntrunc=2)

    np.testing.assert_allclose(spec_out[0], spec_in[0], rtol=0.0, atol=1e-6)
    np.testing.assert_allclose(spec_out[1:], 0.0, rtol=0.0, atol=1e-6)
    backend.close()


def test_reduced_gaussian_spectogrd_multi_field_shape_and_values():
    backend = ReducedGaussianSpharmt(np.array([20, 24, 28, 24, 20], dtype=np.int32))
    spec = np.zeros((6, 2), dtype=np.complex128)
    spec[0, 0] = np.sqrt(2.0)
    spec[0, 1] = 2.0 * np.sqrt(2.0)

    grid = backend.spectogrd(spec)

    assert grid.shape == (116, 2)
    np.testing.assert_allclose(grid[:, 0], 1.0, rtol=0.0, atol=1e-6)
    np.testing.assert_allclose(grid[:, 1], 2.0, rtol=0.0, atol=1e-6)
    backend.close()


def test_reduced_gaussian_psichi_and_smoothing_match_bridge_reference():
    backend = ReducedGaussianSpharmt(np.array([20, 24, 28, 24, 20], dtype=np.int32))
    bridge = backend._get_bridge_spharmt()

    base_spec = np.zeros((6,), dtype=np.complex128)
    base_spec[1] = 0.25 - 0.1j
    base_spec[2] = -0.2 + 0.05j
    base_spec[4] = 0.15 + 0.12j

    full_u, full_v = bridge.getgrad(base_spec)
    packed_u = backend._bridge_to_packed_grid(full_u, "test packed u")
    packed_v = backend._bridge_to_packed_grid(full_v, "test packed v")

    psi_spec_expected = bridge.getpsispec(full_u, full_v, ntrunc=2)
    chi_spec_expected = bridge.getchispec(full_u, full_v, ntrunc=2)
    psi_spec_actual = backend.getpsispec(packed_u, packed_v, ntrunc=2)
    chi_spec_actual = backend.getchispec(packed_u, packed_v, ntrunc=2)
    np.testing.assert_allclose(psi_spec_actual, psi_spec_expected, rtol=0.0, atol=5e-8)
    np.testing.assert_allclose(chi_spec_actual, chi_spec_expected, rtol=0.0, atol=5e-8)

    psi_expected = backend._bridge_to_packed_grid(
        bridge.spectogrd(psi_spec_expected),
        "test packed psi",
    )
    chi_expected = backend._bridge_to_packed_grid(
        bridge.spectogrd(chi_spec_expected),
        "test packed chi",
    )
    psi_actual = backend.getpsi(packed_u, packed_v, ntrunc=2)
    chi_actual = backend.getchi(packed_u, packed_v, ntrunc=2)
    both_psi, both_chi = backend.getpsichi(packed_u, packed_v, ntrunc=2)
    np.testing.assert_allclose(psi_actual, psi_expected, rtol=0.0, atol=5e-8)
    np.testing.assert_allclose(chi_actual, chi_expected, rtol=0.0, atol=2e-7)
    np.testing.assert_allclose(both_psi, psi_expected, rtol=0.0, atol=5e-8)
    np.testing.assert_allclose(both_chi, chi_expected, rtol=0.0, atol=2e-7)

    smooth = np.array([1.0, 0.8, 0.6, 0.4, 0.2], dtype=np.float64)
    full_scalar = bridge.spectogrd(base_spec)
    packed_scalar = backend._bridge_to_packed_grid(
        full_scalar, "test packed scalar smooth"
    )
    smoothed_expected = backend._bridge_to_packed_grid(
        bridge.specsmooth(full_scalar, smooth),
        "test packed scalar smooth expected",
    )
    smoothed_actual = backend.specsmooth(packed_scalar, smooth)
    np.testing.assert_allclose(smoothed_actual, smoothed_expected, rtol=0.0, atol=5e-8)

    backend.close()
