"""Focused regression tests for the compiled geostrophic backend module."""

from __future__ import annotations

import numpy as np
import pytest

backend = pytest.importorskip("skyborn.calc.geostrophic.geostrophicwind")


def test_compiled_backend_accepts_keyword_signature_2d():
    """The direct backend should preserve the keyword contract used by core.py."""

    nlat = 5
    mlon = 8
    z = np.arange(nlat * mlon, dtype=np.float32).reshape(nlat, mlon)
    glon = np.linspace(0.0, 315.0, mlon, dtype=np.float64)
    glat = np.linspace(-60.0, 60.0, nlat, dtype=np.float64)

    ug, vg = backend.z2geouv(
        z,
        zmsg=-999.0,
        glon=glon,
        glat=glat,
        iopt=1,
        nlat=nlat,
        mlon=mlon,
    )

    assert ug.shape == (nlat, mlon)
    assert vg.shape == (nlat, mlon)
    assert ug.dtype == np.float64
    assert vg.dtype == np.float64


def test_compiled_backend_accepts_keyword_signature_3d():
    """The direct backend should accept shape-validation keywords for 3D calls."""

    nlat = 4
    mlon = 6
    n3rd = 3
    z = np.arange(nlat * mlon * n3rd, dtype=np.float32).reshape(nlat, mlon, n3rd)
    glon = np.linspace(0.0, 300.0, mlon, dtype=np.float64)
    glat = np.linspace(-45.0, 45.0, nlat, dtype=np.float64)

    ug, vg = backend.z2geouv_3d(
        z,
        zmsg=-999.0,
        glon=glon,
        glat=glat,
        iopt=1,
        nlat=nlat,
        mlon=mlon,
        n3rd=n3rd,
    )

    assert ug.shape == (nlat, mlon, n3rd)
    assert vg.shape == (nlat, mlon, n3rd)
    assert ug.dtype == np.float64
    assert vg.dtype == np.float64
