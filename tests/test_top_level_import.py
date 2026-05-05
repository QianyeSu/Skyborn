"""Tests for the lightweight package root."""

import importlib
import sys


def _drop_skyborn_modules():
    for name in list(sys.modules):
        if name == "skyborn" or name.startswith("skyborn."):
            sys.modules.pop(name, None)


def test_import_skyborn_does_not_preload_heavy_submodules():
    _drop_skyborn_modules()

    skyborn = importlib.import_module("skyborn")

    assert skyborn.__version__
    assert "skyborn.calc" not in sys.modules
    assert "skyborn.gridfill" not in sys.modules
    assert "skyborn.plot" not in sys.modules
    assert "skyborn.spharm" not in sys.modules
    assert "skyborn.windspharm" not in sys.modules


def test_legacy_top_level_export_resolves_lazily():
    _drop_skyborn_modules()
    skyborn = importlib.import_module("skyborn")

    from skyborn.causality import liang_causality

    assert skyborn.liang_causality is liang_causality
    assert "skyborn.causality" in sys.modules


def test_common_top_level_exports_do_not_import_entire_calc_package():
    _drop_skyborn_modules()
    skyborn = importlib.import_module("skyborn")

    from skyborn.calc.calculations import linear_regression

    assert skyborn.linear_regression is linear_regression
    assert "skyborn.calc.calculations" in sys.modules
    assert "skyborn.calc.geostrophic" not in sys.modules
    assert "skyborn.spharm" not in sys.modules
    assert "skyborn.windspharm" not in sys.modules


def test_import_skyborn_calc_does_not_preload_extension_subpackages():
    _drop_skyborn_modules()

    calc = importlib.import_module("skyborn.calc")

    assert "linear_regression" in dir(calc)
    assert "skyborn.calc.geostrophic" not in sys.modules
    assert "skyborn.calc.GPI" not in sys.modules
    assert "skyborn.calc.troposphere" not in sys.modules
    assert "skyborn.spharm" not in sys.modules
    assert "skyborn.windspharm" not in sys.modules
