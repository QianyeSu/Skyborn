"""Tests for the lightweight package root."""

import importlib
import sys
from contextlib import contextmanager


@contextmanager
def _isolated_skyborn_imports():
    original_modules = {
        name: module
        for name, module in sys.modules.items()
        if name == "skyborn" or name.startswith("skyborn.")
    }

    for name in list(original_modules):
        sys.modules.pop(name, None)

    try:
        yield
    finally:
        for name in list(sys.modules):
            if name == "skyborn" or name.startswith("skyborn."):
                sys.modules.pop(name, None)
        sys.modules.update(original_modules)


def test_import_skyborn_does_not_preload_heavy_submodules():
    with _isolated_skyborn_imports():
        skyborn = importlib.import_module("skyborn")

        assert skyborn.__version__
        assert "skyborn.calc" not in sys.modules
        assert "skyborn.gridfill" not in sys.modules
        assert "skyborn.plot" not in sys.modules
        assert "skyborn.spharm" not in sys.modules
        assert "skyborn.windspharm" not in sys.modules


def test_legacy_top_level_export_resolves_lazily():
    with _isolated_skyborn_imports():
        skyborn = importlib.import_module("skyborn")

        from skyborn.causality import liang_causality

        assert skyborn.liang_causality is liang_causality
        assert "skyborn.causality" in sys.modules


def test_common_top_level_exports_do_not_import_entire_calc_package():
    with _isolated_skyborn_imports():
        skyborn = importlib.import_module("skyborn")

        from skyborn.calc.calculations import linear_regression

        assert skyborn.linear_regression is linear_regression
        assert "skyborn.calc.calculations" in sys.modules
        assert "skyborn.calc.geostrophic" not in sys.modules
        assert "skyborn.spharm" not in sys.modules
        assert "skyborn.windspharm" not in sys.modules


def test_import_skyborn_calc_does_not_preload_extension_subpackages():
    with _isolated_skyborn_imports():
        calc = importlib.import_module("skyborn.calc")

        assert "linear_regression" in dir(calc)
        assert "skyborn.calc.geostrophic" not in sys.modules
        assert "skyborn.calc.GPI" not in sys.modules
        assert "skyborn.calc.troposphere" not in sys.modules
        assert "skyborn.spharm" not in sys.modules
        assert "skyborn.windspharm" not in sys.modules


def test_calc_lazy_exports_resolve_submodules_and_objects():
    with _isolated_skyborn_imports():
        calc = importlib.import_module("skyborn.calc")

        geostrophic_module = calc.geostrophic
        assert geostrophic_module.__name__ == "skyborn.calc.geostrophic"
        assert "skyborn.calc.geostrophic" in sys.modules

        from skyborn.calc.calculations import pearson_correlation

        assert calc.pearson_correlation is pearson_correlation
        assert "skyborn.calc.calculations" in sys.modules


def test_calc_unknown_attribute_raises_attribute_error():
    with _isolated_skyborn_imports():
        calc = importlib.import_module("skyborn.calc")

        try:
            calc.this_does_not_exist
        except AttributeError as exc:
            assert "this_does_not_exist" in str(exc)
        else:
            raise AssertionError("Expected AttributeError for unknown calc attribute")
