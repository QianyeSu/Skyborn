import importlib.util
import sys
import types
from pathlib import Path

import pytest

MODULE_PATH = (
    Path(__file__).resolve().parents[1] / "src" / "skyborn" / "cdo" / "__init__.py"
)


def _load_cdo_shim(module_name: str):
    spec = importlib.util.spec_from_file_location(module_name, MODULE_PATH)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_cdo_namespace_reexports_skyborn_cdo(monkeypatch):
    fake = types.ModuleType("skyborn_cdo")

    class FakeCdo:
        pass

    class FakeCdoError(Exception):
        pass

    fake.Cdo = FakeCdo
    fake.CdoError = FakeCdoError
    fake.get_cdo_path = lambda: "/tmp/cdo"
    fake.get_cdo_version = lambda: "2.6.0"

    monkeypatch.setitem(sys.modules, "skyborn_cdo", fake)

    module = _load_cdo_shim("test_skyborn_cdo_shim_success")

    assert module.Cdo is FakeCdo
    assert module.CdoError is FakeCdoError
    assert module.get_cdo_path() == "/tmp/cdo"
    assert module.get_cdo_version() == "2.6.0"
    assert module.__all__ == ["Cdo", "CdoError", "get_cdo_path", "get_cdo_version"]


def test_cdo_namespace_raises_helpful_import_error(monkeypatch):
    monkeypatch.setitem(sys.modules, "skyborn_cdo", None)

    with pytest.raises(ImportError, match="skyborn-cdo"):
        _load_cdo_shim("test_skyborn_cdo_shim_missing")
