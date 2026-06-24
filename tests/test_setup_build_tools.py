"""Regression tests for build-tool command resolution."""

from __future__ import annotations

import importlib.util
from pathlib import Path


def _load_setup_module():
    setup_path = Path(__file__).resolve().parents[1] / "setup.py"
    spec = importlib.util.spec_from_file_location("skyborn_setup", setup_path)
    if spec is None or spec.loader is None:
        raise RuntimeError("Unable to load setup.py for testing")

    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_ninja_command_prefers_executable_when_available(monkeypatch):
    setup_module = _load_setup_module()
    monkeypatch.setattr(setup_module.shutil, "which", lambda name: "/opt/bin/ninja")

    assert setup_module.MesonBuildExt._ninja_command("--version") == [
        "/opt/bin/ninja",
        "--version",
    ]


def test_ninja_command_falls_back_to_python_module_when_needed(monkeypatch):
    setup_module = _load_setup_module()
    monkeypatch.setattr(setup_module.shutil, "which", lambda name: None)
    monkeypatch.setattr(setup_module.sys, "executable", "/opt/conda/bin/python")

    assert setup_module.MesonBuildExt._ninja_command("--version") == [
        "/opt/conda/bin/python",
        "-m",
        "ninja",
        "--version",
    ]


def test_meson_command_prefers_executable_when_available(monkeypatch):
    setup_module = _load_setup_module()
    monkeypatch.setattr(setup_module.shutil, "which", lambda name: "/opt/bin/meson")

    assert setup_module.MesonBuildExt._meson_command("--version") == [
        "/opt/bin/meson",
        "--version",
    ]


def test_meson_command_falls_back_to_python_module_when_needed(monkeypatch):
    setup_module = _load_setup_module()
    monkeypatch.setattr(setup_module.shutil, "which", lambda name: None)
    monkeypatch.setattr(setup_module.sys, "executable", "/opt/conda/bin/python")

    assert setup_module.MesonBuildExt._meson_command("--version") == [
        "/opt/conda/bin/python",
        "-m",
        "mesonbuild.mesonmain",
        "--version",
    ]
