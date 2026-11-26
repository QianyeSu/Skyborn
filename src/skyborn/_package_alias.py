"""Zero-copy subpackage aliasing for extracted Skyborn packages."""

from __future__ import annotations

import importlib
import importlib.util
import sys
from importlib.abc import Loader, MetaPathFinder
from importlib.machinery import ModuleSpec
from types import ModuleType
from typing import Dict


class _SubpackageAliasLoader(Loader):
    """Loader that maps ``alias.*`` directly onto existing ``target.*`` modules."""

    def __init__(self, alias_name: str, target_name: str) -> None:
        self.alias_name = alias_name
        self.target_name = target_name
        self._saved_attrs: Dict[int, tuple[str, str | None, object]] = {}

    def create_module(self, spec: ModuleSpec) -> ModuleType:
        module = importlib.import_module(self.target_name)
        self._saved_attrs[id(module)] = (
            module.__name__,
            module.__package__,
            module.__spec__,
        )
        return module

    def exec_module(self, module: ModuleType) -> None:
        saved = self._saved_attrs.pop(id(module), None)
        if saved is not None:
            orig_name, orig_package, orig_spec = saved
            module.__name__ = orig_name
            module.__package__ = orig_package
            module.__spec__ = orig_spec
        sys.modules[self.alias_name] = module


class _SubpackageAliasFinder(MetaPathFinder):
    """Meta-path finder redirecting legacy package subtrees/modules to external packages."""

    def __init__(self) -> None:
        self.mappings: Dict[str, str] = {}
        self.exact_mappings: Dict[str, str] = {}

    def register(self, alias_prefix: str, target_prefix: str) -> None:
        self.mappings[alias_prefix] = target_prefix

    def register_exact(self, alias_name: str, target_name: str) -> None:
        self.exact_mappings[alias_name] = target_name

    def find_spec(self, fullname: str, path=None, target=None) -> ModuleSpec | None:
        target_name = self.exact_mappings.get(fullname)
        if target_name is None:
            for alias_prefix, target_prefix in self.mappings.items():
                prefix = f"{alias_prefix}."
                if fullname.startswith(prefix):
                    suffix = fullname[len(prefix) :]
                    target_name = f"{target_prefix}.{suffix}"
                    break
        if target_name is None:
            return None

        target_spec = importlib.util.find_spec(target_name)
        if target_spec is None:
            return None
        loader = _SubpackageAliasLoader(fullname, target_name)
        is_pkg = target_spec.submodule_search_locations is not None
        spec = ModuleSpec(fullname, loader, is_package=is_pkg)
        if is_pkg:
            spec.submodule_search_locations = list(
                target_spec.submodule_search_locations
            )
        return spec


_ALIAS_FINDER = _SubpackageAliasFinder()


def register_subpackage_alias(alias_prefix: str, target_prefix: str) -> None:
    """Register ``alias_prefix.*`` to resolve transparently to ``target_prefix.*``."""
    _ALIAS_FINDER.register(alias_prefix, target_prefix)
    if _ALIAS_FINDER not in sys.meta_path:
        sys.meta_path.insert(0, _ALIAS_FINDER)


def register_module_alias(alias_name: str, target_name: str) -> None:
    """Register exact module ``alias_name`` to resolve to ``target_name``."""
    _ALIAS_FINDER.register_exact(alias_name, target_name)
    if _ALIAS_FINDER not in sys.meta_path:
        sys.meta_path.insert(0, _ALIAS_FINDER)
