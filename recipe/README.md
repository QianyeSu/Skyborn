## Submission Recipe

This directory now tracks the `conda-forge/staged-recipes` oriented recipe.

- `meta.yaml`
  - uses the public GitHub tag tarball for `v0.3.20`
- `conda_build_config.yaml`
  - keeps only the compiler override needed to stay on a consistent GNU
    toolchain on Windows for this Meson + F2PY build

For local path-based validation inside the repository, use
`D:\Skyborn\recipe-local` instead.
