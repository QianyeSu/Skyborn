## Submission Recipe

This directory now tracks the `conda-forge/staged-recipes` oriented recipe.

- `meta.yaml`
  - uses the public GitHub source archive for the current `0.4.4` commit
    (the `v0.4.4` tag will be used after the release is created)
- `conda_build_config.yaml`
  - keeps only the compiler override needed to stay on a consistent GNU
    toolchain on Windows for this Meson + F2PY build
