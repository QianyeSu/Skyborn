## Summary

- This PR adds `skyborn` as a new `conda-forge` package.
- Skyborn is an atmospheric-science utilities package with compiled Python,
  Fortran, and Cython extensions.
- The current submission targets version `0.3.20`.

## Upstream

- Home: https://github.com/QianyeSu/Skyborn
- Documentation: https://skyborn.readthedocs.io/
- Source used by this recipe:
  - `https://github.com/QianyeSu/Skyborn/archive/refs/tags/v0.3.20.tar.gz`

## Packaging Notes

- This package is not `noarch: python` because it builds compiled extensions.
- The recipe uses Meson + F2PY and keeps explicit C / Fortran compiler macros.
- Local validation was completed on Windows with the `conda-forge` GNU
  toolchain.

## Local Validation Snapshot

- `conda build recipe-local -c conda-forge`
- `conda build recipe-local --output -c conda-forge`
- `conda install --use-local skyborn=0.3.20`
- Import smoke checks passed for:
  - `skyborn`
  - `skyborn.calc`
  - `skyborn.interp`
  - `skyborn.spharm`

## Package Content Notes

- The packaging cleanup for this release removed development-only artifacts
  from the install payload, including:
  - `*.dll.a`
  - `*.pyf`
  - `meson.build`
  - `gridfill/_gridfill.pyx`

## Maintainer

- GitHub: `@QianyeSu`
