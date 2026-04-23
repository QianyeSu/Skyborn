# Skyborn `conda-forge/staged-recipes` Submission Checklist

This checklist is specific to the current Skyborn `0.3.20` packaging state.

## 1. Prepare the staged-recipes working tree

1. Fork `conda-forge/staged-recipes`.
2. Create a branch from `main`.
3. Create a new directory:

```text
recipes/skyborn/
```

4. Copy these files from this repository into that new directory:

```text
D:\Skyborn\recipe\meta.yaml
D:\Skyborn\recipe\conda_build_config.yaml
```

## 2. Keep the submission recipe, not the local recipe

Use the submission-oriented files under:

```text
D:\Skyborn\recipe\
```

Do not submit the local debug recipe under:

```text
D:\Skyborn\recipe-local\
```

That local recipe uses `source.path: ..` and is only for local Windows
validation.

## 3. Source archive details for `0.3.20`

Current release source in `meta.yaml`:

```yaml
source:
  url: https://github.com/QianyeSu/Skyborn/archive/refs/tags/v0.3.20.tar.gz
  sha256: f1d6b063283a50ae58181d6ac25627a0a25b685c782232d21be58edf424b7c95
```

If the source artifact changes, regenerate `sha256` before submission.

## 4. Recipe shape currently chosen

This package currently uses:

- `meta.yaml`, not `recipe.yaml`
- compiled extensions, so this is **not** `noarch: python`
- explicit compiler macros:
  - `{{ compiler('c') }}`
  - `{{ compiler('fortran') }}`

This matches the current package structure better than a simple noarch recipe.

## 5. Windows-specific note for this package

Skyborn currently needs a consistent GNU-style toolchain on Windows because the
package builds multiple Meson + F2PY extensions and the mixed MSVC + GNU
combination was not stable for local validation.

The submission recipe therefore keeps:

```yaml
c_compiler:
  - gcc

cxx_compiler:
  - gxx

fortran_compiler:
  - gfortran
```

inside `conda_build_config.yaml`.

## 6. Minimum checks before opening the PR

Local checks already completed in this repository:

- `conda render` passed for the submission recipe
- `conda build` passed locally on Windows using the local validation recipe
- `conda install --use-local skyborn=0.3.20` succeeded
- smoke import passed:
  - `skyborn`
  - `skyborn.calc`
  - `skyborn.interp`
  - `skyborn.spharm`

Recommended last pass before PR:

```bash
conda render recipe -m recipe/conda_build_config.yaml -c conda-forge --override-channels --python 3.12 --numpy 1.26
```

Inside the staged-recipes fork, this becomes:

```bash
conda render recipes/skyborn -m recipes/skyborn/conda_build_config.yaml -c conda-forge --override-channels --python 3.12 --numpy 1.26
```

## 7. What to expect from the PR

The `staged-recipes` CI will:

- lint the recipe
- try to build it on multiple platforms
- review source/license/test metadata

For Skyborn, the likely review focus is:

- compiled-extension complexity
- Windows toolchain choice
- import tests being sufficient for a first submission
- dependency breadth

## 8. What to say in the PR description

Useful concise points:

- Skyborn is an atmospheric-science package with Python, Fortran, and Cython
  extensions.
- The recipe uses the public GitHub `v0.3.20` tag tarball.
- The package is not `noarch` because it builds compiled extensions.
- Local Windows validation was completed with Meson + F2PY extensions and a
  successful `conda install --use-local` smoke test.

## 9. After merge

After the first submission is merged:

- `conda-forge` will create the `skyborn-feedstock`
- future version bumps should go to the feedstock, not back to
  `staged-recipes`

## 10. Current authoritative references

- Staged-recipes:
  - https://conda-forge.org/docs/maintainer/understanding_conda_forge/staged_recipes/
- Contributing packages:
  - https://conda-forge.org/docs/maintainer/adding_pkgs/
- Maintaining packages:
  - https://conda-forge.org/docs/maintainer/updating_pkgs/
- Infrastructure and compiler defaults:
  - https://conda-forge.org/docs/maintainer/infrastructure
