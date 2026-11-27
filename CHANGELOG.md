# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

Planned as 1.0.0.

### Added

- VACF sampling in the generated PoreAna analysis scripts.
- `CITATION.cff`, which drives GitHub's "Cite this repository" button.
- CodeQL code scanning and Dependabot updates for GitHub Actions and uv.

### Changed

- `construct.py`: extracted `_write_positions(path, positions)`, which writes
  an (N, 3) NumPy position array to a GROMACS position file, replacing repeated
  inline loops in all molecule insertion paths.
- `construct.py`: extracted `_gmx_insert_cmd(folder_gro, file_box, mol, nmol,
  dr, pos_file)`, which builds the GROMACS `insert-molecules` command with the
  optional `-dr` and `-ip` flags for the pore, box and slit code paths.
- `construct.py`: the position files iterate the pore shapes by name instead of
  relying on the shapes being the first keys of the pore yml.
- `actuate.py`: the shared folder and file names (`_gro/`, `_top/`, `_mdp/`,
  `box.gro`, `topol.top`, `index.ndx`) are instance attributes set once in
  `Actuate.__init__`.
- `simulate.py`: the Jinja2 rendering of the analysis scripts collects the
  molecule metadata in one loop and selects the template in one branch.
- The analysis templates `auto_dens.py` and `auto_dens_box.py` are trimmed.
- `tests/data/forhlr.sh` loads GROMACS 2024.3 instead of 2016.5 and requests
  GPU nodes with `#SBATCH --gres`.
- Packaging moved from `setup.py` and `MANIFEST.in` to a single PEP 621
  `pyproject.toml` with the version taken from `poresim.__version__`, and
  from `requirements.txt` to PEP 735 dependency groups (`dev`, `docs`)
  resolved with [uv](https://docs.astral.sh/uv/).
- Python 3.12 or newer is required; 3.12 and 3.13 are tested.
- CI runs on uv: ruff lint and format check, the full test suite with
  coverage, a package build with `twine check`, the docs build with
  `sphinx-build -W` and its deployment to GitHub Pages, and a weekly
  `pip-audit` of the runtime dependencies.
- Documentation sources moved from `docsrc/` via `docs/` to `docs/source/`,
  converted from reStructuredText to MyST Markdown, with the furo theme and an
  API reference generated with autosummary; `docs/DESIGN.md` documents the
  amber colour palette.
- The test suite moved from unittest to pytest and is split into
  `test_unit.py` and `test_integration.py`, with new unit tests for the `Box`
  getters and setters, `add_mol` validation, `add_topol`, `_write_positions`
  and `_gmx_insert_cmd`.
- The code base is formatted with `ruff format` and passes `ruff check`.

### Removed

- The unused `from re import A` import in `actuate.py`.
- The unittest suite `tests/test_simple.py`, a duplicate of the pytest suite,
  and `tests/venv.txt`.

## [0.3.0] - 2025-10-17

### Changed

- New version due to a change of GitHub organisation.

## [0.2.0] - 2024-11-02

### Added

- Targeted placement of molecules in the pore.
- Molecules near the wall of a slit pore, to build a layer on the wall.
- Insertion of the molecules with GROMACS position files.

### Changed

- Update for PoreMS 0.3.0 and PoreAna 0.2.3.

## [0.1.2] - 2023-09-23

### Added

- Pores can be filled with mixtures, with a target density for more than one
  molecule.

## [0.1.1] - 2022-02-07

### Changed

- Update for PoreMS 0.2.5 and PoreAna 0.2.2.
- Improved templates.
- Updated simulation docs.

## [0.1.0] - 2021-03-08

- Initial version.
