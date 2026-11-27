# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

Planned as 1.0.0.

### Added

- VACF sampling in the generated PoreAna analysis scripts
  (`auto_dens.py`, `auto_dens_box.py`), with the correlation length, the new
  time origin, the sample step and the frame length as editable inputs at the
  top of the script.
- `__version__` in the `poresim` package, read by the build and the docs.
- `CITATION.cff`, which drives GitHub's "Cite this repository" button.
- `docs/DESIGN.md`, documenting the amber colour palette of the logo and docs.
- GitHub Actions workflows:
  - `ci.yml`: ruff lint and format check, the test suite with coverage on
    Python 3.12 and 3.13 with upload to Codecov, and a package build checked
    with `twine check`;
  - `docs.yml`: the documentation built with `sphinx-build -W` and deployed to
    GitHub Pages from `main`;
  - `security.yml`: a `pip-audit` of the runtime dependencies on every change
    and weekly;
  - `codeql.yml`: CodeQL code scanning;
  - Dependabot updates for GitHub Actions and uv.
- Tests:
  - `tests/test_unit.py` with fast tests for the `Box` getters and setters,
    the `add_mol` validation, `add_topol` for all types, the topology and itp
    generation, `Construct._write_positions` byte by byte and
    `Construct._gmx_insert_cmd` with all flag combinations;
  - `tests/test_integration.py` with full simulation folders for pore, box,
    benchmark and two-phase systems;
  - `tests/test_construct.py` with one test per placement branch of the box
    construction: fixed numbers and fills in the reservoir, pore, wall and
    box, slit and cylinder shapes, box systems with and without areas.

  Line coverage is 100 percent, branch coverage 99 percent.

### Changed

- `construct.py`: extracted `_write_positions(path, positions)`, which writes
  an (N, 3) NumPy position array to a GROMACS position file and replaces the
  repeated inline loops of all molecule insertion paths.
- `construct.py`: extracted `_gmx_insert_cmd(folder_gro, file_box, mol, nmol,
  dr, pos_file)`, which builds the GROMACS `insert-molecules` command with the
  optional `-dr` and `-ip` flags for the pore, box and slit code paths.
- `construct.py`: the position files iterate the pore shapes by name instead of
  relying on the shapes being the first keys of the pore yml.
- `actuate.py`: the shared folder and file names (`_gro/`, `_top/`, `_mdp/`,
  `box.gro`, `topol.top`, `index.ndx`) are instance attributes set once in
  `Actuate.__init__` and used by `_equilibration` and `_simulation`.
- `simulate.py`: the Jinja2 rendering of the analysis scripts collects the
  molecule metadata in a single loop, derives the fill flag once and selects
  the template in one branch.
- The analysis templates `auto_dens.py` and `auto_dens_box.py` are trimmed and
  pass the VACF direction as an axis index (`direction=2`).
- The `sort.py` template drops its commented-out ionic liquid code.
- `Benchmark` calls `super().__init__()` without arguments.
- `analyze.py` and `construct.py` open their text files with
  `encoding="utf-8"`.
- `tests/data/forhlr.sh` loads GROMACS 2024.3 instead of 2016.5 and requests
  GPU nodes with `#SBATCH --gres`.
- Packaging moved from `setup.py`, `MANIFEST.in` and `requirements.txt` to a
  single PEP 621 `pyproject.toml`:
  - the version is taken from `poresim.__version__`;
  - the licence is declared as the SPDX expression `GPL-3.0-or-later`;
  - the dependencies have lower bounds (`jinja2>=3.1.6`, `numpy>=1.26`,
    `pyyaml>=6.0`), jinja2 is no longer pinned;
  - test and docs tooling live in the PEP 735 dependency groups `dev` and
    `docs`, resolved with [uv](https://docs.astral.sh/uv/), instead of a
    published `dev` extra;
  - the pytest, coverage and ruff configuration is part of the file.
- Python 3.12 or newer is required; 3.12 and 3.13 are tested.
- The code base is formatted with `ruff format` and passes `ruff check` with
  the rule sets E, F, I, UP, B, C4 and RUF.
- Documentation:
  - the sources moved from `docsrc/` via `docs/` to `docs/source/` and were
    converted from reStructuredText to MyST Markdown;
  - the theme is furo, with the brand colour for links and a white plate
    behind the figures in dark mode;
  - the API reference is generated with autosummary;
  - the illustrations are SVG instead of PDF, which browsers cannot show
    inline, and the logo is no longer clipped at the top;
  - the start page reuses the feature list and installation section of the
    README.
- The README has badges for PyPI, CI, docs, CodeQL, coverage, DOI, Python and
  licence, a feature list, installation, usage and development instructions
  with uv, and the citation.
- The test suite moved from unittest to pytest, runs from `tests/` regardless
  of the working directory and does not collect the files it writes to
  `tests/output/`.

### Fixed

- `Box.add_mol` used mutable default arguments for `area`, `box` and
  `kwargs_gmx`, so values set for one molecule could leak into the next.
- A wrong import, `from re import A`, in `actuate.py`.
- A line of code that had been pasted into a comment in `simulate.py`.

### Removed

- `setup.py`, `MANIFEST.in`, `requirements.txt`, `docs/requirements.txt`,
  `docs/html.sh` and the template that forced the docs into light mode.
- The committed HTML build of the old documentation (`docs/v_old/`) and the
  `docsrc/` folder; the documentation of an older version is built from its
  release.
- The unittest suite `tests/test_simple.py`, a duplicate of the pytest suite,
  and `tests/venv.txt`.
- The workflows `workflow.yml` and `lint.yml`, replaced by `ci.yml` and
  `docs.yml`.

### Known issues

- A molecule added with `inp="fill"` and a mass but without `auto_dens` makes
  the box construction fail with a `TypeError`, in pore and in box systems,
  although `Box.add_mol` accepts it and `fill.sh` handles it.
- `Topology` keeps only the first block of a section that appears twice in a
  topology, so a second `[ dihedrals ]` block, as acpype writes for the
  impropers, is missing from the generated itp file.
- `Benchmark` with an iterator other than `"nodes"` or `"np"` fails only after
  the first files are written, with an `UnboundLocalError`.

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

[Unreleased]: https://github.com/PoreMS/PoreSim/compare/v0.3.0...HEAD
[0.3.0]: https://github.com/PoreMS/PoreSim/releases/tag/v0.3.0
