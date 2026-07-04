# v1.0.0

### Code quality
* `construct.py` — extracted `_write_positions(path, positions)` helper: writes an (N, 3) NumPy position array to a GROMACS position file, replacing repeated inline loops; used for all molecule-insertion position files
* `construct.py` — extracted `_gmx_insert_cmd(folder_gro, file_box, mol, nmol, dr, pos_file)` helper: centralises GROMACS `insert-molecules` command construction with optional `-dr` and `-ip` flags, eliminating duplicated string-building across pore/box/slit code paths
* `actuate.py` — shared folder/file path constants (`_gro/`, `_top/`, `_mdp/`, `box.gro`, `topol.top`, `index.ndx`) extracted from local variables in each method to instance attributes on `Actuate.__init__`; both `_equilibration` and `_simulation` methods now reference `self._folder_*` / `self._file_*`
* `actuate.py` — removed unused `from re import A` import
* `simulate.py` — Jinja2 template rendering for analysis scripts (`auto_dens.py` / `auto_dens_box.py`) refactored: molecule metadata collected into a single `jinja2_dict` loop; `has_fill` flag derived once instead of re-evaluated; template selection unified in a single branch
* Template scripts `auto_dens.py` and `auto_dens_box.py` trimmed; `sort.py` template removed (dead code)
* `tests/data/forhlr.sh` — GROMACS module updated from 2016.5 to 2024.3; `#SBATCH --gres` line added for GPU node support

### Tests
* Converted `tests/test_simple.py` from unittest to pytest
* Split into `test_unit.py` (fast, no file I/O) and `test_integration.py` (full output generation)
* `test_unit.py` — new unit tests: `Box` getter/setter round-trips, `add_mol` validation, `add_topol` all types, `Construct._write_positions` byte-level output, `Construct._gmx_insert_cmd` all flag combinations

### Documentation
* RST source files migrated to MyST Markdown; Sphinx theme updated to furo; API docs via sphinx-autoapi
* Docs source moved from `docsrc/` to `docs/`; previous built HTML preserved at `docs/v_old/`
* Copyright year updated to 2026; DESIGN.md added documenting the yellow/amber color palette

### CI / tooling
* GitHub Actions: added ruff linting workflow (`lint.yml`)
* GitHub Actions: added pip-audit security scan workflow (`security.yml`)
* CI matrix updated: Python 3.12–3.13; `python_requires` bumped to `>=3.12`

### Administrative
* `setup.py`: version 1.0.0, `python_requires='>=3.12'`, author email updated, jinja2 unpinned
* README: updated image paths, Python version, PyPI badge, testing and installation instructions


# v0.3.0
* New version due to a change of GitHub organisation.

# v0.2.0
* Update for PoreMS 0.3.0 and PoreAna 0.2.3
* Possibility to place molecules in a targeted manner in the pore
* For a slit pore you can set molecules near the wall to generate a layer on the wall
* Insertion of the molecules with Gromacs position files

# v0.1.2
* Filling of pores also possible for mixtures (target density can be defined for more than one molecule)

# v0.1.1
* Update for PoreMS 0.2.5 and PoreAna 0.2.2
* Improved templates
* Update simulation docs

# v0.1.0
* Initial version
