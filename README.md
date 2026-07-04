<img src="https://github.com/PoreMS/PoreSim/blob/main/docs/pics/logo_text_sub.svg" width="60%">

--------------------------------------

[![PyPI Version](https://img.shields.io/badge/PyPI-1.0.0-orange)](https://pypi.org/project/poresim/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://github.com/PoreMS/PoreSim/blob/main/LICENSE)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17395962.svg)](https://doi.org/10.5281/zenodo.17395962)
[![Build Status](https://github.com/PoreMS/PoreSim/actions/workflows/workflow.yml/badge.svg)](https://github.com/PoreMS/PoreSim/actions/workflows/workflow.yml)
[![codecov](https://codecov.io/gh/PoreMS/PoreSim/branch/main/graph/badge.svg)](https://codecov.io/gh/PoreMS/PoreSim)

## Documentation

Online documentation is available at [porems.github.io/PoreSim](https://porems.github.io/PoreSim/).

The docs include an example for creating [simulation boxes](https://porems.github.io/PoreSim/simulation.html) and an [API reference](https://porems.github.io/PoreSim/autoapi/index.html).


## Dependencies

PoreSim requires Python 3.12+.

Installation requires [pyyaml](https://pypi.org/project/PyYAML/), [numpy](https://numpy.org/), and [jinja2](https://pypi.org/project/Jinja2/).


## Installation

The latest stable release can be installed from PyPI:

    pip install poresim

Or install the development version directly from GitHub:

    pip install git+https://github.com/PoreMS/PoreSim.git#egg=poresim

Or download the repository and install in the top directory via:

    pip install .


## Testing

Install in editable mode with test dependencies:

    pip install -e ".[dev]"

Then run the tests:

    pytest tests/test_unit.py          # fast unit tests
    pytest tests/test_integration.py   # full integration tests (slow)


## Development

PoreSim development takes place on GitHub: [www.github.com/PoreMS/PoreSim](https://github.com/PoreMS/PoreSim)

Please submit any reproducible bugs you encounter to the [issue tracker](https://github.com/PoreMS/PoreSim/issues).


## How to Cite PoreSim

When citing PoreSim please use the current **Zenodo DOI** corresponding to the used PoreSim version. (Current DOI is listed in the badges.)

## Legacy Code Notice

This repository contains the actively maintained and current codebase.
Earlier development stages and archived versions are stored in the following legacy repository:

[https://github.com/Ajax23/PoreSim](https://github.com/Ajax23/PoreSim)

Please note:
- Legacy versions are no longer actively maintained
- APIs, file formats, and dependencies may differ
