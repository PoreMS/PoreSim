import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from poresim import __version__

project = "PoreSim"
copyright = "2021-2026, Hamzeh Kraus, Marc Högler"
author = "Hamzeh Kraus, Marc Högler"
release = __version__
version = __version__

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.intersphinx",
    "sphinx.ext.mathjax",
    "sphinx_autodoc_typehints",
    "sphinx.ext.autosummary",
    "myst_parser",
]

source_suffix = {".rst": "restructuredtext", ".md": "markdown"}

myst_enable_extensions = ["colon_fence"]

autosummary_generate = True

autodoc_default_options = {
    "members": True,
    "undoc-members": True,
    "show-inheritance": True,
}
autodoc_member_order = "bysource"

napoleon_google_docstring = False
napoleon_numpy_docstring = True

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable", None),
}

exclude_patterns = []

html_theme = "furo"
html_static_path = ["_static"]
html_css_files = ["style.css"]
html_logo = "_static/logo_text.svg"
html_favicon = "_static/favicon.ico"
