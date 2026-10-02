"""Sphinx configuration for the nctpy documentation.

Matches the documentation setup of the lab's other packages (snaplab_tools): pydata-sphinx-theme, MyST Markdown
and myst-nb, NumPy-style docstrings through autosummary. Notebooks are committed with their outputs and only
rendered here; CI is what executes them.
"""

import re
from pathlib import Path

_ROOT = Path(__file__).resolve().parents[2]

# -- Project information -----------------------------------------------------------------------
project = "nctpy"
author = "Linden Parkes, Jason Z. Kim, Jennifer Stiso"
copyright = "2022-%Y, Linden Parkes, Jason Z. Kim, Jennifer Stiso"

# Read the version from the package source rather than importing it, so the docs build does not depend on the
# package being importable at config time.
_init = (_ROOT / "src" / "nctpy" / "__init__.py").read_text()
release = re.search(r'^__version__ = ["\'](.*)["\']', _init, re.M).group(1)
version = release

# -- General configuration ---------------------------------------------------------------------
extensions = [
    "myst_nb",  # MyST Markdown + notebooks (pulls in myst_parser)
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",  # NumPy-style docstrings, as used throughout the package
    "sphinx.ext.intersphinx",
    "sphinx.ext.mathjax",
    "sphinx.ext.viewcode",
    "sphinx_copybutton",
    "sphinx_design",
]
templates_path = ["_templates"]
exclude_patterns = [
    "_build",
    "**.ipynb_checkpoints",
    "Thumbs.db",
    ".DS_Store",
    # The example notebooks are rendered through their .rst copies until each page is converted.
    "pages/examples/*.ipynb",
]

# -- Autodoc / autosummary ---------------------------------------------------------------------
autosummary_generate = True
autodoc_default_options = {
    "members": True,
    "inherited-members": False,
}
autodoc_typehints = "description"
# Types come from the docstrings; annotations only add a type to a documented parameter that has none.
autodoc_typehints_description_target = "documented"
autodoc_member_order = "bysource"
napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_use_param = True
napoleon_use_rtype = True
napoleon_preprocess_types = True
# The optional dependencies of nctpy.plotting ([plot]) and nctpy.optimize ([optimize]) are not part of the docs
# install, so they are mocked: the API reference is built from the docstrings and never runs this code.
autodoc_mock_imports = ["matplotlib", "nibabel", "nilearn", "seaborn", "torch"]

# -- Intersphinx -------------------------------------------------------------------------------
intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "scipy": ("https://docs.scipy.org/doc/scipy/", None),
    "matplotlib": ("https://matplotlib.org/stable/", None),
    "nilearn": ("https://nilearn.github.io/stable/", None),
}
nitpicky = False  # missing objects in third-party packages should not fail a -W build

# -- MyST / notebooks --------------------------------------------------------------------------
myst_enable_extensions = [
    "colon_fence",
    "deflist",
    "dollarmath",
    "substitution",
]
myst_heading_anchors = 3

# Notebooks are NOT executed here: they are committed with their outputs and this build only renders them, so the
# published documentation never depends on the scientific stack or on downloads succeeding. CI executes them.
nb_execution_mode = "off"
nb_merge_streams = True

# -- HTML output -------------------------------------------------------------------------------
html_theme = "pydata_sphinx_theme"
html_static_path = ["_static"]
html_css_files = ["custom.css"]
html_title = f"{project} {release}"
html_theme_options = {
    "github_url": "https://github.com/LindenParkesLab/nctpy",
    "icon_links": [
        {
            "name": "GitHub",
            "url": "https://github.com/LindenParkesLab/nctpy",
            "icon": "fa-brands fa-github",
        },
    ],
    "navbar_end": ["theme-switcher", "navbar-icon-links"],
    "show_prev_next": True,
    "navigation_with_keys": False,
    "show_toc_level": 2,
    "footer_start": ["copyright"],
    "footer_end": ["sphinx-version"],
}
html_sidebars = {
    "index": [],
}
html_context = {
    "github_user": "LindenParkesLab",
    "github_repo": "nctpy",
    "github_version": "main",
    "doc_path": "docs/source",
    "default_mode": "auto",
}
