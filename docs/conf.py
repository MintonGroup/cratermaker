import datetime
import os
import sys
from contextlib import suppress

import sphinx_autosummary_accessors
import yaml
from sphinx.application import Sphinx
from sphinx.util import logging
import os

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information
project = 'Cratermaker'
copyright = f'{datetime.datetime.now().year}, David A. Minton'
author = 'David A. Minton'
with open(os.path.join(os.path.dirname(__file__), os.pardir, "version.txt"), 'r') as file:
    version = file.read().strip() 
release = version

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.intersphinx",
    "sphinx.ext.mathjax",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
]


templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store','**.ipynb_checkpoints']


# The name of the Pygments (syntax highlighting) style to use.
pygments_style = "sphinx"


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output


# Use Autodoc and Napolean for extracting docstrings
autosummary_generate = True
autodoc_typehints = "none"
napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_use_param = False
napoleon_use_rtype = False

intersphinx_mapping = {
    "python": ("https://docs.python.org/3/", None),
    "numpy": ("https://numpy.org/doc/stable", None),
    "xarray" : ("https://docs.xarray.dev/en/stable/", None),
    "uxarray" : ("https://uxarray.readthedocs.io/", None),
}

templates_path = ["_templates"]

html_theme = 'sphinx_book_theme'
html_title =""
html_static_path = ["_static"]

html_context = {
    "github_user": "profminton",
    "github_repo": "cratermaker",
    "github_version": "main",
    "doc_path": "docs",
}

# Theme options are theme-specific and customize the look and feel of a theme
# further.  For a list of options available for each theme, see the
# documentation.
html_theme_options = dict(
    # analytics_id=''  this is configured in rtfd.io
    # canonical_url="",
    repository_url="https://github.com/profminton/cratermaker",
    repository_branch="main",
    navigation_with_keys=False,  # pydata/pydata-sphinx-theme#1492
    path_to_docs="docs",
    use_edit_page_button=True,
    use_repository_button=True,
    use_issues_button=True,
    home_page_in_toc=False,
    extra_footer="""<p>Development of Cratermaker is supported by NASA Lunar Data Analysis Program Grant #80NSSC21K1719<br>
    Theme by the <a href="https://ebp.jupyterbook.org">Executable Book Project</a></p>""",
    icon_links=[],  # workaround for pydata/pydata-sphinx-theme#1220
    announcement="🍾 <a href='https://github.com/profminton/cratermaker/discussions/1'>Cratermaker is currently under development</a> 🎉",
)


# The name of an image file (relative to this directory) to place at the top
# of the sidebar.
html_logo = "_static/logos/Cratermaker_Social_Preview.svg"

# The name of an image file (within the static path) to use as favicon of the
# docs.  This file should be a Windows icon file (.ico) being 16x16 or 32x32
# pixels large.
html_favicon = "_static/logos/Cratermaker_Icon.svg"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]
html_css_files = ["style.css"]


