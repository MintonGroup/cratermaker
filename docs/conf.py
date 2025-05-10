import datetime
import inspect
import os
import sys

import cratermaker

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information
project = "Cratermaker"
copyright = f"{datetime.datetime.now().year}, David A. Minton"
author = "David A. Minton"
version = cratermaker.__version__
release = version

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

sys.path.insert(0, os.path.abspath("_exts"))
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.intersphinx",
    "sphinx.ext.extlinks",
    "sphinx.ext.mathjax",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "IPython.sphinxext.ipython_directive",
    "IPython.sphinxext.ipython_console_highlighting",
    "sphinx_autosummary_accessors",
    "sphinx.ext.linkcode",
    "sphinxext.opengraph",
    "sphinx_copybutton",
    "sphinx_design",
    "sphinx_inline_tabs",
    "sphinx_gallery.gen_gallery",
    "cratermaker_autodoc",
]

extlinks = {
    "issue": ("https://github.com/MintonGroup/cratermaker/issues/%s", "GH%s"),
    "pull": ("https://github.com/MintonGroup/cratermaker/pull/%s", "PR%s"),
    "discussion": ("https://github.com/MintonGroup/cratermaker/discussions/%s", "D%s"),
    "release": ("https://github.com/MintonGroup/cratermaker/releases/tag/%s", "%s"),
}


templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store", "**.ipynb_checkpoints"]


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

autodoc_default_options = {
    "members": True,
    "undoc-members": False,
    "show-inheritance": True,
}

intersphinx_mapping = {
    "python": ("https://docs.python.org/3/", None),
    "numpy": ("https://numpy.org/doc/stable", None),
    "xarray": ("https://docs.xarray.dev/en/stable/", None),
}

templates_path = ["_templates"]

sphinx_gallery_conf = {
    "examples_dirs": "../examples",
    "gallery_dirs": "auto_examples",
}

html_theme = "sphinx_book_theme"
html_title = ""
html_static_path = ["_static"]

html_context = {
    "github_user": "MintonGroup",
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
    repository_url="https://github.com/MintonGroup/cratermaker",
    repository_branch="main",
    navigation_with_keys=False,  # pydata/pydata-sphinx-theme#1492
    path_to_docs="docs",
    use_edit_page_button=True,
    use_repository_button=True,
    use_issues_button=True,
    home_page_in_toc=False,
    extra_footer="""<p>Development of Cratermaker was supported by NASA Lunar Data Analysis Program Grants <a href="https://www.usaspending.gov/award/ASST_NON_80NSSC21K1719_8000">#80NSSC21K1719</a> and <a href="https://www.usaspending.gov/award/ASST_NON_80NSSC25K7050_8000">#80NSSC25K7050</a><br>
    Theme by the <a href="https://ebp.jupyterbook.org">Executable Book Project</a></p>""",
    icon_links=[],  # workaround for pydata/pydata-sphinx-theme#1220
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

# configuration for sphinxext.opengraph
ogp_site_url = "https://cratermaker.readthedocs.io/en/latest/"
ogp_image = "https://cratermaker.readthedocs.io/en/stable/_static/logos/Cratermaker_Social_Preview.png"
ogp_custom_meta_tags = [
    '<meta name="image" property="og:image" content="https://cratermaker.readthedocs.io/en/stable/_static/logos/Cratermaker_Social_Preview.png" />',
]


# based on numpy doc/source/conf.py
def linkcode_resolve(domain, info):
    """
    Determine the URL corresponding to Python object
    """
    if domain != "py":
        return None

    modname = info["module"]
    fullname = info["fullname"]

    submod = sys.modules.get(modname)
    if submod is None:
        return None

    obj = submod
    for part in fullname.split("."):
        try:
            obj = getattr(obj, part)
        except AttributeError:
            return None

    try:
        fn = inspect.getsourcefile(inspect.unwrap(obj))
    except TypeError:
        fn = None
    if not fn:
        return None

    try:
        source, lineno = inspect.getsourcelines(obj)
    except OSError:
        lineno = None

    if lineno:
        linespec = f"#L{lineno}-L{lineno + len(source) - 1}"
    else:
        linespec = ""

    fn = os.path.relpath(fn, start=os.path.dirname(cratermaker.__file__))


def html_page_context(app, pagename, templatename, context, doctree):
    # Disable edit button for docstring generated pages
    if "generated" in pagename:
        context["theme_use_edit_page_button"] = False
