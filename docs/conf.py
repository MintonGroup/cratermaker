import datetime
import inspect
import os
import re
import sys
from pathlib import Path

import pyvista
from pyvista.plotting.utilities.sphinx_gallery import DynamicScraper

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

# necessary when building the sphinx gallery
pyvista.BUILDING_GALLERY = True
pyvista.OFF_SCREEN = True
os.environ["PYVISTA_BUILDING_GALLERY"] = "true"

sys.path.insert(0, os.path.abspath("_exts"))
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.intersphinx",
    "sphinx.ext.extlinks",
    "sphinx.ext.mathjax",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.doctest",
    "IPython.sphinxext.ipython_directive",
    "IPython.sphinxext.ipython_console_highlighting",
    "sphinx_autosummary_accessors",
    "sphinx_copybutton",
    "sphinx_design",
    "sphinx_inline_tabs",
    "sphinx_gallery.gen_gallery",
    "cratermaker_autodoc",
    "pyvista.ext.plot_directive",
    "pyvista.ext.viewer_directive",
    "sphinxcontrib.video",
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
    "member-order": "bysource",
    "no-index-entry": True,
}
autoapi_generate_api_docs = False

intersphinx_mapping = {
    "python": ("https://docs.python.org/3/", None),
    "numpy": ("https://numpy.org/doc/stable", None),
    "xarray": ("https://docs.xarray.dev/en/stable/", None),
    "uxarray": ("https://uxarray.readthedocs.io/en/latest/", None),
}

templates_path = ["_templates"]

sphinx_gallery_conf = {
    "examples_dirs": "../examples",
    "gallery_dirs": "auto_examples",
    "filename_pattern": r"\.py",  # Allows all Python files to be found, not just ones that begin with `plot_`
    "image_scrapers": (DynamicScraper(), "matplotlib"),
    "parallel": True,
}
html_theme = "sphinx_book_theme"
html_title = ""

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_css_files = ["style.css"]
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
html_theme_options = {
    "repository_url": "https://github.com/MintonGroup/cratermaker",
    "repository_branch": "main",
    "path_to_docs": "docs",
    "use_edit_page_button": True,
    "use_repository_button": True,
    "use_issues_button": True,
    "home_page_in_toc": False,
    "extra_footer": """<p>Development of Cratermaker was supported by NASA Lunar Data Analysis Program Grants <a href="https://www.usaspending.gov/award/ASST_NON_80NSSC21K1719_8000">#80NSSC21K1719</a> and <a href="https://www.usaspending.gov/award/ASST_NON_80NSSC25K7050_8000">#80NSSC25K7050</a><br>
    Theme by the <a href="https://ebp.jupyterbook.org">Executable Book Project</a></p>""",
    "logo": {
        "image_light": "_images/logos/Cratermaker_Social_Preview_light.svg",
        "image_dark": "_images/logos/Cratermaker_Social_Preview_dark.svg",
    },
    "show_toc_level": 4,
}


def html_page_context(app, pagename, templatename, context, doctree):
    # Disable edit button for docstring generated pages
    if "generated" in pagename:
        context["theme_use_edit_page_button"] = False


rst_prolog = """
.. |simdir| replace:: The main project simulation directory. Default is the current working directory if None.
.. |rng| replace:: A numpy random number generator. If None, a new generator is created using the rng_seed if it is provided.
.. |rng_seed| replace:: The rng_seed for the RNG. If None, a new RNG is created.
.. |rng_state| replace:: The state of the random number generator. If None, a new state is created.
.. |kwargs| replace:: Additional keyword arguments that are either ignored or passed to internal functions as needed.
.. |interval_export| replace:: The interval number to export. If None, all intervals currently saved will be exported. Default is None.
.. |ask_overwrite_methods| replace:: If True, the user will be prompted to confirm before overwriting any existing files. If False, existing files will be overwritten without confirmation. If None, the default behavior of the class will be used. This will only persist for the duration of the export, and will be reset to its original value afterwards.
.. |Simulation| replace:: :py:class:`~cratermaker.core.simulation.Simulation`
.. |sim.run| replace:: :py:class:`Simulation.run() <cratermaker.core.simulation.Simulation.run>`
.. |sim.populate| replace:: :py:class:`Simulation.populate() <cratermaker.core.simulation.Simulation.populate>`
.. |sim.emplace| replace:: :py:class:`Simulation.emplace() <cratermaker.core.simulation.Simulation.emplace>`
.. |sim.smallest_crater| replace:: :py:attr:`~cratermaker.core.simulation.Simulation.smallest_crater`
.. |Production| replace:: :py:class:`~cratermaker.components.production.Production`
.. |NPF| replace:: :py:class:`~cratermaker.components.production.neukum.NeukumProduction`
.. |PowerLawProduction| replace:: :py:class:`~cratermaker.components.production.powerlaw.PowerLawProduction`
.. |production.function| replace:: :py:meth:`Production.function() <cratermaker.components.production.Production.function>`
.. |production.sample| replace:: :py:meth:`Production.sample() <cratermaker.components.production.Production.sample>`
.. |production.age_from_D_N| replace:: :py:meth:`Production.age_from_D_N() <cratermaker.components.production.Production.age_from_D_N>`
.. |production.N_D_units| replace:: :py:attr:`~cratermaker.components.production.Production.N_D_units`
.. |production.N_conversion_factor| replace:: :py:attr:`~cratermaker.components.production.Production.N_conversion_factor`
.. |production.D_conversion_factor| replace:: :py:attr:`~cratermaker.components.production.Production.D_conversion_factor`
.. |Surface| replace:: :py:class:`~cratermaker.components.surface.Surface`
.. |LocalSurface| replace:: :py:class:`~cratermaker.components.surface.LocalSurface`
.. |Target| replace:: :py:class:`~cratermaker.components.target.Target`
.. |Crater| replace:: :py:class:`~cratermaker.components.crater.Crater`
.. |Crater.maker| replace:: :py:meth:`Crater.maker() <cratermaker.components.crater.Crater.maker>`
.. |Crater.diameter| replace:: :py:attr:`~cratermaker.components.crater.CraterFixed.diameter`
.. |Crater.projectile_diameter| replace:: :py:attr:`~cratermaker.components.crater.CraterFixed.projectile_diameter`
.. |Crater.production_time| replace:: :py:attr:`~cratermaker.components.crater.CraterVariable.production_time`
.. |Crater.production_ND| replace:: :py:attr:`~cratermaker.components.crater.CraterVariable.production_ND`
.. |Crater.production_sequence| replace:: :py:attr:`~cratermaker.components.crater.CraterVariable.production_sequence`
.. |Morphology| replace:: :py:class:`~cratermaker.components.morphology.Morphology`
.. |BasicMoon| replace:: :py:class:`~cratermaker.components.morphology.basicmoon.BasicMoonMorphology`
.. |Counting| replace:: :py:class:`~cratermaker.components.counting.Counting`
.. |Projectile| replace:: :py:class:`~cratermaker.components.projectile.Projectile`
.. |Scaling| replace:: :py:class:`~cratermaker.components.scaling.Scaling`
"""
