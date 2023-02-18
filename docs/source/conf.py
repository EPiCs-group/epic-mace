# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import sys, os
sys.path.insert(0, os.path.abspath('../../'))

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'epic-mace'
copyright = '2023, Ivan Yu. Chernyshov'
author = 'Ivan Yu. Chernyshov'
release = '0.4.2'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'nbsphinx',
    'sphinx_gallery.load_style',
    'myst_parser'
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']

nbsphinx_thumbnails = {
    'notebooks/general_pipeline.ipynb': 'images/mace.png',
    'notebooks/complex_init.ipynb': 'images/mace.png',
    'notebooks/stereomers.ipynb': 'images/mace.png',
    'notebooks/embedding.ipynb': 'images/mace.png',
    'notebooks/features.ipynb': 'images/mace.png'
}


html_theme_options = {
    'logo_only': False,
}
html_logo = 'images/mace.png'
html_favicon = 'images/mace.png'


# -- Options for exts ----------------------------------------------

# myst-nb
nb_execution_mode = "off"
nb_execution_cache_path = "source/notebooks/cache"

