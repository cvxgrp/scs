# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))

import subprocess
import sphinx_rtd_theme

# -- Project information -----------------------------------------------------

project = "SCS"
copyright = "2021, Brendan O'Donoghue"
author = "Brendan O'Donoghue"

# The full version, including alpha/beta/rc tags
__version__ = "3.2.9"

release = __version__
version = __version__

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = ["sphinx.ext.mathjax", "breathe", "sphinx_rtd_theme"]

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# sphinx pygments style uses ugly green boxes for code blocks
# pygments_style = 'sphinx'
pygments_style = "default"

# html_sidebars = {
#   '**': [
#       'about.html', 'navigation.html', 'searchbox.html',
#   ]
# }

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.

# html_theme = 'alabaster'
html_theme = "sphinx_rtd_theme"
html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]


def setup(app):
    app.add_css_file("css/scs_theme.css")


html_logo = "_static/scs_logo_transparent.png"
html_favicon = "_static/favicon.ico"
html_theme_options = {
    "logo_only": True,
    "display_version": True,
    #'github_banner': True,
    #'github_user': 'cvxgrp',
    #'github_repo': 'scs',
    #'logo': 'scs_logo_transparent.png',
    #'logo_name': False,
    #'github_button': False,
    #'github_type': 'star',
    "analytics_id": "G-9CY7R8S5N2",
}

rst_epilog = ".. |version| replace:: %s" % __version__

# Breathe docs
subprocess.call("doxygen Doxyfile", shell=True)

breathe_projects = {"scs": "doxygen_out/xml/"}
breathe_default_project = "scs"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]
