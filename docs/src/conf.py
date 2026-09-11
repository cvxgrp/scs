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

# -- Project information -----------------------------------------------------

project = "SCS"
copyright = "2021, Brendan O'Donoghue"
author = "Brendan O'Donoghue"

# The full version, including alpha/beta/rc tags
__version__ = "3.3.1"

release = __version__
version = __version__

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = ["sphinx.ext.mathjax", "breathe", "sphinx_rtd_theme"]

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

html_theme = "sphinx_rtd_theme"


def setup(app):
    app.add_css_file("css/scs_theme.css")


html_logo = "_static/scs_logo_transparent.png"
html_favicon = "_static/favicon.ico"
html_theme_options = {
    "logo_only": True,
    #'github_banner': True,
    #'github_user': 'cvxgrp',
    #'github_repo': 'scs',
    #'logo': 'scs_logo_transparent.png',
    #'logo_name': False,
    #'github_button': False,
    #'github_type': 'star',
}

# Google Analytics (GA4). Previously configured with the sphinx_rtd_theme
# "analytics_id" theme option, which the theme deprecated.
html_js_files = [
    (
        "https://www.googletagmanager.com/gtag/js?id=G-9CY7R8S5N2",
        {"async": "async"},
    ),
    "js/analytics.js",
]

# Breathe docs
subprocess.call("doxygen Doxyfile", shell=True)

breathe_projects = {"scs": "doxygen_out/xml/"}
breathe_default_project = "scs"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]

# -- Options for the PDF user guide (LaTeX builder) ---------------------------
#
# `make guide` builds docs/src/_build/latex/scs_user_guide.pdf from the same
# sources as the HTML site, so the two cannot drift; the release number above
# appears on the title page, which is what versions the guide. Its master
# document, guide/index.rst, selects the pages the PDF contains.
latex_engine = "pdflatex"
latex_documents = [
    ("guide/index", "scs_user_guide.tex", "SCS User Guide", author, "manual"),
]
latex_logo = "_static/scs_logo_transparent.png"
latex_show_urls = "footnote"
latex_elements = {
    "papersize": "letterpaper",
    "pointsize": "11pt",
    # the site's front page repeats the title; the PDF gets a table of contents
    "tableofcontents": r"\sphinxtableofcontents",
    # pdflatex cannot typeset these characters on its own; they appear in code
    # comments (the JavaScript example, the aa.h docstrings) that reach the
    # PDF through literalinclude and Breathe.
    "preamble": r"""
\DeclareUnicodeCharacter{2208}{\ensuremath{\in}}
\DeclareUnicodeCharacter{2264}{\ensuremath{\leq}}
\DeclareUnicodeCharacter{2265}{\ensuremath{\geq}}
\DeclareUnicodeCharacter{2192}{\ensuremath{\rightarrow}}
\DeclareUnicodeCharacter{00D7}{\ensuremath{\times}}
\DeclareUnicodeCharacter{03B3}{\ensuremath{\gamma}}
% Sphinx typesets every struct member and function as its own list
% environment with generous vertical padding, which stretches the C API
% reference over pages. Tighten it: no extra space above or below each
% entry, and none between a signature and its one-line description.
\makeatletter
\renewenvironment{fulllineitems}{%
  \begin{list}{}{\labelwidth \leftmargin
                 \rightmargin \z@ \topsep \z@ \partopsep \z@
                 \itemsep \z@ \parsep \z@
                 \let\makelabel=\py@itemnewline}%
}{\end{list}}
\makeatother
% Tables in a manual read better a size down, and the settings and compile
% flag tables have long monospace names that do not fit a column at 11pt.
\usepackage{etoolbox}
\AtBeginEnvironment{longtable}{\small}
\AtBeginEnvironment{tabulary}{\small}
\AtBeginEnvironment{tabular}{\small}
""",
}
