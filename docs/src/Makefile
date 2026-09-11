# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line, and also
# from the environment for the first two.
# -W turns warnings into errors, so a broken cross reference or a missing
# include fails the build instead of shipping to the site. Override on the
# command line (SPHINXOPTS=) to build a work in progress that still warns.
SPHINXOPTS    ?= -W
SPHINXBUILD   ?= sphinx-build
SOURCEDIR     = .
BUILDDIR      = _build

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

# --- Generated example output ------------------------------------------------
#
# The captured output shown on the example pages is generated, not committed,
# so it can never fall out of step with the solver (see issue #448). The
# Python examples need `scs`, `numpy` and `scipy`, plus `cvxpy` (entropy.py
# cross-checks its answer against CVXPY); the C example is built against this
# tree, so it needs a compiler and BLAS/LAPACK.
#
# MATLAB is deliberately excluded: regenerating `examples/qp.m.out` needs a
# MATLAB licence, which CI does not have, so that file stays committed and is
# refreshed by hand. See examples/matlab.rst.
#
# Sphinx has no optional `literalinclude`, so a missing .out file is a hard
# build error. `docs` therefore depends on this target, and the dependency
# check below fails with an actionable message rather than letting Sphinx
# report a confusing missing-file error.

EXAMPLEDIR = examples
SCSROOT    = ../..
# Named DOCS_PYTHON rather than PYTHON: a command-line PYTHON=... would be
# inherited by the sub-make of the C library below, where PYTHON means
# "build the Python extension" and adds -DPYTHON to the C flags.
DOCS_PYTHON ?= python3

.PHONY: example_outputs check_example_deps guide

check_example_deps:
	@$(DOCS_PYTHON) -c "import scs, numpy, scipy, cvxpy" 2>/dev/null || { \
		echo "error: building the docs runs the examples, which need scs."; \
		echo "       install the dependencies with:"; \
		echo "           pip install scs numpy scipy cvxpy"; \
		exit 1; \
	}

example_outputs: check_example_deps
	@echo "Generating example output..."
	@set -e; for f in $(EXAMPLEDIR)/python/*.py; do \
		echo "  $(DOCS_PYTHON) $${f##*/}"; \
		( cd $(EXAMPLEDIR)/python && $(DOCS_PYTHON) "$${f##*/}" > "$${f##*/}.out" ); \
	done
	@echo "  cc     qp.c"
	@$(MAKE) -C $(SCSROOT) --no-print-directory out/libscsdir.a
	@$(CC) -I$(SCSROOT)/include -o $(EXAMPLEDIR)/qp.out $(EXAMPLEDIR)/qp.c \
		$(SCSROOT)/out/libscsdir.a -lm -llapack -lblas
	@( cd $(EXAMPLEDIR) && ./qp.out > qp.c.out )

# --- PDF user guide ---------------------------------------------------------
#
# The same sources as the site, rendered through Sphinx's LaTeX builder, so the
# guide cannot drift from the site; the release number in conf.py appears on
# the title page. pdflatex is run directly (three passes, as for any LaTeX
# document with a table of contents) rather than through latexmk, so the only
# requirement beyond the docs toolchain is a TeX installation.
guide: example_outputs
	@$(SPHINXBUILD) -b latex "$(SOURCEDIR)" "$(BUILDDIR)/latex" $(SPHINXOPTS) $(O)
	@cd "$(BUILDDIR)/latex" && for pass in 1 2 3; do \
		pdflatex -interaction=nonstopmode -halt-on-error scs_user_guide.tex > pdflatex.out 2>&1 || \
			{ tail -40 pdflatex.out; exit 1; }; \
	done
	@echo "PDF user guide written to $(BUILDDIR)/latex/scs_user_guide.pdf"

docs: example_outputs
	@make html
	@cp -a _build/html/. ..
	@touch ../.nojekyll
	@rm -rf doxygen_out/*
	@make clean

