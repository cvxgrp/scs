# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line, and also
# from the environment for the first two.
SPHINXOPTS    ?=
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
PYTHON    ?= python3

.PHONY: example_outputs check_example_deps

check_example_deps:
	@$(PYTHON) -c "import scs, numpy, scipy, cvxpy" 2>/dev/null || { \
		echo "error: building the docs runs the examples, which need scs."; \
		echo "       install the dependencies with:"; \
		echo "           pip install scs numpy scipy cvxpy"; \
		exit 1; \
	}

example_outputs: check_example_deps
	@echo "Generating example output..."
	@set -e; for f in $(EXAMPLEDIR)/python/*.py; do \
		echo "  $(PYTHON) $${f##*/}"; \
		( cd $(EXAMPLEDIR)/python && $(PYTHON) "$${f##*/}" > "$${f##*/}.out" ); \
	done
	@echo "  cc     qp.c"
	@$(MAKE) -C $(SCSROOT) --no-print-directory out/libscsdir.a
	@$(CC) -I$(SCSROOT)/include -o $(EXAMPLEDIR)/qp.out $(EXAMPLEDIR)/qp.c \
		$(SCSROOT)/out/libscsdir.a -lm -llapack -lblas
	@( cd $(EXAMPLEDIR) && ./qp.out > qp.c.out )

docs: example_outputs
	@make html
	@cp -a _build/html/. ..
	@touch ../.nojekyll
	@rm -rf doxygen_out/*
	@make clean

