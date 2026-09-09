# Build the SCS paper. Needs pdflatex (any TeX Live / MacTeX install).
#
# Three passes: the first writes the .aux and .toc, the second resolves the
# cross-references, and the third settles them again in case the table of
# contents changed length and moved everything by a page. There is no .bib --
# the bibliography is a `thebibliography` environment inside scs.tex -- so
# bibtex is not needed.

LATEX ?= pdflatex
LATEXFLAGS ?= -interaction=nonstopmode -halt-on-error

.PHONY: all clean distclean

all: scs.pdf

scs.pdf: scs.tex
	$(LATEX) $(LATEXFLAGS) scs.tex
	$(LATEX) $(LATEXFLAGS) scs.tex
	$(LATEX) $(LATEXFLAGS) scs.tex

clean:
	rm -f scs.aux scs.log scs.out scs.toc

distclean: clean
	rm -f scs.pdf
