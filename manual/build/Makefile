# Makefile for the PBeFS manual (single source -> web + PDF).
# Run from manual/build/.
#
#   make pdf       # PDF via Pandoc + pdflatex       -> ../PBeFS-manual.pdf
#   make html      # static MkDocs site              -> ../../manual-site/
#   make serve     # live MkDocs preview at 127.0.0.1:8000
#   make figures   # (re)build ../figures/NAME.pdf and .svg from NAME.tex
#   make generated # the excerpts of code (snarf.py) and the tables
#                  # (tables.py, from ../results if it is there)
#   make clean     # remove build artifacts
#
# Prerequisites:
#   pandoc, a LaTeX distribution (pdflatex + the usual packages),
#   python3 with mkdocs + mkdocs-material + pymdown-extensions,
#   poppler-utils (pdftocairo) for figure conversion.

MANUAL  := ..
FIGDIR  := $(MANUAL)/figures
PDF_OUT := $(MANUAL)/PBeFS-manual.pdf
PANDOC  ?= pandoc
MKDOCS  ?= mkdocs
FIGTEX  := $(filter-out $(FIGDIR)/pbefs-tikz.tex $(FIGDIR)/tikz-pbefs.tex,\
                        $(wildcard $(FIGDIR)/*.tex))
# lamalunga, which computes the palette of the figures from one primary:
# the copy of its author if there is one, else the one in this directory
LAMALUNGA ?= $(firstword $(wildcard $(HOME)/GDrive/lamalunga) $(CURDIR))

.PHONY: pdf html serve figures generated clean

# the .pdf of a figure is a product, not an intermediate file of its .svg
.PRECIOUS: $(FIGDIR)/%.pdf

pdf: figures generated
	$(PANDOC) --defaults=pandoc.yaml -o $(PDF_OUT)
	@echo "PDF written to $(PDF_OUT)"

html: generated
	$(MKDOCS) build --strict -f mkdocs.yml
	@echo "Static site written to ../../manual-site/"

serve: generated
	$(MKDOCS) serve -f mkdocs.yml

generated:
	python3 snarf.py
	python3 tables.py

# Each TikZ source NAME.tex (a standalone document that inputs the shared
# styles of pbefs-tikz.tex) gives NAME.pdf for the PDF and NAME.svg for the
# web; a figure shipped only as a .pdf gets its .svg as well.
figures: $(FIGTEX:.tex=.svg)
	@for pdf in $(FIGDIR)/*.pdf; do \
	  [ -e "$$pdf" ] || continue; \
	  svg="$${pdf%.pdf}.svg"; \
	  if [ ! -f "$$svg" ]; then pdftocairo -svg "$$pdf" "$$svg"; fi; \
	done

$(FIGDIR)/%.pdf: $(FIGDIR)/%.tex $(FIGDIR)/tikz-pbefs.tex
	cd $(FIGDIR) && TEXINPUTS=$(LAMALUNGA)//: \
	    pdflatex -interaction=nonstopmode -halt-on-error \
	    $*.tex > /dev/null && rm -f $*.aux $*.log

$(FIGDIR)/%.svg: $(FIGDIR)/%.pdf
	pdftocairo -svg $< $@

clean:
	rm -rf ../../manual-site ../generated
	rm -f $(PDF_OUT)
