# --------------------------------------------------------------------------- #
#    Makefile of the manual of PBeFS                                          #
#                                                                             #
#        make              pbefs-manual.pdf, by latexmk and biber             #
#        make chapter CH=<file of chapters/, without .tex> [OUT=<dir>]        #
#                          that chapter alone, OUT/chapter-<CH>.pdf           #
#        make html         the pages of HTML in html/, by make4ht             #
#        make clean        everything that the targets above write            #
#                                                                             #
#    Each target first writes generated/ by snarf.py (the excerpts of the     #
#    headers and the version) and by tables.py (the tables of the             #
#    benchmarks, from results/ if it is there).                               #
#                                                                             #
#                                Donato Meoli                                 #
#                         Dipartimento di Informatica                         #
#                             Università di Pisa                              #
# --------------------------------------------------------------------------- #

LATEXMK = latexmk -pdf -bibtex -interaction=nonstopmode -halt-on-error
CH ?= introduction
OUT ?= build-$(CH)

all: pbefs-manual.pdf

generated: FORCE
	python3 snarf.py
	python3 tables.py

pbefs-manual.pdf: generated FORCE
	mkdir -p build/chapters
	$(LATEXMK) -outdir=build pbefs-manual.tex
	cp build/pbefs-manual.pdf .

chapter: generated
	mkdir -p $(OUT)/chapters
	$(LATEXMK) -jobname=chapter-$(CH) -outdir=$(OUT) \
	    -usepretex='\def\onlychapter{$(CH)}' pbefs-manual.tex

html: generated
	mkdir -p html
	make4ht -c pbefs-manual.cfg -f html5+common_domfilters -d html \
	    -l -B build-html -e pbefs-manual.mk4 pbefs-manual.tex "fn-in,2"

clean:
	rm -rf generated html build build-html build-* pbefs-manual.pdf *.html *.css *.svg *.4ct *.4tc *.idv \
	    *.lg *.tmp *.xref *.dvi *.run.xml *.bbl *.bcf

.PHONY: all chapter html clean FORCE
FORCE:

# --------------------------------------------------------------------------- #
# ----------------------------- END Makefile -------------------------------- #
# --------------------------------------------------------------------------- #
