EXTRA_DIST = theory.lyx theory.tex theory.bbl \
	quadrisection.eps \
	mean-cortex-top.ps \
	mean-sphere-top.ps \
	fractional-depth.eps \
	approximate-path.eps \
	sphere-top.ps sphere-bottom.ps sphere-left.ps \
	sphere-right.ps \
	mappings.eps \
	trianglemap-sphere.eps \
	00111_smooth_128.ps \
	00111_smooth_32.ps \
	00111_smooth_8.ps \
	00111_smooth_2.ps \
	similarity-neighbourhood.eps \
	search-neighbourhood.eps \
	surfacepoint.eps

SUFFIXES = .tex .dvi .pdf

TEXINPUTS = .:$(srcdir):
export TEXINPUTS

.PHONY: all
all: theory.ps

.tex.dvi:
	latex $<
	latex $<

theory.ps: theory.dvi
	dvips $? -o $@

CLEANFILES = theory.log theory.aux theory.dvi theory.ps
