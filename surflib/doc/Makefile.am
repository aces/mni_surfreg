
if HAVE_DOXYGEN
dox: html/index.html
else
dox:
endif

all-local: dox

# The dependencies are incomplete, since doc files depend on the
# source files.  Until that is fixed, do "make clean" before "make"
# if you want docs for updated sources.
html/index.html: Doxyfile
	@DOXYGEN@


clean-local:
	rm -rf html man

