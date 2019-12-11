#
#				Makefile.am
#
# Makefile for Sphinx documentation.
# Process this file with autoconf to generate a Makefile
#

EXTRA_DIST	= requirements.txt TODO src theme

# You can set these variables from the command line.
SPHINXOPTS    =
SPHINXBUILD   = sphinx-build
SPHINXPROJ    = @PACKAGE_NAME@
SOURCEDIR     = src
BUILDDIR      = build

.PHONY: all clean-am mostlyclean-am help html dirhtml singlehtml pickle json htmlhelp qthelp devhelp epub latex latexpdf latexpdfja text man texinfo info gettext changes xml pseudoxml linkcheck doctest coverage

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).

clean:
	rm -rf $(BUILDDIR)/*

pdf:	latexpdf

html dirhtml singlehtml pickle json htmlhelp qthelp devhelp epub latex latexpdf latexpdfja text man texinfo info gettext changes xml pseudoxml linkcheck doctest coverage: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

