# Minimal makefile for Sphinx documentation
#

SPHINXOPTS    =
ifdef OS
   SPHINXBUILD   = python3 -msphinx
else
   SPHINXBUILD   = python3 -msphinx
endif
SPHINXPROJ    = PoreSim
SOURCEDIR     = .
BUILDDIR      = _build

help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile

%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
