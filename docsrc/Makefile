# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line.
SPHINXOPTS    =
ifdef OS
   SPHINXBUILD   = python -msphinx
   RM = cmd.exe /c rd /s /q
   CP = cmd.exe /c xcopy /s /e /h /i /y
   TOUCH = fsutil file createnew
   TOUCHLEN = 0
   DOCS = ..\docs
   BUILDHTML = _build\html\.
else
   SPHINXBUILD   = python3 -msphinx
   RM = rm -r
   CP = cp -a
   TOUCH = touch
   TOUCHLEN =
   DOCS = ../docs
   BUILDHTML = _build/html/.
endif
# SPHINXBUILD   = python3 -msphinx
SPHINXPROJ    = MolDyn
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

# Github documentation copy
github:
	@$(RM) $(DOCS)
	@mkdir $(DOCS)
	@$(TOUCH) $(DOCS)/.nojekyll $(TOUCHLEN)
	@make html
	@$(CP) $(BUILDHTML) $(DOCS)
