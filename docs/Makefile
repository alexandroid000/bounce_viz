# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line.
SPHINXOPTS    =
SPHINXBUILD   = uv run sphinx-build
SPHINXAPIDOC  = uv run sphinx-apidoc
SPHINXPROJ    = BounceVisualization
SOURCEDIR     = source
BUILDDIR      = build

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

apidoc: clean
	@rm -f ./source/*rst*
	@echo "Removing all .rst files under 'source'..."
	@$(SPHINXAPIDOC) -o ./source/ ../../bounce_viz
	@rm -f ./source/modules.rst
	@find ./source -name '*.rst' -type f -exec vim -c '%s/^\([a-zA-Z0-9\_]\+[.]\)*\([a-zA-Z0-9\_]\+\) \(module\|package\)/\2/g' -c 'wq' -- {} \;
	@make -j html 
