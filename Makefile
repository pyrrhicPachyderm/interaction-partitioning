SHELL := /bin/bash
LATEXMK_FLAGS = --pdf --cd
RM := rm -f

maindoc := proposal
supporting_tex_files := bibliography/references.bib reference-styles/authoryear.tex

all: $(maindoc).pdf
.PHONY: all

$(maindoc).pdf: $(maindoc).tex $(supporting_tex_files)
	latexmk $(LATEXMK_FLAGS) --jobname="$(basename $@)" $<

clean:
	@(\
		shopt -s globstar;\
		$(RM) **/*.aux **/*.log **/*.fls **/*.fdb_latexmk;\
		$(RM) **/*.out **/*.bbl **/*.bcf **/*.blg **/*.run.xml;\
	)
	@$(RM) $(maindoc).pdf
.PHONY: clean

spellcheck: $(maindoc).tex
	@for file in $^; do \
		aspell check --per-conf=./aspell.conf "$$file" ;\
	done
.PHONY: spellcheck

#Submodules
bibliography/references.bib reference-styles/authoryear.tex &:
	git submodule update --init
