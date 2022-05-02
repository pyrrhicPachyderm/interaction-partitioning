SHELL := /bin/bash
LATEXMK_FLAGS = --pdf --cd
RM := rm -f

doc_raws := proposal.tex article.rnw
supporting_tex_files := references.bib bibliography/references.bib reference-styles/ecology.tex reference-styles/authoryear.tex

doc_pdfs := $(patsubst %.tex,%.pdf,$(patsubst %.rnw,%.pdf,$(doc_raws)))

all: $(doc_pdfs)
.PHONY: all

%-dedented.rnw: dedent-noweb %.rnw
	./$< <$(word 2,$^) >$@
%.tex: %-dedented.rnw
	R -e 'library(knitr);knit("$<","$@")'
%.pdf: %.tex $(supporting_tex_files)
	latexmk $(LATEXMK_FLAGS) --jobname="$(basename $@)" $<

clean:
	@(\
		shopt -s globstar;\
		$(RM) **/*.aux **/*.log **/*.fls **/*.fdb_latexmk;\
		$(RM) **/*.out **/*.bbl **/*.bcf **/*.blg **/*.run.xml;\
		$(RM) **/*.fff **/*.lof;\
		$(RM) **/*.o **/*.d **/*.out;\
		$(RM) **/*-dedented.rnw;\
	)
	@$(RM) $(patsubst %.rnw,%.tex,$(filter %.rnw,$(doc_raws)))
	@$(RM) $(doc_pdfs)
	@$(RM) output/*
.PHONY: clean

spellcheck: $(doc_raws)
	@for file in $^; do \
		aspell check --per-conf=./aspell.conf "$$file" ;\
	done
.PHONY: spellcheck

raw_data_file := TCL_DrosMCT/Data/d_both.csv
processed_data_files := output/focal-vector.data output/response-vector.data output/design-matrix.data
output_file := output/brute.data

test_data_files := output/test-focal-vector.data output/test-response-vector.data output/test-design-matrix.data

$(processed_data_files) &: scripts/reshape-DrosMCT $(raw_data_file)
	./$< $(raw_data_file) $(processed_data_files)

#output_template takes the output file name, the program file name, and the flags.
define output_template =
output/$(1).data: src/$(2).out $$(processed_data_files)
	./$$< $$(processed_data_files) $$@ -p $(3)
endef

$(eval $(call output_template,brute,brute,))
$(eval $(call output_template,rjmcmc-flat,rjmcmc,))
$(eval $(call output_template,rjmcmc-aic,rjmcmc,-a))

article.tex: $(output_file)

$(test_data_files) &: scripts/generate-test-data
	./$< $(test_data_files)
testdata: $(test_data_files)
.PHONY: testdata

brutetest: src/brute.out $(test_data_files)
	./$< $(test_data_files) -
.PHONY: brutetest

article.tex: data/species.csv r/parameters.R r/input-data.R r/post-process.R r/brute-post-process.R r/coclassification-table.R r/mantel-test.R r/dist-matrix.R

#Submodules
bibliography/% reference-styles/% TCL_DrosMCT/% &:
	git submodule update --init

#Secondary with no targets prevents deletion of intermediate files.
.SECONDARY:

##########################################################################################################
#Stuff for making C/C++ files. Making extensive use of Peter Miller's "Recursive Make Considered Harmful".
##########################################################################################################
#`CC` should be `gcc` or `g++`, and `CEXT` should correspondingly `c` or `cpp`
CC := g++
CEXT := cpp
CFLAGS += -std=c++17 -Wall -Wpedantic -Wextra -Wno-unused-parameter -O3 -I/usr/include/eigen3
CLIBS := -lm

#These two should work for most projects, but keep an eye on them.
CDIRS := $(shell find . -name '.git' -prune -o -name '*.$(CEXT)' -printf '%h\n' | sort -u) #All directories containing .c/cpp files, not searching .git
CMAINSRC := $(shell grep -lr --exclude-dir='.git' --include='*.$(CEXT)' -E '^(int|void)\s*main') #All .c/cpp files with a main function, not searching .git

###############################
#A big ugly pile of make hacks.
#All the stuff that might need editing as the project changes should be above here.
###############################
CFLAGS += $(patsubst %,-I%,$(CDIRS))
#Find all C files
CSRC := $(shell find -name '*.$(CEXT)')
#Find all .c/.cpp files except those that lead to excutables (all the ones without mains)
CSUPPSRC := $(shell find -name '*.$(CEXT)' $(shell printf '! -name %s ' $(notdir $(CMAINSRC))))
#Determine the object files
COBJ := $(CSRC:.$(CEXT)=.o)
CSUPPOBJ := $(CSUPPSRC:.$(CEXT)=.o)
#Compile the object files
#This can't be done with the built-in implicit rule as we want to use CFLAGS even for C++
%.o: %.$(CEXT)
	$(CC) $(CFLAGS) -c $< -o $@
#Link the object files
#Note that this assumes all executables depend on all object files, but I can't think how to fix that nicely.
%.out: %.o $(CSUPPOBJ)
	$(CC) $(CFLAGS) $^ -o $@ $(CLIBS)
#Include the autogenerated dependencies files
include $(COBJ:.o=.d)
#Rule to build dependencies files
%.d: %.$(CEXT) depend.sh
	./depend.sh $(CC) $(shell dirname $<) $(CFLAGS) $< > $@
