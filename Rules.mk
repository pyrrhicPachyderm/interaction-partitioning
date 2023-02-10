$(from_root)-clean:
	@(\
		shopt -s globstar &&\
		cd $(from_root) &&\
		$(RM) **/*.o **/*.d **/*.out;\
	)
	@$(RM) $(from_root)/output/*/*
	@for i in $(from_root)/output/*; do if [ -f "$$i" ]; then $(RM) "$$i"; fi; done
.PHONY: $(from_root)-clean

#'indv' for individual response or 'time' for time series.
processed_indv_data_files = $(patsubst %,$(from_root)/output/$(1)/%,focal-vector.data response-vector.data design-matrix.data)
processed_time_data_files = $(patsubst %,$(from_root)/output/$(1)/%,id-vector.data time-vector.data density-matrix.data)
priors_file = $(from_root)/output/$(1)/priors.data

tcl_raw_data := $(from_root)/TCL_DrosMCT/Data/d_both.csv
goldberg_raw_data := $(from_root)/data/goldberg/species.csv $(from_root)/data/goldberg/figure2
carrara_raw_data := $(from_root)/data/carrara/data.xls

cxr_additional_output := $(from_root)/output/cxr/species.csv

#tcl, test, and carrara use standard per capita fecundity response or time series data, so the default 1 (doubling each generation/day) is a reasonable guess.
cxr_r_guess := 1000 #The response is per capita seed production. Some of the species involved could produce thousands of seeds, so 1000 is a reasonable guess.
goldberg_r_guess := 100 #The response is final mass in milligrams. The largest of the species could grow to about 130 mg, so 100 is a reasonable guess.

#process_data_template takes the dataset abbreviation, the dataset type (indv or time), the raw data file(s), additional output files for the processing script, and options to the prior guessing script.
define process_data_template =
$$(call processed_$(2)_data_files,$(1)) $(4) &: $(from_root)/scripts/process-$(1) $(3)
	./$$< $(3) $$(call processed_$(2)_data_files,$(1)) $(4)
$$(call priors_file,$(1)): $(from_root)/scripts/guess-priors $$(call processed_$(2)_data_files,$(1))
	./$$< $$@ $(2) $$(call processed_$(2)_data_files,$(1)) $(5)
endef

$(eval $(call process_data_template,tcl,indv,$(tcl_raw_data),,))
$(eval $(call process_data_template,cxr,indv,,$(cxr_additional_output),-r $(cxr_r_guess)))
$(eval $(call process_data_template,goldberg,indv,$(goldberg_raw_data),,-r $(goldberg_r_guess)))
$(eval $(call process_data_template,carrara,time,$(carrara_raw_data),,))
$(eval $(call process_data_template,test,indv,,,))

#output_template takes the dataset abbreviation, the dataset type (indv or time), output file name, the program file name, and the flags.
define output_template =
$(from_root)/output/$(1)/$(3).data: $(from_root)/src/$(4).out $$(call processed_$(2)_data_files,$(1))
	./$$< $$(call processed_$(2)_data_files,$(1)) $$@ $(5)
endef

$(eval $(call output_template,test,indv,brute,brute,))
$(eval $(call output_template,tcl,indv,brute,brute,))
$(eval $(call output_template,tcl,indv,rjmcmc,rjmcmc,))
$(eval $(call output_template,cxr,indv,rjmcmc,rjmcmc,))
$(eval $(call output_template,goldberg,indv,rjmcmc,rjmcmc,))

define article_analysis_template =
$(from_root)/output/article-data-$(1).rda: $(from_root)/article-analysis-$(1) $(shell grep -oE '"/[^"]*\.((R)|(csv)|(data))"' $(from_root)/article-analysis-$(1) | sed 's/"\(.*\)"/$(from_root)\1/' | tr '\n' ' ')
	./$$< $$@
endef

$(eval $(call article_analysis_template,1))
$(eval $(call article_analysis_template,2))

#Test data analysis.

brutetest: $(from_root)/output/test/brute.data
.PHONY: brutetest

#Submodules
$(from_root)/TCL_DrosMCT/%:
	git -C "$(from_root)" submodule update --init

#Secondary with no targets prevents deletion of intermediate files.
.SECONDARY:

##########################################################################################################
#Stuff for making C/C++ files. Making extensive use of Peter Miller's "Recursive Make Considered Harmful".
##########################################################################################################
#`CC` should be `gcc` or `g++`, and `CEXT` should correspondingly `c` or `cpp`
CC := g++
CEXT := cpp
CFLAGS += -std=c++20 -Wall -Wpedantic -Wextra -Wno-unused-parameter -Wno-non-template-friend -O3 -I/usr/include/eigen3
CLIBS := -lm

#These two should work for most projects, but keep an eye on them.
CDIRS := $(shell find '$(from_root)' -name '.git' -prune -o -name '*.$(CEXT)' -printf '%h\n' | sort -u) #All directories containing .c/cpp files, not searching .git
CMAINSRC := $(shell grep -lr --exclude-dir='.git' --include='*.$(CEXT)' -E '^(int|void)\s*main' '$(from_root)') #All .c/cpp files with a main function, not searching .git

###############################
#A big ugly pile of make hacks.
#All the stuff that might need editing as the project changes should be above here.
###############################
CFLAGS += $(patsubst %,-I%,$(CDIRS))
#Find all C files
CSRC := $(shell find '$(from_root)' -name '*.$(CEXT)')
#Find all .c/.cpp files except those that lead to excutables (all the ones without mains)
CSUPPSRC := $(shell find '$(from_root)' -name '*.$(CEXT)' $(shell printf '! -name %s ' $(notdir $(CMAINSRC))))
#Determine the object files
COBJ := $(CSRC:.$(CEXT)=.o)
CSUPPOBJ := $(CSUPPSRC:.$(CEXT)=.o)
#Compile the object files
#This can't be done with the built-in implicit rule as we want to use CFLAGS even for C++
$(from_root)/%.o: $(from_root)/%.$(CEXT)
	$(CC) $(CFLAGS) -c $< -o $@
#Link the object files
#Note that this assumes all executables depend on all object files, but I can't think how to fix that nicely.
$(from_root)/%.out: $(from_root)/%.o $(CSUPPOBJ)
	$(CC) $(CFLAGS) $^ -o $@ $(CLIBS)
#Include the autogenerated dependencies files
include $(COBJ:.o=.d)
#Rule to build dependencies files
$(from_root)/%.d: $(from_root)/%.$(CEXT) $(from_root)/depend.sh
	./$(word 2,$^) $(CC) $(from_root) $(shell dirname $<) $(CFLAGS) $< > $@
