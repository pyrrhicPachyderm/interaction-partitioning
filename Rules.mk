HERE-clean:
	@(\
		shopt -s globstar &&\
		cd "HERE" &&\
		$(RM) **/*.o **/*.d **/*.out;\
	)
	@$(RM) HERE/output/*/*
	@for i in HERE/output/*; do if [ -f "$$i" ]; then $(RM) "$$i"; fi; done
	@$(RM) HERE/Rules-*.mk
clean: HERE-clean
.PHONY: clean HERE-clean

#Note that this if for GNU time (/usr/bin/time) not bash built-in time.
HERE_time_format := "real\t%e\nuser\t%U\nsys\t%S"

#'indv' for individual response, 'pop' for population size response, or 'time' for time series.
HERE_processed_indv_data_files = $(patsubst %,HERE/output/$(1)/%,indv-focal-vector.data indv-response-vector.data indv-design-matrix.data)
HERE_processed_pop_data_files = $(patsubst %,HERE/output/$(1)/%,pop-focal-vector.data pop-response-vector.data pop-design-matrix.data)
HERE_processed_time_data_files = $(patsubst %,HERE/output/$(1)/%,time-id-vector.data time-time-vector.data time-density-matrix.data)
HERE_priors_file = HERE/output/$(1)/priors.data

HERE_tcl_raw_data := HERE/TCL_DrosMCT/Data/d_both.csv
HERE_goldberg_raw_data := HERE/data/goldberg/species.csv HERE/data/goldberg/figure2
HERE_carrara_raw_data := HERE/data/carrara/data.xls

HERE_cxr_additional_output := HERE/output/cxr/species.csv HERE/output/cxr/min-obs.data

#test and carrara use standard per capita fecundity response or time series data, so the default 1 (doubling each generation/day) is a reasonable guess.
HERE_tcl_r_guess := 10 #From Terry et al. 2021.
HERE_cxr_r_guess := 1000 #The response is per capita seed production. Some of the species involved could produce thousands of seeds, so 1000 is a reasonable guess.
HERE_goldberg_r_guess := 100 #The response is final mass in milligrams. The largest of the species could grow to about 130 mg, so 100 is a reasonable guess.

#process_data_template takes the dataset abbreviation, the dataset type (e.g. indv, time), the raw data file(s), additional output files for the processing script, and the options.
define HERE_process_data_template =
$$(call HERE_processed_$(2)_data_files,$(1)) $(4) &: HERE/scripts/process-$(1) $(3)
	./$$< $(3) $(2) $$(call HERE_processed_$(2)_data_files,$(1)) $(4) $(5)
endef

$(eval $(call HERE_process_data_template,tcl,indv,$(HERE_tcl_raw_data),,))
$(eval $(call HERE_process_data_template,tcl,pop,$(HERE_tcl_raw_data),,))
$(eval $(call HERE_process_data_template,cxr,indv,,$(HERE_cxr_additional_output),-m 4))
$(eval $(call HERE_process_data_template,goldberg,indv,$(HERE_goldberg_raw_data),,))
$(eval $(call HERE_process_data_template,carrara,time,$(HERE_carrara_raw_data),,))
$(eval $(call HERE_process_data_template,test,indv,,))

#priors_template takes the dataset abbreviation, the dataset type (e.g. indv, time), the error distribution, and options to the prior guessing script.
define HERE_priors_template =
$$(call HERE_priors_file,$(1)): HERE/scripts/guess-priors $$(call HERE_processed_$(2)_data_files,$(1))
	./$$< $$@ $(3) $(2) $$(call HERE_processed_$(2)_data_files,$(1)) $(4)
endef

$(eval $(call HERE_priors_template,tcl,pop,negativebinomial,-r $(HERE_tcl_r_guess)))
$(eval $(call HERE_priors_template,cxr,indv,negativebinomial,-r $(HERE_cxr_r_guess)))
$(eval $(call HERE_priors_template,goldberg,indv,gamma,-r $(HERE_goldberg_r_guess)))
$(eval $(call HERE_priors_template,carrara,time,normal,))

#output_template takes the dataset abbreviation, the dataset type (e.g. indv, time), output file name, the program file name, the model, the error distribution, and the flags.
define HERE_output_template =
HERE/output/$(1)/$(3).data HERE/output/$(1)/$(3)-runtime.data &: HERE/src/$(4).out $$(call HERE_processed_$(2)_data_files,$(1)) $$(if $$(findstring rjmcmc,$(4)),$$(call HERE_priors_file,$(1)))
	$(TIME) -f $(HERE_time_format) -o HERE/output/$(1)/$(3)-runtime.data ./$$< HERE/output/$(1)/$(3).data $(5) $(6) $(2) $$(filter-out $$<,$$^) $(7)
endef

$(eval $(call HERE_output_template,test,indv,brute,brute,lotkavolterra,normal,))
$(eval $(call HERE_output_template,tcl,indv,brute,brute,lotkavolterra,normal,))
$(eval $(call HERE_output_template,tcl,pop,rjmcmc,rjmcmc,bevertonholt,negativebinomial,-g -s1000000 -t10))
$(eval $(call HERE_output_template,cxr,indv,rjmcmc,rjmcmc,bevertonholt,negativebinomial,-g -d15 -s10000000 -t100))
$(eval $(call HERE_output_template,goldberg,indv,rjmcmc,rjmcmc,bevertonholt,gamma,-g -s10000000 -t100))
$(eval $(call HERE_output_template,carrara,time,rjmcmc,rjmcmc,lotkavolterra,normal,-g -d15))

define HERE_article_analysis_template =
HERE/output/article-data-$(1).rda: HERE/article-analysis-$(1) $(shell grep -oE '"/[^"]*\.((R)|(csv)|(data))"' HERE/article-analysis-$(1) | sed 's@"\(.*\)"@HERE\1@' | tr '\n' ' ')
	./$$< $$@
endef

$(eval $(call HERE_article_analysis_template,1))
$(eval $(call HERE_article_analysis_template,2))

#Submodules
HERE/TCL_DrosMCT/%:
	git -C "HERE" submodule update --init

#Secondary with no targets prevents deletion of intermediate files.
.SECONDARY:

##########################################################################################################
#Stuff for making C/C++ files. Making extensive use of Peter Miller's "Recursive Make Considered Harmful".
##########################################################################################################
#`cc` should be `gcc` or `g++`, and `cext` should correspondingly `c` or `cpp`
HERE_cc := g++
HERE_cext := cpp
HERE_cflags := $(CFLAGS) -std=c++20 -Wall -Wpedantic -Wextra -Wno-unused-parameter -Wno-non-template-friend -O3 -fopenmp -I/usr/include/eigen3
HERE_clibs := $(CLIBS) -lm

#These two should work for most projects, but keep an eye on them.
HERE_cdirs := $(shell find 'HERE' -name '.git' -prune -o -name '*.$(HERE_cext)' -printf '%h\n' | sort -u) #All directories containing .c/cpp files, not searching .git
HERE_cmainsrc := $(shell grep -lr --exclude-dir='.git' --include='*.$(HERE_cext)' -E '^(int|void)\s*main' 'HERE') #All .c/cpp files with a main function, not searching .git

###############################
#A big ugly pile of make hacks.
#All the stuff that might need editing as the project changes should be above here.
###############################
HERE_cflags += $(patsubst %,-I%,$(HERE_cdirs))
#Find all C files
HERE_csrc := $(shell find 'HERE' -name '*.$(HERE_cext)')
#Find all .c/.cpp files except those that lead to excutables (all the ones without mains)
HERE_csuppsrc := $(shell find 'HERE' -name '*.$(HERE_cext)' $(shell printf '! -name %s ' $(notdir $(HERE_cmainsrc))))
#Determine the object files
HERE_cobj := $(HERE_csrc:.$(HERE_cext)=.o)
HERE_csuppobj := $(HERE_csuppsrc:.$(HERE_cext)=.o)
#Compile the object files
HERE/%.o: HERE/%.$(HERE_cext)
	$(HERE_cc) $(HERE_cflags) -c $< -o $@
#Link the object files
#Note that this assumes all executables depend on all object files, but I can't think how to fix that nicely.
HERE/%.out: HERE/%.o $(HERE_csuppobj)
	$(HERE_cc) $(HERE_cflags) $^ -o $@ $(HERE_clibs)
#Include the autogenerated dependencies files
#They are entirely targets and prerequisites, which are expanded immediately, so they will interpret $(here) correctly.
here := HERE
include $(HERE_cobj:.o=.d)
#Rule to build dependencies files
HERE/%.d: HERE/%.$(HERE_cext) HERE/depend.sh
	./$(word 2,$^) $(HERE_cc) HERE $(shell dirname $<) $(HERE_cflags) $< > $@
