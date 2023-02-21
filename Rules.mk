$(from_root)-clean:
	@(\
		shopt -s globstar &&\
		cd $(from_root) &&\
		$(RM) **/*.o **/*.d **/*.out;\
	)
	@$(RM) $(from_root)/output/*/*
	@for i in $(from_root)/output/*; do if [ -f "$$i" ]; then $(RM) "$$i"; fi; done
clean: $(from_root)-clean
.PHONY: clean $(from_root)-clean

#Note that this if for GNU time (/usr/bin/time) not bash built-in time.
time_format := "real\t%e\nuser\t%U\nsys\t%S"

#'indv' for individual response, 'pop' for population size response, or 'time' for time series.
processed_indv_data_files = $(patsubst %,$(from_root)/output/$(1)/%,indv-focal-vector.data indv-response-vector.data indv-design-matrix.data)
processed_pop_data_files = $(patsubst %,$(from_root)/output/$(1)/%,pop-focal-vector.data pop-response-vector.data pop-design-matrix.data)
processed_time_data_files = $(patsubst %,$(from_root)/output/$(1)/%,time-id-vector.data time-time-vector.data time-density-matrix.data)
priors_file = $(from_root)/output/$(1)/priors.data

tcl_raw_data := $(from_root)/TCL_DrosMCT/Data/d_both.csv
goldberg_raw_data := $(from_root)/data/goldberg/species.csv $(from_root)/data/goldberg/figure2
carrara_raw_data := $(from_root)/data/carrara/data.xls

cxr_additional_output := $(from_root)/output/cxr/species.csv $(from_root)/output/cxr/min-obs.data

#tcl, test, and carrara use standard per capita fecundity response or time series data, so the default 1 (doubling each generation/day) is a reasonable guess.
cxr_r_guess := 1000 #The response is per capita seed production. Some of the species involved could produce thousands of seeds, so 1000 is a reasonable guess.
goldberg_r_guess := 100 #The response is final mass in milligrams. The largest of the species could grow to about 130 mg, so 100 is a reasonable guess.

#process_data_template takes the dataset abbreviation, the dataset type (e.g. indv, time), the raw data file(s), additional output files for the processing script.
define process_data_template =
$$(call processed_$(2)_data_files,$(1)) $(4) &: $(from_root)/scripts/process-$(1) $(3)
	./$$< $(3) $(2) $$(call processed_$(2)_data_files,$(1)) $(4)
endef

$(eval $(call process_data_template,tcl,indv,$(tcl_raw_data),))
$(eval $(call process_data_template,tcl,pop,$(tcl_raw_data),))
$(eval $(call process_data_template,cxr,indv,,$(cxr_additional_output)))
$(eval $(call process_data_template,goldberg,indv,$(goldberg_raw_data),))
$(eval $(call process_data_template,carrara,time,$(carrara_raw_data),))
$(eval $(call process_data_template,test,indv,,))

#priors_template takes the dataset abbreviation, the dataset type (e.g. indv, time), the error distribution, and options to the prior guessing script.
define priors_template =
$$(call priors_file,$(1)): $(from_root)/scripts/guess-priors $$(call processed_$(2)_data_files,$(1))
	./$$< $$@ $(3) $(2) $$(call processed_$(2)_data_files,$(1)) $(4)
endef

$(eval $(call priors_template,tcl,pop,negativebinomial,))
$(eval $(call priors_template,cxr,indv,negativebinomial,-r $(cxr_r_guess)))
$(eval $(call priors_template,goldberg,indv,gamma,-r $(goldberg_r_guess)))
$(eval $(call priors_template,carrara,time,normal,))

#output_template takes the dataset abbreviation, the dataset type (e.g. indv, time), output file name, the program file name, the model, the error distribution, and the flags.
define output_template =
$(from_root)/output/$(1)/$(3).data $(from_root)/output/$(1)/$(3)-runtime.data &: $(from_root)/src/$(4).out $$(call processed_$(2)_data_files,$(1)) $$(if $$(findstring rjmcmc,$(4)),$$(call priors_file,$(1)))
	$(TIME) -f $(time_format) -o $(from_root)/output/$(1)/$(3)-runtime.data ./$$< $(from_root)/output/$(1)/$(3).data $(5) $(6) $(2) $$(filter-out $$<,$$^) $(7)
endef

$(eval $(call output_template,test,indv,brute,brute,lotkavolterra,normal,))
$(eval $(call output_template,tcl,indv,brute,brute,lotkavolterra,normal,))
$(eval $(call output_template,tcl,pop,rjmcmc,rjmcmc,bevertonholt,negativebinomial,-g -s1000000 -t10))
$(eval $(call output_template,cxr,indv,rjmcmc,rjmcmc,bevertonholt,negativebinomial,-g -d15 -s50000000 -t500))
$(eval $(call output_template,goldberg,indv,rjmcmc,rjmcmc,bevertonholt,gamma,-g -s10000000 -t100))
$(eval $(call output_template,carrara,time,rjmcmc,rjmcmc,lotkavolterra,normal,-g -d15))

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
#`cc` should be `gcc` or `g++`, and `cext` should correspondingly `c` or `cpp`
cc := g++
cext := cpp
cflags := $(CFLAGS) -std=c++20 -Wall -Wpedantic -Wextra -Wno-unused-parameter -Wno-non-template-friend -O3 -fopenmp -I/usr/include/eigen3
clibs := $(CLIBS) -lm

#These two should work for most projects, but keep an eye on them.
cdirs := $(shell find '$(from_root)' -name '.git' -prune -o -name '*.$(cext)' -printf '%h\n' | sort -u) #All directories containing .c/cpp files, not searching .git
cmainsrc := $(shell grep -lr --exclude-dir='.git' --include='*.$(cext)' -E '^(int|void)\s*main' '$(from_root)') #All .c/cpp files with a main function, not searching .git

###############################
#A big ugly pile of make hacks.
#All the stuff that might need editing as the project changes should be above here.
###############################
cflags += $(patsubst %,-I%,$(cdirs))
#Find all C files
csrc := $(shell find '$(from_root)' -name '*.$(cext)')
#Find all .c/.cpp files except those that lead to excutables (all the ones without mains)
csuppsrc := $(shell find '$(from_root)' -name '*.$(cext)' $(shell printf '! -name %s ' $(notdir $(cmainsrc))))
#Determine the object files
cobj := $(csrc:.$(cext)=.o)
csuppobj := $(csuppsrc:.$(cext)=.o)
#Compile the object files
$(from_root)/%.o: $(from_root)/%.$(cext)
	$(cc) $(cflags) -c $< -o $@
#Link the object files
#Note that this assumes all executables depend on all object files, but I can't think how to fix that nicely.
$(from_root)/%.out: $(from_root)/%.o $(csuppobj)
	$(cc) $(cflags) $^ -o $@ $(clibs)
#Include the autogenerated dependencies files
include $(cobj:.o=.d)
#Rule to build dependencies files
$(from_root)/%.d: $(from_root)/%.$(cext) $(from_root)/depend.sh
	./$(word 2,$^) $(cc) $(from_root) $(shell dirname $<) $(cflags) $< > $@
