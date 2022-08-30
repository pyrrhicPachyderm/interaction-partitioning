$(from_root)-clean:
	@(\
		shopt -s globstar &&\
		cd $(from_root) &&\
		$(RM) **/*.o **/*.d **/*.out;\
	)
	@$(RM) $(from_root)/output/*
.PHONY: $(from_root)-clean

raw_data_file := $(from_root)/TCL_DrosMCT/Data/d_both.csv
processed_data_files := $(patsubst %,$(from_root)/output/%,focal-vector.data response-vector.data design-matrix.data)
output_file := $(from_root)/output/brute.data

$(processed_data_files) &: $(from_root)/scripts/reshape-DrosMCT $(raw_data_file)
	./$< $(raw_data_file) $(processed_data_files)

#output_template takes the output file name, the program file name, and the flags.
define output_template =
$(from_root)/output/$(1).data: $(from_root)/src/$(2).out $$(processed_data_files)
	./$$< $$(processed_data_files) $$@ -p $(3)
endef

$(eval $(call output_template,brute,brute,))
$(eval $(call output_template,rjmcmc-flat,rjmcmc,))
$(eval $(call output_template,rjmcmc-aic,rjmcmc,-a))

additional_data_files := $(from_root)/data/species.csv
r_source_files := $(patsubst %,$(from_root)/r/%,parameters.R input-data.R post-process.R brute-post-process.R coclassification-table.R grouped-matrix.R mantel-test.R dist-matrix.R)

$(from_root)/output/article-data.rda: $(from_root)/article-analysis $(processed_data_files) $(output_file) $(additional_data_files) $(r_source_files)
	./$< $@

#Test data analysis.

test_data_files := $(patsubst %,$(from_root)/ouput/%,test-focal-vector.data test-response-vector.data test-design-matrix.data)

$(test_data_files) &: $(from_root)/scripts/generate-test-data
	./$< $(test_data_files)
testdata: $(test_data_files)
.PHONY: testdata

brutetest: $(from_root)/src/brute.out $(test_data_files)
	./$< $(test_data_files) -
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
CFLAGS += -std=c++17 -Wall -Wpedantic -Wextra -Wno-unused-parameter -O3 -I/usr/include/eigen3
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
