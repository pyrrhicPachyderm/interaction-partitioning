SHELL := /bin/bash
RM := rm -f
TIME := /usr/bin/time

all: output/article-data-1.rda output/article-data-2.rda
.PHONY: all

here := .
include Rules-$(here).mk

.SECONDEXPANSION:
Rules-%.mk: $$(subst |,/,$$*)/Rules.mk
	sed 's@HERE@$(subst |,/,$*)@g' <'$<' >'$@'
