SHELL := /bin/bash
RM := rm -f

all: output/article-data.rda
.PHONY: all

from_root := .
include $(from_root)/Rules.mk

clean: $(from_root)-clean
.PHONY: clean
