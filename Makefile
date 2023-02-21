SHELL := /bin/bash
RM := rm -f
TIME := /usr/bin/time

all: output/article-data-1.rda output/article-data-2.rda
.PHONY: all

from_root := .
include $(from_root)/Rules.mk
