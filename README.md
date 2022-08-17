# Code for simplifying species interaction models by grouping species

## Contents
- `src` contains the C++ code for the program to perform model-fitting.
- `scripts` contains standalone R scripts for pre-processing the data for the model-fitting program.
- `r` contains R code defining functions helpful for post-processing the output of the model-fitting program.
- `article-analysis` is an R script that post-processes the output of the model-fitting program into data and figures for an article, and saves them as a .rda file.
- `Makefile` defines the relationship between the data, the code, and the code's various outputs.

## Usage
- Simply running `make` will compile all the code and run all the analyses, producing `output/article-data.rda`, which may then be loaded in R for examination.
