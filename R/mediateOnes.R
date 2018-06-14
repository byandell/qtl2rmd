#! /path/to/Rscript --vanilla --default-packages=utils
args <- commandArgs(TRUE)
target_name <- args[[1]]
chr_id <- args[[2]]

cat(target_name, chr_id, "\n", file = stderr())

source("R/mediateOne.R")