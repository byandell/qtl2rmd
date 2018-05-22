#! /path/to/Rscript --vanilla --default-packages=utils
args <- commandArgs(TRUE)
target_name <- args[[1]]
chr_id <- args[[2]]
simMediate <- as.numeric(args[[3]])
resultpath <- "sims"

cat(target_name, chr_id, simMediate, resultpath, "\n", file = stderr())

source("R/mediateQtl2.R")