setwd("~/Documents/Research/rqtl/qtl2rmd")
RmdFilename <- "Rmd/mediateOne.Rmd"

doBatch <- TRUE

targets <- read.csv("data/targets.csv", stringsAsFactors = FALSE)

for(i in seq_len(nrow(targets))) {
  cat(targets$target[i], targets$chr[i], targets$mediator[i], "run analysis\n",
      file = stderr())

  # Preset any `params` parameters you like
  param_vals <- list(
    dataSetup = "dataJaxMadison.R",
    target_name = targets$target[i],
    chr_id = targets$chr[i],
    pos_Mbp = targets$pos[i],
    mediator_name = targets$mediator[i],
    datapath = "../data",
    resultpath = "../one")
    
  out_filename <- paste0(param_vals$resultpath, "/",
                         targets$target[i], "_",
                         targets$chr[i], "_",
                         targets$mediator[i], ".html")
    
  source("R/runOne.R")
}

