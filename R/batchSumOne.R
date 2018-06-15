setwd("~/Documents/Research/rqtl/qtl2rmd")
RmdFilename <- "Rmd/summaryOne.Rmd"

doBatch <- TRUE

targets <- read.table("data/exvivo_peak_summary.txt", stringsAsFactors = FALSE,
                      header = TRUE, fill = TRUE)

for(i in seq_len(nrow(targets))) {
  cat(targets$lodcolumn[i], targets$chr[i], "run analysis\n",
      file = stderr())

  # Preset any `params` parameters you like
  param_vals <- list(
    dataSetup = "dataJaxMadison.R",
    target_name = targets$lodcolumn[i],
    chr_id = targets$chr[i],
    datapath = "../results/med_qtl2",
    resultpath = "../results/med_qtl2")
    
  out_filename <- paste0(param_vals$resultpath, "/",
                         targets$lodcolumn[i], "_",
                         targets$chr[i],
                         "_med_qtl2.html")
    
  source("R/runOne.R")
}

