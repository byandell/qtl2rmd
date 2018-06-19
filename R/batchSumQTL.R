setwd("~/Documents/Research/rqtl/qtl2rmd")
RmdFilename <- "Rmd/summaryQTL.Rmd"

doBatch <- TRUE

targets_QTL <- 
  dplyr::filter(
    read.table("data/exvivo_peak_summary.txt",
               stringsAsFactors = FALSE,
               header = TRUE, fill = TRUE),
    !is.na(QTL))

for(i in seq_len(nrow(targets_QTL))) {
  cat("chr", targets_QTL$chr[i], "QTL", targets_QTL$QTL[i], "interactive plot\n",
      file = stderr())

  # Preset any `params` parameters you like
  param_vals <- list(
    QTL = targets_QTL$QTL[i],
    datapath = "../results/med_qtl2",
    resultpath = "../results/med_qtl2/plotly")
    
  out_filename <- paste0(param_vals$resultpath, "/",
                         "chr_", targets_QTL$chr[i],
                         "_QTL_", targets_QTL$QTL[i],
                         "_med_qtl2.html")
    
  source("R/runOne.R")
}

