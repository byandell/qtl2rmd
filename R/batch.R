setwd("~/Documents/Research/rqtl/qtl2rmd")
RmdFilename <- "mediateQtl2.Rmd"

doBatch <- TRUE

covar_info <- read.csv("local/Jax/islet_secr_for_brian.csv")
if(!file.exists(exvivo <- "local/Jax/exvivo_pheno.rds")) {
  load("~/Documents/Research/attie_alan/DO/AttieDOv2/DerivedData/Attie_islet_secr_data_v5.Rdata")
  saveRDS(dataset.exvivo$pheno, file = exvivo)
}
exvivo_pheno <- readRDS(exvivo)

pheno_names <- covar_info$data_name[-(1:2)]

for(pheno_name in pheno_names) {
  cat(pheno_name, "run analysis\n", file = stderr())

  # Preset any `params` parameters you like
  param_vals <- list(
    pheno_name = pheno_name,
    dataSetup = "dataJaxMadison.R",
    coefType = "coef",
    resultpath = paste0("batch/", pheno_name))
    
  out_filename <- paste0("batch/", pheno_name, ".html")
    
  source("R/runMediation.R")
}

