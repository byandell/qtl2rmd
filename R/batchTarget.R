setwd("~/Documents/Research/rqtl/qtl2rmd")
RmdFilename <- "mediateQtl2.Rmd"

doBatch <- TRUE

datapath <- "data"
covar_info <- read.csv(file.path(datapath, "islet_secr_for_brian.csv"))
peak_info <- 
  readxl::read_excel(
    file.path(datapath, "Ex vivo peak summary with QTL count[1].xlsx"),
    col_types = c(rep("text", 2), rep("numeric", 12), "text"))
    
if(!file.exists(exvivo <- file.path(datapath, "exvivo_pheno.rds"))) {
  load("~/Documents/Research/attie_alan/DO/AttieDOv2/DerivedData/Attie_islet_secr_data_v5.Rdata")
  saveRDS(dataset.exvivo$pheno, file = exvivo)
}
exvivo_pheno <- readRDS(exvivo)

pheno_names <- peak_info$lodcolumn
chr_names <- peak_info$chr

for(i in seq_along(pheno_names)) {
  cat(pheno_names[i], chr_names[i], "run analysis\n", file = stderr())

  # Preset any `params` parameters you like
  param_vals <- list(
    target_name = pheno_names[i],
    dataSetup = "R/dataJaxMadison.R",
    coefType = "coef",
    chrID = chr_names[i],
    resultpath = file.path("batch", pheno_names[i]))
    
  out_filename <- paste0(param_vals$resultpath, "/",
                         pheno_names[i], "_", chr_names[i], ".html")
    
  source("R/runMediation.R")
}

