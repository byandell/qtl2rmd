## Set upAttie data

## This assumes below that covariates and data transform are fixed.
## See dataMadison.R for alternative

# Load Attie data.
cat("load Attie data\n", file = stderr())

# Only getting pheno and covar from Jax.

if(!file.exists(exvivo <- "local/Jax/exvivo_pheno.rds")) {
  cat("one time load of Jax data\n", file = stderr())
  load("~/Documents/Research/attie_alan/DO/AttieDOv2/DerivedData/Attie_islet_secr_data_v5.Rdata")
  saveRDS(dataset.exvivo$pheno, file = exvivo)
}

# Assume data already transformed.
pheno_data <- readRDS(exvivo)

covar_info <- read.csv("local/Jax/islet_secr_for_brian.csv")
covar_names <- unlist(stringr::str_split(
  (covar_info %>% dplyr::filter(data_name == target_name))$covar,":"))
form <- formula(paste("~", paste(covar_names, collapse = "+")))
covar = model.matrix(form, data = pheno_data)[,-1]

project_info <- data.frame(project = "AttieDOv2",
                           taxa = "CCmouse",
                           directory = datapath,
                           stringsAsFactors = FALSE)
taxa_dir <- file.path(project_info$directory, 
                      project_info$taxa)
cat(project_dir <- file.path(taxa_dir, 
                          project_info$project), "\n", file = stderr())

kinship <- readRDS(file.path(project_dir, "kinship.rds"))

cat("kinship\n", file = stderr())

if(params$showPeaks > 0)
  peaks <- readRDS(file.path(project_dir, "peaks.rds"))

#############################################################################

## Query functions

## Set up SNP variant query functions.

cc_dbfile = file.path(taxa_dir, "cc_variants.sqlite")
query_variant <- qtl2::create_variant_query_func(cc_dbfile, filter = "type=='snp'")

## Set up gene query functions.

## MGI query

gene_dbfile = file.path(taxa_dir, "mouse_genes.sqlite")
query_gene <- qtl2::create_gene_query_func(dbfile = gene_dbfile,
                                           filter="(source=='MGI')")

## List environment for query_gene

#print(ls.str(environment(query_gene)))

## Query genoprobs

query_probs <- qtl2pattern::create_probs_query_func(project_dir, method = "fast")

## query mRNA

query_mrna <- qtl2pattern::create_mrna_query_func(project_dir)
