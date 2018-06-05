params <- list(
  target_name = "GLP1_G83_ins_secrete_gm",
  chr_id = "16",
  pos_Mbp = 90,
  dataSetup = "R/dataJaxMadison.R",
  snpScan = 10,
  offset = 2,
  datapath = 'data',
  resultpath = 'results',
  useIntcov = FALSE)

library(dplyr)
library(ggplot2)

# Take care of any overridden names
for(i in names(params)) {
  if(exists(i))
    params[[i]] <- get(i)
}

target_name <- params$target_name
chr_id <- as.character(params$chr_id)

datapath <- params$datapath

if((resultpath <- params$resultpath) != "") {
  resultpath <- file.path(resultpath, target_name)
  if(!dir.exists(resultpath))
    dir.create(resultpath)
  cat(resultpath, "\n", file = stderr())
  list.files(resultpath)
}

source(params$dataSetup)
kinship <- kinship[[chr_id]]
pheno_data <- pheno_data[, target_name, drop = FALSE]
addcovar <- qtl2pattern::covar_df_mx(covar)
intcovar <- NULL

### Genome scan information

cat("Genome scan information\n", file = stderr())

genoprobs <- query_probs()$probs
map <- query_probs()$map
peaks <- 
  readxl::read_excel(
    file.path("data", "Ex vivo peak summary with QTL count[1].xlsx"),
    col_types = c(rep("text", 2), rep("numeric", 12), "text")) %>%
  dplyr::filter(lodcolumn == target_name,
                chr == chr_id)

pos_Mbp <- peaks$pos[1]
if(is.na(pos_Mbp)) {
  pos_Mbp <- params$pos_Mbp
}
cat("using peak position", pos_Mbp, "\n", file = stderr())  

peak_mar <- qtl2::find_marker(map, chr_id, pos_Mbp)
cat(target_name, chr_id, pos_Mbp, file = stderr())
cat(peak_mar, "\n", file = stderr())

snpScan <- as.numeric(params$snpScan)

### mRNA information

cat("mRNA information\n", file = stderr())

mrna <- query_mrna(chr_id)

# Find peaks for local mRNA that are close to SNP peak.

mrna_ids <- 
  (dplyr::filter(
    mrna$peaks,
    qtl_chr == chr_id,
    abs(qtl_pos - pos_Mbp) <= snpScan,
    gene_start >= pos_Mbp - snpScan,
    gene_end <= pos_Mbp + snpScan))$gene_id

mrna_symbols <- 
  (dplyr::filter(
    mrna$annot,
    id %in% mrna_ids))$symbol

tmpfn <- function(mrna_id) {
  
}