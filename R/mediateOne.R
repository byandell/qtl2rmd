params <- list(
  target_name = "AA_G83_ins_secrete_gm",
  chr_id = "1",
  pos_Mbp = NULL,
  dataSetup = "R/dataJaxMadison.R",
  snpScan = 10,
  SNPlevels = 2,
  offset = 2,
  datapath = 'data',
  resultpath = 'results',
  useIntcov = FALSE)

library(dplyr)

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

dataSetup <- params$dataSetup
source(dataSetup)

target <- as.matrix(pheno_data[, target_name, drop = FALSE])
rm(pheno_data)
addcovar <- qtl2pattern::covar_df_mx(covar)
kinship <- kinship[[chr_id]]

### Peak position for target

peaks <- readxl::read_excel(file.path(datapath,"Ex vivo peak summary with QTL count[1].xlsx"),
                            col_types = c(rep("text", 2), rep("numeric", 12), "text")) %>%
  dplyr::filter(lodcolumn == target_name,
                chr == chr_id)
pos_Mbp <- peaks$pos[1]
if(is.na(pos_Mbp)) {
  pos_Mbp <- params$pos_Mbp
  stopifnot(!is.null(pos_Mbp))
}
cat("using peak position", pos_Mbp, "\n", file = stderr())  

snpScan <- as.numeric(params$snpScan)
if(is.na(snpScan) || snpScan <= 0) {
  start <- end <- NULL
} else {
  start <- pos_Mbp - snpScan
  end <- pos_Mbp + snpScan
}

# Reduce to common IDs
m <- qtl2::get_common_ids(target,
                          addcovar,
                          complete.cases = TRUE)
target <- target[m,, drop = FALSE]
kinship <- kinship[m,m]
addcovar <- addcovar[m,, drop=FALSE]

### Mediator information

mrna <- query_mrna(chr_id)

# Find peaks for local mRNA that are close to SNP peak.

med_signif <- 
  (dplyr::filter(
    mrna$peaks,
    qtl_chr == chr_id,
    abs(qtl_pos - pos_Mbp) <= params$snpScan,
    gene_start >= pos_Mbp - params$snpScan,
    gene_end <= pos_Mbp + params$snpScan))$gene_id

mrna.expr <- mrna$expr[, med_signif, drop = FALSE]
m <- match(med_signif, mrna$annot$id)
mrna.annot <- mrna$annot[m,, drop = FALSE]
rm(mrna)

#################################################

if(params$SNPlevels == 3) {
  genoprobs <- query_probs(chr = chr_id, start = start, stop = end, allele = FALSE)
  map <- genoprobs$map
  genoprobs <- genoprobs$probs
} else {
  genoprobs <- query_probs(chr = chr_id, start = start, stop = end)
  map <- genoprobs$map
  genoprobs <- genoprobs$probs
}

target_scan <-
  qtl2::scan1snps(
    genoprobs = genoprobs[,chr_id],
    map = map, 
    pheno = target,
    kinship = kinship,
    addcovar = addcovar,
    chr = chr_id, start = start, end = end,
    query_func = query_variant,
    cores = cores,
    keep_all_snps = FALSE)

out <- list()

for(mediator_id in med_signif) {
  mediator_name <- (mrna.annot %>% filter(id == mediator_id))$symbol
  cat(mediator_id, mediator_name, "\n", file = stderr())

  mediator <- as.matrix(mrna.expr[, mediator_id, drop = FALSE])
  annotation <- 
    dplyr::filter(
      mrna.annot,
      id == mediator_id)
  
  # Mediation test on SNPs in peak region.
  # Association mapping for causal models
  
  out[[mediator_id]] <- 
    intermediate::mediation_qtl2(target,
                                 mediator,
                                 annotation, covar, covar, kinship,
                                 genoprobs, map,
                                 drop_lod = 1.5, query_variant,
                                 cores = 4, target_scan)
}

class(out) <- c("listof_mediation_index", class(out))
saveRDS(out, 
        file = file.path(resultpath,
                         paste0(target_name, "_", chr_id, "_med_qtl2.rds")))
