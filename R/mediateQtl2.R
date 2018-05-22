params <- list(
  target_name = "GLP1_G83_ins_secrete_gm",
  chrID = "1",
  dataSetup = "R/dataJaxMadison.R",
  snpScan = 10,
  offset = 2,
  datapath = 'data',
  resultpath = 'results',
  useIntcov = FALSE)

library(dplyr)
library(ggplot2)

cat(target_name, chr_id, "\n", file = stderr())

# target_name and chr_id might be supplied by calling script
if(!exists("target_name")) {
  target_name <- params$target_name
}
if(!exists("chr_id")) {
  chr_id <- as.character(params$chrID)
}
cat(target_name, chr_id, "\n", file = stderr())

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
target <- pheno_data[, target_name, drop = FALSE]
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
peak_mar <- qtl2::find_marker(map, chr_id, pos_Mbp)
ap <- qtl2::pull_genoprobpos(genoprobs, peak_mar)

### mRNA information

cat("mRNA information\n", file = stderr())

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
mrna.annot <- 
  dplyr::filter(
    mrna$annot,
    id %in% med_signif)
rm(mrna)

### Mediation with allele probs

cat("Mediation with allele probs\n", file = stderr())

mrna.annot$driver <- qtl2::find_marker(map, chr_id, mrna.annot$qtl_pos)
driver_med <- genoprobs[[chr_id]][,,unique(mrna.annot$driver), drop = FALSE]

med_test <- intermediate::mediation_test(
  target   = target,
  mediator = mrna.expr,
  annotation = mrna.annot,
  covar_tar = covar,
  covar_med = covar,
  kinship = kinship,
  intcovar = intcovar,
  driver = NULL,
  driver_med = driver_med)
sum_med <- 
   dplyr::arrange(
     summary(med_test),
     pvalue)

if(resultpath != "")
  write.csv(sum_med, file = file.path(resultpath,
                                      paste0(target_name, "_", chr_id, "_mediation.csv")))
if(resultpath != "")
  write.csv(med_test$test, file = file.path(resultpath,
                                      paste0(target_name, "_", chr_id, "_effect.csv")))

### Association mapping

cat("Association mapping\n", file = stderr())

# Map all of Chr `r chr_id` unless `snpScan` (`r params$snpScan`) is positive.

snpScan <- as.numeric(params$snpScan)
if(is.na(snpScan) || snpScan <= 0) {
  start <- end <- NULL
} else {
  start <- pos_Mbp - snpScan
  end <- pos_Mbp + snpScan
}

pairprobs <- query_probs(chr = chr_id, start = start, stop = end, allele = FALSE)
pairmap <- pairprobs$map
pairprobs <- pairprobs$probs

assoc_ins = qtl2::scan1snps(genoprobs = pairprobs[,chr_id], map = pairmap, 
                      pheno = target, kinship = kinship,
                      addcovar = addcovar, intcovar = intcovar, chr = chr_id, start = start,
                      end = end, query_func = query_variant, cores = 4,
                      keep_all_snps = FALSE)

ts <- qtl2::top_snps(assoc_ins$lod, assoc_ins$snpinfo) %>%
  arrange(desc(lod))

if(resultpath != "")
  write.csv(ts, file = file.path(resultpath,
                                 paste0(target_name, "_", chr_id, "_topsnps.csv")))

### Mediation with SNP probs

cat("Mediation with SNP probs\n", file = stderr())

# Find best SNP for each mediator.

assoc_med = qtl2::scan1snps(genoprobs = pairprobs[,chr_id], map = pairmap, 
                      pheno = mrna.expr,
                      kinship = kinship,
                      addcovar = addcovar, intcovar = intcovar, chr = chr_id, start = start,
                      end = end, query_func = query_variant, cores = 4,
                      keep_all_snps = FALSE)
ts_med <- qtl2pattern::top_snps_all(assoc_med$lod, assoc_med$snpinfo) %>%
  arrange(desc(lod))

ts_med <- ts_med %>%
   group_by(pheno, chr) %>%
   summarize(
     pos = pos[which.max(lod)],
     snp_id = snp_id[which.max(lod)],
     lod = max(lod)) %>%
   ungroup %>%
   arrange(desc(lod))

# Set up drivers for mediators.

peak_snp <- ts$snp_id[1]

m <- match(mrna.annot$id, ts_med$pheno)
mrna.annot$driver <- ts_med$snp_id[m]
driver_med_snp <-
  qtl2::genoprob_to_snpprob(
    pairprobs, 
    assoc_ins$snpinfo)[[chr_id]][,, c( peak_snp, unique(ts_med$snp_id)), drop = FALSE]

med2_test <- intermediate::mediation_test(
  target   = target,
  mediator = mrna.expr,
  annotation = mrna.annot,
  covar_tar = covar,
  covar_med = covar,
  kinship = kinship,
  intcovar = intcovar,
  driver = NULL,
  driver_med = driver_med_snp)
sum_med2 <- summary(med2_test)

if(resultpath != "")
  write.csv(sum_med2, file = file.path(resultpath,
                                       paste0(target_name, "_", chr_id, "_mediation_snp.csv")))
if(resultpath != "")
  write.csv(med2_test$test, file = file.path(resultpath,
                                       paste0(target_name, "_", chr_id, "_effect_snp.csv")))
