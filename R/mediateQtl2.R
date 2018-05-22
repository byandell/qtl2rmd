params <- list(
  target_name = "GLP1_G83_ins_secrete_gm",
  chr_id = "1",
  dataSetup = "R/dataJaxMadison.R",
  snpScan = 10,
  offset = 2,
  datapath = 'data',
  resultpath = 'results',
  useIntcov = FALSE)

library(dplyr)
library(ggplot2)

cat(target_name, chr_id, "\n", file = stderr())

# Take care of any overridden names
for(i in names(params)) {
  if(exists(i))
    params[[i]] <- get(i)
}

target_name <- params$target_name
chr_id <- as.character(params$chr_id)

cat(target_name, chr_id, "\n", file = stderr())


datapath <- params$datapath

if((resultpath <- params$resultpath) != "") {
  resultpath <- file.path(resultpath, target_name)
  if(!dir.exists(resultpath))
    dir.create(resultpath)
  cat(resultpath, "\n", file = stderr())
  list.files(resultpath)
}

if(!exists("simMediate")) {
  simMediate <- 0
  csvEnd <- ".csv"
} else {
  csvEnd <- paste0("_", simMediate, ".csv")
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
peak_mar <- qtl2::find_marker(map, chr_id, pos_Mbp)
cat(target_name, chr_id, pos_Mbp, file = stderr())
cat(peak_mar, "\n", file = stderr())
driver_tar <- qtl2::pull_genoprobpos(genoprobs, peak_mar)

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
m <- match(med_signif, mrna$annot$id)
mrna.annot <- mrna$annot[m,, drop = FALSE]
rm(mrna)

### Mediation with allele probs

cat("Mediation with allele probs\n", file = stderr())

mrna.annot$driver <- qtl2::find_marker(map, chr_id, mrna.annot$qtl_pos)
driver_med <- genoprobs[[chr_id]][,,unique(mrna.annot$driver), drop = FALSE]

target <- pheno_data
mediator <- mrna.expr

if(simMediate) {
  cat("simulate by shuffling residuals\n", file = stderr())
  source("R/shuffleQtl2.R")
  target <- shuffleQtl2(driver_tar, target, kinship, covar)
  mediator <- shuffleQtl2M(driver_med, mediator, kinship, covar, mrna.annot$driver)
}

med_test <- intermediate::mediation_test(
  target   = target,
  mediator = mediator,
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
                                      paste0(target_name, "_", chr_id, "_med_test_allele", csvEnd)))
if(resultpath != "" & !simMediate)
  write.csv(med_test$fit, file = file.path(resultpath,
                                      paste0(target_name, "_", chr_id, "_med_fit_allele.csv")))

tar_test <- intermediate::mediation_test(
  target   = target,
  mediator = mediator,
  annotation = mrna.annot,
  covar_tar = covar,
  covar_med = covar,
  kinship = kinship,
  intcovar = intcovar,
  driver = driver_tar,
  driver_med = driver_med)
sum_tar <- 
  dplyr::arrange(
    summary(tar_test),
    pvalue)

if(resultpath != "")
  write.csv(sum_tar, file = file.path(resultpath,
                                      paste0(target_name, "_", chr_id, "_tar_test_allele", csvEnd)))
if(resultpath != "" & !exists("simMediate"))
  write.csv(tar_test$test, file = file.path(resultpath,
                                            paste0(target_name, "_", chr_id, "_tar_fit_allele.csv")))

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

#pairprobs <- query_probs(chr = chr_id, start = start, stop = end, allele = FALSE)
#pairmap <- pairprobs$map
#pairprobs <- pairprobs$probs

assoc_ins = qtl2::scan1snps(genoprobs = genoprobs[,chr_id], map = map, 
                      pheno = target, kinship = kinship,
                      addcovar = addcovar, intcovar = intcovar, chr = chr_id, start = start,
                      end = end, query_func = query_variant, cores = 4,
                      keep_all_snps = FALSE)

ts <- qtl2::top_snps(assoc_ins$lod, assoc_ins$snpinfo) %>%
  arrange(desc(lod))

### Mediation with SNP probs

cat("Mediation with SNP probs\n", file = stderr())

# Find best SNP for each mediator.

assoc_med = qtl2::scan1snps(genoprobs = genoprobs[,chr_id], map = map, 
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
tmp <- qtl2::genoprob_to_snpprob(genoprobs, assoc_ins$snpinfo)
driver_tar_snp <- tmp[[chr_id]][,, peak_snp]
driver_med_snp <- tmp[[chr_id]][,, unique(ts_med$snp_id), drop = FALSE]

target <- pheno_data
mediator <- mrna.expr

if(exists("simMediate")) {
  cat("simulate by shuffling residuals\n", file = stderr())
  source("R/shuffleQtl2.R")
  target <- shuffleQtl2(driver_tar_snp, target, kinship, covar)
  mediator <- shuffleQtl2M(driver_med_snp, mediator, kinship, covar, mrna.annot$driver)
}

med2_test <- intermediate::mediation_test(
  target   = target,
  mediator = mediator,
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
                                       paste0(target_name, "_", chr_id, "_med_test_snp", csvEnd)))
if(resultpath != "" & !exists("simMediate"))
  write.csv(med2_test$fit, file = file.path(resultpath,
                                       paste0(target_name, "_", chr_id, "_med_fit_snp.csv")))

tar2_test <- intermediate::mediation_test(
  target   = target,
  mediator = mediator,
  annotation = mrna.annot,
  covar_tar = covar,
  covar_med = covar,
  kinship = kinship,
  intcovar = intcovar,
  driver = driver_tar_snp,
  driver_med = driver_med_snp)
sum_tar2 <- summary(tar2_test)

if(resultpath != "")
  write.csv(sum_tar2, file = file.path(resultpath,
                                       paste0(target_name, "_", chr_id, "_tar_test_snp", csvEnd)))
if(resultpath != "" & !exists("simMediate"))
  write.csv(med2_test$test, file = file.path(resultpath,
                                             paste0(target_name, "_", chr_id, "_tar_fit_snp.csv")))
