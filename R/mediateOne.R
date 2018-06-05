params <- list(
  target_name = "GLP1_G83_ins_secrete_gm",
  chr_id = "16",
  mediator_name = "",
  mediator_id = "",
  pos_Mbp = 90,
  dataSetup = "R/dataJaxMadison.R",
  snpScan = 10,
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
mediator_id <- params$mediator_id
mediator_name <- params$mediator_name

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

map <- query_probs(chr_id)
genoprobs <- map$probs
map <- map$map
kinship <- kinship[[chr_id]]

### Peak position for target

peaks <- readxl::read_excel(file.path(datapath,"Ex vivo peak summary with QTL count[1].xlsx"),
                            col_types = c(rep("text", 2), rep("numeric", 12), "text")) %>%
  dplyr::filter(lodcolumn == target_name,
                chr == chr_id)
pos_Mbp <- peaks$pos[1]
if(is.na(pos_Mbp)) {
  pos_Mbp <- params$pos_Mbp
}
cat("using peak position", pos_Mbp, "\n", file = stderr())  

snpScan <- as.numeric(params$snpScan)
if(is.na(snpScan) || snpScan <= 0) {
  start <- end <- NULL
} else {
  start <- pos_Mbp - snpScan
  end <- pos_Mbp + snpScan
}

### Mediator information

mrna <- query_mrna(chr_id)

# Must supply one of mediator id or name.
if(mediator_id != "") {
  mrna_id <- mediator_id
  if(mediator_name != "") {
    mrna_name <- 
      (mrna$annot %>%
         filter(id == mediator_id))$symbol
    stopifnot(mediator_name == mrna_name)
  }
} else {
  stopifnot(mediator_name != "")
  mrna_id <- 
    (mrna$annot %>%
       filter(symbol == mediator_name))$id
  stopifnot(length(mrna_id) == 1)
}

mrna_pos <- 
  (mrna$peaks %>%
     filter(gene_id == mrna_id, 
            qtl_chr == chr_id))$qtl_pos
mrna.expr <- as.matrix(mrna$expr[, mrna_id, drop = FALSE])
mrna.annot <- 
  dplyr::filter(
    mrna$annot,
    id == mrna_id)
rm(mrna)

# Association mapping for causal models

# Reduce to common IDs
m <- qtl2::get_common_ids(genoprobs,
                          target,
                          addcovar,
                          mrna.expr,
                          complete.cases = TRUE)
genoprobs <- subset(genoprobs, ind = m)
target <- target[m,, drop = FALSE]
mediator <- mrna.expr[m,, drop = FALSE]
kinship <- qtl2::decomp_kinship(kinship[m,m])
addcovar <- addcovar[m,, drop=FALSE]

t.m <- qtl2::fit1(cbind(1, as.matrix(mediator)),
                  target,
                  kinship,
                  addcovar)
nind <- length(t.m$ind_lod)
t.m <- t.m$lod
m.t <- qtl2::fit1(cbind(1, as.matrix(target)),
                  mediator,
                  kinship,
                  addcovar)$lod

assoc_tar = qtl2::scan1snps(genoprobs = genoprobs[,chr_id], map = map, 
                      pheno = target,
                      kinship = kinship,
                      addcovar = addcovar, chr = chr_id, start = start,
                      end = end, query_func = query_variant, cores = 4,
                      keep_all_snps = FALSE)
assoc_med = qtl2::scan1snps(genoprobs = genoprobs[,chr_id], map = map, 
                      pheno = mediator,
                      kinship = kinship,
                      addcovar = addcovar, chr = chr_id, start = start,
                      end = end, query_func = query_variant, cores = 4,
                      keep_all_snps = FALSE)
assoc_tar_med = qtl2::scan1snps(genoprobs = genoprobs[,chr_id], map = map, 
                      pheno = target,
                      kinship = kinship,
                      addcovar = cbind(addcovar, mediator), chr = chr_id, start = start,
                      end = end, query_func = query_variant, cores = 4,
                      keep_all_snps = FALSE)
assoc_cmst <- assoc_tar
assoc_cmst$lod <- cbind(assoc_med$lod, assoc_tar$lod, assoc_tar$lod,
                        assoc_tar_med$lod, assoc_tar$lod)
colnames(assoc_cmst$lod) <- c("causal","reactive","independent","undecided","mediation")
assoc_cmst$lod[,1] <- t.m + assoc_med$lod[,1]
assoc_cmst$lod[,2] <- m.t + assoc_tar$lod[,1]
assoc_cmst$lod[,3] <- assoc_tar$lod[,1] + assoc_med$lod[,1] 
assoc_cmst$lod[,4] <- assoc_tar_med$lod[,1] + t.m + assoc_med$lod[,1] 
assoc_cmst$lod[,5] <- assoc_tar$lod[,1] - assoc_tar_med$lod[,1]

rm(assoc_tar, assoc_med, assoc_tar_med)
invisible(gc(verbose = FALSE))

nmar_snp <- 3
penalty_snp <- c(nmar_snp, nmar_snp, 2*(nmar_snp-1), 2*nmar_snp-1,0) * log10(nind) / 2

tmp <- assoc_cmst$lod -
  matrix(penalty_snp, nrow(assoc_cmst$lod), ncol(assoc_cmst$lod), byrow = TRUE) +
  penalty_snp[1]

ts <- qtl2pattern::top_snps_all(tmp, assoc_cmst$snpinfo) %>%
   group_by(pheno) %>%
   summarize(minpos = min(pos[lod == max(lod)]),
             maxpos = max(pos[lod == max(lod)]),
             snp_id = snp_id[which.max(lod)][1],
             sdp = paste(sort(unique(sdp[lod == max(lod)], collapse = ","))),
             lod = max(lod)) %>%
  arrange(desc(lod))

## Mediation test with best SNPs

sp <- qtl2::genoprob_to_snpprob(genoprobs, assoc_cmst$snpinfo)[[chr_id]][,, ts$snp_id[1]]

tar_test <- intermediate::mediation_test(
  target = target,
  mediator = mediator,
  annotation = mrna.annot,
  covar_tar = covar,
  covar_med = covar,
  kinship = kinship,
  driver = sp)

tar_best <- summary(tar_test) %>%
  select(symbol,chr,pvalue,triad,mediation,LR) %>%
  rename(lod = "LR") %>%
  mutate(lod = lod / log(10),
         mediation = mediation / log(10))

# Mediation test over interval

sp <- qtl2::genoprob_to_snpprob(genoprobs, assoc_cmst$snpinfo)[[chr_id]]

mediators <- mediator[, rep(mrna_id, dim(sp)[3]), drop = FALSE]
colnames(mediators) <- dimnames(sp)[[3]]

m <- match(mrna_id, mrna.annot$id)
annotation <- mrna.annot[rep(m, ncol(mediators)),]
annotation$id <- colnames(mediators)
annotation$driver <- colnames(mediators)
annotation$pos <- assoc_cmst$snpinfo$pos
annotation$qtl_pos <- NULL
annotation$qtl_lod <- NULL

# Within interval already selected, 
#   pick the region of maximum joint LOD (within 1.5 of peak).
#   run mediation test and find the best models (using BIC among the four models)
#   report out pvalues.

m <- which(assoc_cmst$lod[,"undecided"] >= max(assoc_cmst$lod[,"undecided"]) - 1.5)
rng <- range(assoc_cmst$snpinfo$pos[m])
m <- which(annotation$pos >= rng[1] & annotation$pos <= rng[2])

med_test_scan <- intermediate::mediation_test(
    target = target,
    mediator = mediators[,m],
    annotation = annotation[m,],
    covar_tar = covar,
    covar_med = covar,
    kinship = kinship,
    driver = NULL,
    driver_med = sp[,,m])

best_snp <-
  (med_test_scan$test %>%
     filter(model == "undecided",
            LR >= max(LR) - 1.5 * log(10)))$id
tmp <- med_test_scan$best
m <- match(best_snp, tmp$id)
best_snp <- tmp[m,]
m <- match(best_snp$id, assoc_cmst$snpinfo$snp_id)
best_snp$sdp <- assoc_cmst$snpinfo$sdp[m]
best_snp <- 
  best_snp %>%
  arrange(pvalue) %>%
  mutate(pattern = qtl2pattern::sdp_to_pattern(sdp, LETTERS[1:8]))

topsnps <- 
  qtl2::index_snps(map,
    query_variant(chr_id, start, end)) %>%
  filter(pos >= min(best_snp$pos),
         pos <= max(best_snp$pos),
         sdp %in% unique(best_snp$sdp))
m <- match(best_snp$id, topsnps$snp_id)
best_snp$index <- topsnps$index[m]
topsnps <- left_join(topsnps,
                     best_snp %>%
                       select(-lod) %>%
                       rename(lod = "LR") %>%
                       mutate(lod = lod / log(10)) %>%
                       select(index, triad, lod, pvalue, symbol),
                     by = "index")

if(resultpath != "") {
  write.csv(
    topsnps, 
    file = file.path(
      resultpath, 
      paste0(target_name, "_", chr_id, "_", mediator_name, "_topsnps.csv")))
}