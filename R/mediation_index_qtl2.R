mediation_index_qtl2 <- function(target, mediator,
                             annotation, covar_tar, covar_med, kinship,
                             genoprobs, map,
                             drop_lod = 1.5, query_variant) {

  ############################################################
  # Association mapping for joint mediator and target
  # Goal: find region of maximal LOD for interval mediation tests.
  
  # Reduce to common IDs
  m <- qtl2::get_common_ids(genoprobs,
                            target,
                            covar_tar,
                            covar_med,
                            mediator,
                            complete.cases = TRUE)
  genoprobs <- subset(genoprobs, ind = m)
  target <- target[m,, drop = FALSE]
  mediator <- mrna.expr[m,, drop = FALSE]
  kinship <- qtl2::decomp_kinship(kinship[m,m])
  covar_tar <- covar_tar[m,, drop=FALSE]
  covar_med <- covar_med[m,, drop=FALSE]
  
  # Regress target on mediator.
  t.m <- qtl2::fit1(cbind(1, as.matrix(mediator)),
                    target,
                    kinship,
                    covar_tar)
  nind <- length(t.m$ind_lod)
  t.m <- t.m$lod

  # Association for mediators.
  assoc_cmst <- qtl2::scan1snps(genoprobs = genoprobs[,chr_id], map = map, 
                              pheno = mediator,
                              kinship = kinship,
                              addcovar = covar_med, chr = chr_id, start = start,
                              end = end, query_func = query_variant, cores = 4,
                              keep_all_snps = FALSE)
  med_lod <- assoc_med$lod[,1]
  # Association for target adjusted by mediators.
  assoc_cmst = qtl2::scan1snps(genoprobs = genoprobs[,chr_id], map = map, 
                               pheno = target,
                               kinship = kinship,
                               addcovar = cbind(covar_tar, mediator),
                               chr = chr_id, start = start, end = end,
                               query_func = query_variant, cores = 4,
                               keep_all_snps = FALSE)
  assoc_cmst$lod[,1] <- assoc_cmst$lod[,1] + t.m + med_lod 

  #########################################################
  # Mediation test on interval
  m <- which(assoc_cmst$lod[,1] >= max(assoc_cmst$lod[,1]) - drop_lod)
  driver_med <- 
    qtl2::genoprob_to_snpprob(
      genoprobs,
      assoc_cmst$snpinfo[m,, drop = FALSE])[[chr_id]]
  driver_index <- assoc_cmst$snpinfo$pos[m]
  
  med_test_scan <- 
    intermediate::mediation_index(target, mediators, driver = NULL,
                            annotation, covar_tar, covar_med, kinship,
                            driver_med, driver_index)
    
  #########################################################
  # Find region of joint (undecided) LR that is within drop_lod of maximum.
  best_snp <-
    (dplyr::filter(
      med_test_scan$test,
      model == "undecided",
      LR >= max(LR) - drop_lod * log(10)))$id
  tmp <- med_test_scan$best
  m <- match(best_snp, tmp$id)
  best_snp <- tmp[m,]
  m <- match(best_snp$id, assoc_cmst$snpinfo$snp_id)
  best_snp$sdp <- assoc_cmst$snpinfo$sdp[m]
  best_snp <- 
    dplyr::mutate(
      dplyr::arrange(
        best_snp, 
        pvalue),
      pattern = qtl2pattern::sdp_to_pattern(sdp, LETTERS[1:8]))
  
  # Find all top SNPs in region that match sdp.
  topsnps <- 
    dplyr::filter(
      qtl2::index_snps(
        map,
        query_variant(chr_id, min(best_snp$pos), max(best_snp$pos))),
      sdp %in% unique(best_snp$sdp))
  
  # Add index column to best_snp for joining.
  m <- match(best_snp$id, topsnps$snp_id)
  best_snp$index <- topsnps$index[m]
  
  dplyr::left_join(
    topsnps,
    dplyr::select(
      dplyr::mutate(
        dplyr::rename(
          dplyr::select(
            best_snp,
            -lod),
          lod = "LR"),
        lod = lod / log(10)),
      index, triad, lod, pvalue, symbol),
    by = "index")
}

