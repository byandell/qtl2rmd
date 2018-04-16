## Set upAttie data

## This assumes below that covariates and data transform are fixed.
## See dataMadison.R for alternative

# Load Attie data.
load("D:/Attie_Islet_Secretion/Attie_islet_secr_data_v2.Rdata")

pheno_data <- pheno_clin

covar = model.matrix(~sex + DOwave, data = pheno_data)[,-1]
pheno_data[, pheno_name] = log(pheno_data[, pheno_name])

#############################################################################

## Query functions

## Set up SNP variant query functions.

cc_dbfile <- "C:/Users/dgatti/Documents/Sanger/cc_variants.sqlite"
query_variant <- qtl2::create_variant_query_func(cc_dbfile, filter = "type=='snp'")

## Set up gene query functions.

## MGI query

gene_dbfile <- "C:/Users/dgatti/Documents/Sanger/mouse_genes.sqlite"
query_gene_MGI <- qtl2::create_gene_query_func(dbfile = gene_dbfile,
                                                 filter="(source=='MGI')")

## AnotationHub query
create_gene_query_func_AH <- function(pattern = c("ensembl", "gtf", "mus musculus", "90"),
                                      filename = "Mus_musculus.GRCm38.88.gtf",
                                      chr_field = "seqnames",
                                      start_field = "start", stop_field = "end",
                                      filter = NULL) {
  if(is.null(pattern) & is.null(filename))
      stop("must provide pattern and filename")
  
  # Visit AnnotationHub to get ensemble entries.
  hub <- AnnotationHub::AnnotationHub()
  hub <- AnnotationHub::query(hub, pattern)
  ensembl <- hub[[names(hub)[hub$title == filename]]]
  ensembl <- as.data.frame(ensembl[ensembl$type == "gene"])
  colnames(ensembl)[colnames(ensembl) == "gene_id"] <- "ensembl_gene"
  colnames(ensembl)[colnames(ensembl) == "gene_name"] <- "symbol"
  ensembl[[start_field]] <- ensembl[[start_field]] * 10^-6
  ensembl[[stop_field]]  <- ensembl[[stop_field]] * 10^-6
  ensembl <- ensembl[, sapply(ensembl, function(x) !all(is.na(x)))]

  function(chr, start = NULL, end = NULL) {
    subset_ensembl <- ensembl[[chr_field]] == chr
    if(!is.null(start))
      subset_ensembl <- subset_ensembl & (ensembl[[start_field]] >= start)
    if(!is.null(end))
      subset_ensembl <- subset_ensembl & (ensembl[[stop_field]] >= end)
    ensembl[subset_ensembl,]
  }
}
query_gene_AH <- create_gene_query_func_AH()

## Combined query

query_gene <- function(chr, start = NULL, end = NULL) {
  MGI <- query_gene_MGI(chr, start, end)
  AH <- query_gene_AH(chr, start, end)
  ts = ts %>% 
  dplyr::select(
    dplyr::left_join(
      MGI,
      dplyr::select(
        AH, 
        ensembl_gene, symbol)),
    snp_id, chr, pos, alleles, ensembl_gene, symbol, consequence:lod)
}

## List environment for query_gene

ls.str(environment(query_gene))

## Query genoprobs

query_probs <- function(chr = NULL, start = NULL, end = NULL) {
  if(is.null(chr))
    return(list(probs = genoprobs, map = map))
  
  probs <- genoprobs[, chr]
  map <- map[chr]
  
  if(length(chr) > 1)
    return(list(probs = probs, map = map))
  
  keep <- rep(TRUE, length(map[[1]]))
  if(!is.null(start))
    keep <- keep[map[[1]] < start] <- FALSE
  if(!is.null(end))
    keep <- keep[map[[1]] > end] <- FALSE
  
  list(probs = subset(probs, mar = names(map[[1]][keep])),
       map = map[[1]][keep])
}

## query mRNA

query_mrna <- function(chr = NULL, start = NULL, end = NULL) {
  if(is.null(chr)) {
    cat("need to supply chr\n", file = stderr())
    return(NULL)
  }
  annot <- annot.mrna[annot.mrna$chr == chr, ]
  list(expr = expr.mrna[, annot$id],
       annot = annot)
}
