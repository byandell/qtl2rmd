## Set upAttie data

## This assumes below that covariates and data transform are fixed.
## See AH_Madison.R for alternative

project_info <- data.frame(project = "AttieDOv2",
                           taxa = "CCmouse",
                           directory = datapath,
                           stringsAsFactors = FALSE)
taxa_dir <- file.path(project_info$directory, 
                      project_info$taxa)
(project_dir <- file.path(taxa_dir, 
                         project_info$project))

analyses <- readRDS(file.path(project_dir, "analyses.rds")) %>%
  filter(pheno == pheno_name)
pheno_data <- readRDS(file.path(project_dir, "pheno_data.rds")) %>%
  qtl2pattern::pheno_trans(analyses$pheno, 
                           analyses$transf,
                           analyses$offset,
                           analyses$winsorize)
covar <- readRDS(file.path(project_dir, "covar.rds")) %>%
  qtl2pattern::get_covar(analyses)

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
