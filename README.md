# qtl2rmd
Rmarkdown for qtl2 mediation analyses

The directory has a collection of scripts, some for one-off and others for batch production runs on the [SLURM](https://slurm.schedmd.com/) workload manager.

### SLURM folder

- SLURMsetup.sh (run this locally to migrate files)
- mediation runs
    + mediateRuns.sh
    + mediateQtl2.sh
- simulation runs
    + slurmSim.sh
    + mediateSims.sh
    + SimsSummary.sh

### Rmd folder

Rmd are for exploration and reporting (see summary scripts).

- mediateOne.Rmd (near production for one mediator)
- mediateQtl2.Rmd (production for one target, all mediators)
- mediateSummary.Rmd (summary of mediation for multiple targets)
- mediateSimmary.Rmd (summary of mediation simulations)
- mediateTest.Rmd (local use; expore mediator choices)
- analyzeQtl2.Rmd (local use; exploratory)

### R folder

These are for use in Rmd and in SLURM.

SLURM setup and Rscripts

- SLURMsetup.R (run this remotely to install packages)
- mediateRuns.R (production runs across targets and mediators; see mediateRuns.sh)
- mediateSims.R (production runs for simulations; see slurmSim.sh & mediateSims.sh)
- SimsSummary.R (initial summary at submit node to aggregate simulation data)

Data import and setup

- dataJax.R
- dataJaxMadison.R
- dataMadison.R
- dataPeaks.R

Local batch runs to create multiple Rmd files.

- batchOne.R
- batchTarget.R
- runAnalysis.R
- runMediation.R
- runOne.R

Various utilities

- extract_tables.R
- med_comp.R
- mediateQtl2.R (this is the key routine)
- mediator_quants.R
- qtl2_render.R
- qtl2_rmd.R
- shuffleQtl2.R

This is built off a script developed by Dan Gatti. In Rstudio, use `Knit with Parameters` from <kbd>Knit</kbd> script pulldown menu to adjust parameter settings. this can also be done with the following commands (assuming you are in the parent folder of qtl2rmd):

The [analyzeQtl2.Rmd](analyzeQtl2.Rmd) Rmarkdown file depends on the following packages and their dependencies:

- [tidyverse](https://www.tidyverse.org/)
- [plotly](https://plot.ly/)
- [AnnotationHub](https://bioconductor.org/packages/release/bioc/html/AnnotationHub.html) (if using [dataJax.R](dataJax.R) data setup)
- [qtl2](https://github.com/rqtl/qtl2)
- [qtl2ggplot](https://github.com/byandell/qtl2ggplot)
- [qtl2pattern](https://github.com/byandell/qtl2pattern)
- [intermediate](https://github.com/byandell/intermediate)

The latter three packages were developed by Yandell in conjunction with key people from The Jackson Lab and UW-Madison, as indicate in their description files. See their readme files for installation instructions.

### Data Setup

Three data setup files are provided

- [dataMadison.R](dataMadison.R): uses UW-Madison data setup prepared for [qtl2shiny](http://www.stat.wisc.edu/~yandell/software/qtl2shiny)
- [dataJaxMadison.R](dataJaxMadison.R): adapted from `dataJax.R` from Dan Gatti
- [dataJax.R](dataJax.R): extracted from code from Dan Gatti (should work at Jax)

The purpose of these files is to separate data setup from data analysis. Files are used to create or extract the following:

- `pheno_data`: phenotype data
- `covar`: covariates
- `kinship`: kinship object
- query routines
    + `query_variant(chr, start, stop)`
    + `query_gene(chr, start, stop)`
    + `query_probs(chr, start, stop)`
    + `query_mrna(chr, start, stop)`

Notice the query routines from [qtl2](https://github.com/rqtl/qtl2) as well as additional query routines for genotype probabilities and mRNA data. These query routines isolate the way these objects are extracted from databases. In particular, [dataMadison.R](dataMadison.R) uses some routines from [qtl2pattern](https://github.com/byandell/qtl2pattern) to get data stored in [fst](http://www.fstpackage.org/), while [dataJax.R](dataJax.R) uses and `Rdata` object and the [AnnotationHub](https://bioconductor.org/packages/release/bioc/html/AnnotationHub.html) package.

### Paramaters

Many parameters are used to set up customized analysis. Most notable is the `pheno_name`, but also the `dataSetup`. See the YAML at top of [analyzeQtl2.Rmd](analyzeQtl2.Rmd), or use `Knit with Parameters` from <kbd>Knit</kbd> script pulldown menu to see and adjust parameter settings. Details are in the Rmarkdown and data setup files. 

