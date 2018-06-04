# R/SLURMsetup.R

# Set up R libraries data in /workspace2/yandell/Rlib
# Run bash on slurm/SLURMsetup.sh to set up data and scripts

## Put local Rlib at from of libPaths
Rlocal <- "/workspace2/yandell/Rlib"
if(!dir.exists(Rlocal)) dir.create(Rlocal)
m <- match(Rlocal, .libPaths())
if(is.na(m)) {
  .libPaths(c(Rlocal, .libPaths()))
} else {
  .libPaths(c(.libPaths()[m], .libPaths()[-m]))
}

## Install packages (order matters).
install.packages("fst")
devtools::install_github("rqtl/qtl2", dependencies = FALSE)
devtools::install_github("rqtl/qtl2fst", dependencies = FALSE)
devtools::install_github("byandell/qtl2ggplot", dependencies = FALSE)
devtools::install_github("byandell/qtl2pattern", dependencies = FALSE)
devtools::install_github("byandell/intermediate", dependencies = FALSE)
