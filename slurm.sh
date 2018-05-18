#!/bin/bash
# /workspace2/yandell/data
scp local/Jax/*.{csv,rds} yandell@lunchbox.stat.wisc.edu:/workspace2/yandell/data
cd ~/Documents/Research/qtl2shiny
scp *{md,R} yandell@lunchbox.stat.wisc.edu:/workspace2/yandell/data/qtl2shiny
cd qtl2shinyData
scp projects.csv yandell@lunchbox.stat.wisc.edu:/workspace2/yandell/data/qtl2shiny/qtl2shinyData
cd CCmouse
scp *.{rds,sqlite} yandell@lunchbox.stat.wisc.edu:/workspace2/yandell/data/qtl2shiny/qtl2shinyData/CCmouse
scp -r AttieDOv2 yandell@lunchbox.stat.wisc.edu:/workspace2/yandell/data/qtl2shiny/qtl2shinyData/CCmouse

# /workspace2/yandell/R
scp dataJaxMadison.R yandell@lunchbox.stat.wisc.edu:/workspace2/yandell/R
scp mediateQtl2.R yandell@lunchbox.stat.wisc.edu:/workspace2/yandell/R
# make changes if needed

# /workspace2/yandell/Rlib
srun --pty /bin/bash
module load R/R-3.4.4
export OPENMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export R_LIBS_USER=/workspace2/yandell/Rlib

# things to do in R using srun on lunchbox
.libPaths(.libPaths()[2:1])
install.packages("fst")
devtools::install_github("rqtl/qtl2", dependencies = FALSE)
devtools::install_github("rqtl/qtl2fst", dependencies = FALSE)
devtools::install_github("byandell/qtl2ggplot", dependencies = FALSE)
devtools::install_github("byandell/qtl2pattern", dependencies = FALSE)
devtools::install_github("byandell/intermediate", dependencies = FALSE)



