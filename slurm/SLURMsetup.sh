#!/bin/bash
# slurm/SLURMsetup.sh

# Set up data in /workspace2/yandell/data
# Run R on R/SLURMsetup.R to set up R libraries

## Set up local data.
cd ~/Documents/Research/rqtl/qtl2rmd
ssh yandell@lunchbox.stat.wisc.edu mkdir /workspace2/yandell/data
scp data/*.{csv,rds,txt} yandell@lunchbox.stat.wisc.edu:/workspace2/yandell/data
## Get CCmouse stuff (this takes awhile)
ssh yandell@lunchbox.stat.wisc.edu mkdir /workspace2/yandell/data/CCmouse
cd ~/Documents/Research/qtl2shiny/qtl2shinyData/CCmouse
scp *.{rds,sqlite} yandell@lunchbox.stat.wisc.edu:/workspace2/yandell/data/CCmouse
scp -r AttieDOv2 yandell@lunchbox.stat.wisc.edu:/workspace2/yandell/data/CCmouse

## Set up R scripts
cd ~/Documents/Research/rqtl/qtl2rmd/R
ssh yandell@lunchbox.stat.wisc.edu mkdir /workspace2/yandell/R
scp *.R yandell@lunchbox.stat.wisc.edu:/workspace2/yandell/R

## Set up SLURM and BASH scripts
## Some of these are slurm scripts; some are bash to run slurm scripts
## See for instance medateQtl2.sh for SLURM setup
cd ~/Documents/Research/rqtl/qtl2rmd/slurm
ssh yandell@lunchbox.stat.wisc.edu mkdir /workspace2/yandell/slurm
scp *.sh yandell@lunchbox.stat.wisc.edu:/workspace2/yandell/slurm



