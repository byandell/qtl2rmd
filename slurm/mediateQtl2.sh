#!/bin/bash
#SBATCH --mail-type=ALL # send us an email about all events
#SBATCH --mail-user=brian.yandell@wisc.edu
# SBATCH this is a comment
#SBATCH -J test_mediateQtl2 # give the job a name
#SBATCH -t 1:00:00 # run this simulation for 5 days (HOURS:MINUTES:SECONDS)
#SBATCH -p short # send this job to the long partition since we want it to run more than 4 days
#SBATCH --mem-per-cpu=500M
#SBATCH --cpus-per-task=4
#SBATCH -n 1
# Load R version 3.4.4
module load R/R-3.4.4
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export R_LIBS_USER=/workspace2/yandell/Rlib
# Run our R program myjob.R which resides in /workspace/[username]/
cd /workspace2/yandell
R CMD BATCH --no-save R/mediateQtl2.R test_mediateQtl2.out
