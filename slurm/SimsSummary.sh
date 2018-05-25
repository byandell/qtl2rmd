#!/bin/bash
#SBATCH --mail-type=ALL # send us an email about all events
#SBATCH --mail-user=brian.yandell@wisc.edu
#SBATCH -p short
#SBATCH -t 1:00:00
#SBATCH --array=0-59
#SBATCH --mem-per-cpu=500M
#SBATCH --cpus-per-task=4
#SBATCH -n 1
# Load R version 3.4.4
module load R/R-3.4.4
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export R_LIBS_USER=/workspace2/yandell/Rlib
ROW=`expr $SLURM_ARRAY_TASK_ID + 2`
export TARGET=`awk -v line=$ROW '{if(NR==line)print $1}' data/exvivo_peak_summary.txt`
export CHR=`awk -v line=$ROW '{if(NR==line)print $2}' data/exvivo_peak_summary.txt`
Rscript R/SimsSummary.R $TARGET $CHR