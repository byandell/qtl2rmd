#!/bin/bash
cd /workspace2/yandell
for i in {1..100..1}
do 
export SIM=$i
echo $SIM
sbatch slurm/mediateSims.sh
done