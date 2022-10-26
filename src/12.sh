#!/bin/bash
#SBATCH -J 10.ubq
#SBATCH -o 10.ubq.out
#SBATCH -e 10.ubq.err
#SBATCH -p gpgpu-1
#SBATCH --exclusive
#SBATCH --mem=240GB
#SBATCH --exclude=k001,b001,b032,b160,p003,p025
#SBATCH --array=5-10


for VARIABLE in {2..6}
do
python 12.ubq_training.py ${SLURM_ARRAY_TASK_ID}
done
