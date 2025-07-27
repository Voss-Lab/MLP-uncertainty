#!/bin/bash
#SBATCH --job-name=test
#SBATCH --error=err
##SBATCH --ntasks=1
#SBATCH --tasks-per-node=120
#SBATCH --cpus-per-task=1
#SBATCH --mem=400G
#SBATCH --time=01:00:00
#SBATCH --partition=xyz
#SBATCH --account=xyz

export OMP_NUM_THREADS=1
python3.8 nve_md.py > out_md
