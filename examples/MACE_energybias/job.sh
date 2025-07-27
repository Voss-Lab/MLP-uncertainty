#!/bin/bash
#SBATCH -t 00:15:00
#SBATCH -J mace
#SBATCH -e err
#SBATCH --partition=xyz
#SBATCH --account=xyz
#SBATCH -N 1
#SBATCH -c 12
#SBATCH --ntasks-per-node=1
#SBATCH --gpus-per-task=1
#SBATCH --gpu-bind=single:1
##SBATCH --gpus=4
##SBATCH --gres=gpu:4
#SBATCH -q xyz

cd $SLURM_SUBMIT_DIR
mpirun python3.11 dyn.py > out_md
