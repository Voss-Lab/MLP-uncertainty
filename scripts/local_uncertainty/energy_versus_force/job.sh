#!/bin/bash
#SBATCH --job-name=analysis
#SBATCH --error=err
#SBATCH --gpus=a100:1
#SBATCH --time=01:00:00
#SBATCH --qos=xyz
#SBATCH --partition=abc
#SBATCH --account=def

python3.11 map_local_uncertainty.py 
