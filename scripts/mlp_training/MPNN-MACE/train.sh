#!/bin/bash
#SBATCH -t 01:00:00
#SBATCH -J mace_model
#SBATCH -e err
#SBATCH -o out
#SBATCH --partition=xyz
#SBATCH --account=xyz
#SBATCH -N 1
#SBATCH -c 12
#SBATCH --ntasks-per-node=4
#SBATCH --gpus-per-task=1
#SBATCH --gpu-bind=single:1
#SBATCH -q xyz

cd $SLURM_SUBMIT_DIR
export OMP_NUM_THREADS=1

python3  ../mace/cli/run_train.py \
    --name="mace" \
    --train_file="train.xyz" \
    --test_file="test.xyz" \
    --valid_fraction=0.10 \
    --energy_key="energy" \
    --E0s="average" \
    --forces_key="forces" \
    --model="MACE" \
    --num_channels=128 \
    --num_interactions=2 \
    --max_L=0 \
    --r_max=8.0 \
    --max_num_epochs=1000 \
    --patience=25 \
    --batch_size=50 \
    --valid_batch_size=10 \
    --swa \
    --ema \
    --ema_decay=0.99 \
    --amsgrad \
    --restart_latest \
    --device=cuda \
    --default_dtype='float64' \
    --seed=15 \
