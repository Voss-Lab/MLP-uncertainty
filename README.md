# MLP-uncertainty

This repository serves as the central hub for a machine learning potentials (MLP) project, integrating customized versions of MACE and ASE along with a comprehensive suite of utility scripts for data generation, training, uncertainty analysis (global, local), and active learning. The atomic structures used for training are also included as .xyz files for transparency and reproducibility.

## Submodules

This repository includes two modified forks as Git submodules:
- [`external/mace_node`](https://github.com/sumanbhasker89/mace_node) custom version of [MACE](https://github.com/ACEsuit/mace)
- [`external/ase_multiproc_calc`](https://github.com/sumanbhasker89/ase_multiproc_calc) modified [ASE](https://gitlab.com/ase/ase) with multiprocessing calculator support

## Scripts Directory Overview

The `scripts/` folder organizes all supporting code used throughout the MLP workflow:

- `data_generation` structure files used in model initialization and benchmarking
- `mlp_training` training and analysis of machine learning potentials
- `global_uncertainty` analysis of uncertainty in ensemble predicted total energy
- `local_uncertainty` Investigations of uncertainty in atomic energy and forces, and associated structural analysis
- `active_learning` retraining based on local energy uncertainty sampling

## Training Data

The `xyz_files` directory contains the datasets used for training in extended XYZ format. These include the initial seeds and active learning configurations.

## Example Simulations

The `examples/` directory contains usage examples of ensemble-based simulations and uncertainty-driven sampling:

- `ASE_multiproc_calc`: Demonstrates a multiprocessing-enabled ASE interface for evaluating ANN ensemble models (e.g., aenet) on CPUs, enabling fast MD with uncertainty estimation.
- `MACE_energybias`: Shows how to perform uncertainty-biased simulations using the MACE ASE calculator by modifying the potential energy with an ensemble uncertainty-based Gaussian bias.

