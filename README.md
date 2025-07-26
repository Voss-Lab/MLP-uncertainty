# MLP-Project

This repository serves as the central hub for a machine learning potentials (MLP) project, integrating customized versions of MACE and ASE along with a comprehensive suite of utility scripts for data generation, training, active learning, and analysis. The atomic structures used for training are also included as .xyz files for transparency and reproducibility.

## Submodules

This repository includes two modified forks as Git submodules:
- [`external/mace_node`](https://github.com/sumanbhasker89/mace_node) — custom version of [MACE](https://github.com/ACEsuit/mace)
- [`external/ase_multiproc_calc`](https://github.com/sumanbhasker89/ase_multiproc_calc) — modified [ASE](https://gitlab.com/ase/ase) with multiprocessing calculator support

## Scripts Directory Overview

The `scripts/` folder organizes all supporting code used throughout the MLP workflow:

- `data_generation/` — Scripts and structure files used in model initialization and benchmarking:
- `mlp_training/` — Input and analysis scripts for training machine learning potentials:
- `nn_ensembles/` — Scripts for analyzing uncertainty vs. error correlation:
- `active_learning/` — Scripts to select retraining data based on local uncertainty (spikes):
- `uncertainty_analysis/` — Post-MD tools for ensemble uncertainty analysis:

## Training Data

The `xyz_files/` directory contains the datasets used for training in extended XYZ format. These include the initial seeds and active learning configurations.
