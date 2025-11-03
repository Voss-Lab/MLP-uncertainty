# Neural Network Ensemble Tools

This directory contains scripts for creating and analyzing neural network ensemble predictions, including bootstrap sampling, error and uncertainty computation, and visualization.

## Scripts

- `bag_sample.py`: Creates bootstrap samples in extended XYZ format for MLP training in MACE by randomly selecting structures with replacement from the full training dataset. Outputs multiple `.xyz` files representing individual bootstrap samples.

- `get_err_sd.py`: Computes per-atom errors and uncertainties from MPNN ensemble predictions to generate standard deviation (`std`) and error (`err`) files.

- `hexbin.py`: Generates hexbin plots comparing the ensemble-predicted uncertainty (`σ_E`) with the actual error for different ensemble sizes. Computes calibration metrics (R² and mean absolute difference) and saves a multi-panel figure (`hexbin.png`).

