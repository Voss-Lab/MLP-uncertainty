# Examples

This directory contains examples demonstrating the use of customized ASE calculators with machine learning potentials (MLPs) for efficient ensemble simulations and uncertainty-driven sampling.

## Subdirectories

### ASE_multiproc_calc
This example demonstrates a multiprocessing-enabled ASE framework that accelerates MLP ensemble calculations on CPUs. It supports any ML ensemble that can be interfaced as an ASE calculator. The code:
- Evaluates properties like energies and forces using an ensemble of MLPs in parallel
- Computes the ensemble uncertainty (σE) as the standard deviation of predictions
- Enables ensemble-driven MD simulations for configuration space exploration

This framework was used to simulate the aenet ANN ensemble.

### MACE_energybias
This example shows how uncertainty-driven dynamics can be implemented in MACE using an energy bias approach inspired by the UDD-AL method (Kulichenko et al.). The biased potential energy is defined by:

**E_bias = A × exp(−σE² / B²)**

where:
- `A` and `B²` are tunable parameters (`bias_amplitude`, `bias_width`)
- `σE` is the ensemble uncertainty

The `energy_bias` mode in the MACE ASE calculator guides simulations toward high-uncertainty regions, improving the efficiency of active learning.

