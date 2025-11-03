# Uncertainty Analysis

Tools for analyzing uncertainty in MPNN ensemble predictions and visualizing results from active learning simulations.

## Subdirectories

- `configuration_space_sampling`
  - `dyn.py`: Runs MD simulations using the MPNN ensemble potential.
  - `spike.py`, `split_spike.py`: Identify high-uncertainty atom-step combinations from MD trajectories.

- `plotting_tools`
  - `atom_step_plot.py`: plots the first occurence of threshold crossing across the three bonding regimes
  - `sigma_Enode.py`: Generates Figure 2
  - `stats.py`: Tracks the local environment descriptors (local coordination, bonding regime, molecular state) of atoms along the MD trajectory
  - `violin.py`: Generates Figure S4.
  - `get_coord.py`: Computes atomic coordination numbers from a trajectory and saves them in an HDF5 file for further analysis or plotting.
  - `plot_node_sd.py`: Plots node energy standard deviation over time for selected atoms, overlaid with coordination numbers and adsorption states (Figure 3).
  - `plot_force_sd.py`: Plots force standard deviation over time for selected atoms, overlaid with coordination numbers and adsorption states (Figure 4a).

- `energy_versus_force`
  - `map_local_uncertainty.py`: Computes 3σ thresholds for node energy SD and force norm SD, identifies first crossings for atoms exceeding thresholds, and generates a summary text file and a cluster plot (Figure 4b).
  - `local_cluster.py`: Extracts clusters of atoms around a central atom based on first crossing data, and computes PBC-corrected distances from the central atom.

