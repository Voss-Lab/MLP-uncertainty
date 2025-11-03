# Local Uncertainty Analysis

Tools for analyzing uncertainty in MPNN ensemble simulated dynamics.

## Subdirectories

- `configuration_space_sampling`
  - `dyn.py`: Runs MD simulations using the MPNN ensemble potential.
  - `spike.py`, `split_spike.py`: Identify high-uncertainty atom-step combinations from MD trajectories, based on node energies.

- `energy`
  - `atom_step_plot.py`: Plots the first occurrence of threshold crossing across the three bonding regimes, based on node energies.
  - `sigma_Enode.py`: Generates Figure 2, showing node energy evolution in the sampled dynamics.
  - `plot_node_sd.py`: Plots node energy standard deviation over time for selected atoms, overlaid with coordination numbers and adsorption states (Figure 3).
  - `violin.py`: Generates Figure S4, visualizing distributions of node energy SD.

- `force`
  - `plot_force_sd.py`: Plots force standard deviation over time for selected atoms, overlaid with coordination numbers and adsorption states (Figure 4a).
  - `map_local_uncertainty.py`: Computes 3σ thresholds for node energy SD and force norm SD, identifies first crossings for atoms exceeding thresholds, and generates a summary text file and a cluster plot (Figure 4b).
  - Additional files (`force_sd.h5`, `md_data.h5`, `err`, `first_crossings.txt`,`cluster_species.png`) shows the processed data and intermediate results.

- `structural_analysis`
  - `get_coord.py`: Computes atomic coordination numbers from a trajectory and saves them in an HDF5 file for further analysis or plotting.
  - `stats.py`: Tracks local environment descriptors (local coordination, bonding regime, molecular state) of atoms along the MD trajectory.
  - `local_cluster.py`: Extracts clusters of atoms around a central atom based on first threshold crossing data for forces, and computes PBC-corrected distances from the central atom.
