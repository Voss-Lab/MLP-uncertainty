# Uncertainty Analysis

Tools for analyzing uncertainty in MPNN ensemble predictions and visualizing results from active learning simulations.

## Subdirectories

- `configuration_space_sampling/`  
  - `dyn.py`: Runs MD simulations using the MPNN ensemble potential.  
  - `spike.py`, `split_spike.py`: Identify high-uncertainty atom-step combinations from MD trajectories.

- `plotting_tools/`  
  - `atom_step_plot.py`: Generates Figures 5 and 9.  
  - `sigma_Enode.py`: Generates Figure 6.  
  - `stats.py`: Generates Figures 7 and 8.  
  - `violin.py`: Generates Figure S4.

