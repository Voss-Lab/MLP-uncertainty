"""
The script reads MD simualtion data from an HDF5 file to identify atom-steps where the local uncertainty exceeds a specified threshold (spikes)
It prints the first occurrence of spike for each atom and the unique steps where this happens
"""
import h5py
import numpy as np
from ase.io import read, write

with h5py.File('md_data.h5', 'r') as f:
    node_sd = f['node_sd'][()]
    epot = f['epot'][()]
    sd = f['sd'][()]

# Define threshold for spike detection (e.g., 3 standard deviations from the mean)
mean_node_sd = np.mean(node_sd)
std_node_sd = np.std(node_sd)
threshold = mean_node_sd + 3 * std_node_sd
print('threshold=',threshold)

steps = node_sd.shape[0]
atoms = node_sd.shape[1]
atom_step = []
unique_atoms = set()
first_occurrence = {}

for step in range(steps):
    for atom in range(atoms):
        if node_sd[step, atom] > threshold and atom not in first_occurrence:
            first_occurrence[atom] = step
            atom_step.append((atom, step))
print(f'atom:step for first occurrence of threshold exceedance: {first_occurrence}')

# Extract unique steps from the first_occurrence dictionary
unique_steps = set(first_occurrence.values())
print("Unique steps where threshold first exceeds:", sorted(unique_steps))
