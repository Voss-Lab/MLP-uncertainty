"""
The script reads MD simualtion data from an HDF5 file to identify atom-steps where the local uncertainty exceeds a specified threshold (spikes)
It prints all atom-step combinations, all unique steps, total count, and the first 100 unique steps where spikes occurred
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
all_occurrences = []

for step in range(steps):
    for atom in range(atoms):
        if node_sd[step, atom] > threshold:
            all_occurrences.append((atom, step))          

# Extract unique steps from all occurrences
unique_steps = sorted({step for _, step in all_occurrences})
unique_step_count = len(unique_steps)

# Print results
print("Occurrences of spikes (atom, step):", all_occurrences)
print("Unique steps where spikes occurred:", unique_steps)
print("Total count of unique steps:", unique_step_count)
print("First 100 unique steps where spikes occurred:", unique_steps[:100])