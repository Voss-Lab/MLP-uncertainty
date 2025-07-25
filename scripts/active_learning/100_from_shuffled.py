"""
This script reads MD simulation data from an HDF5 file and calculates the average threshold for local uncertainty exceedance (spikes) 
for atoms in different bonding regimes: surface Pt, intermediate Pt, surface H, and gas phase H. 
It collects unique steps from all atom-step combinations where spikes occur, and randomly samples 100 frames to include in retraining.
"""
import h5py
import numpy as np
from ase.io import read, write
import random  

with h5py.File('md_data.h5', 'r') as f:
    node_sd = f['node_sd'][()] 
    epot = f['epot'][()]
    sd = f['sd'][()]

# Define function to calculate threshold for a given atom range
def calculate_average_threshold(atom_range):
    mean_node_sd = np.mean(node_sd[:, atom_range], axis=1)
    std_node_sd = np.std(node_sd[:, atom_range], axis=1)
    threshold = mean_node_sd + 3 * std_node_sd
    average_threshold = np.mean(threshold)
    return average_threshold

# Define atom index ranges for 01-Data/active_learning/inp.xyz
range_1 = list(range(0, 16)) + list(range(48, 64))  # Pt in bottom and top surface layers
range_2 = range(16, 48)  # Intermediate Pt 
range_3 = range(64, 96)  # Surface H
range_4 = range(96, 144)  # Gas phase H

# Create dictionaries to store all occurrences for each range
all_occurrences_range_1 = {}
all_occurrences_range_2 = {}
all_occurrences_range_3 = {}
all_occurrences_range_4 = {}

steps = node_sd.shape[0]
atoms = node_sd.shape[1]

# Loop through steps and atoms to check for threshold exceedance
for step in range(steps):
    for atom in range(atoms):
        # Check for each atom's threshold exceedance in the corresponding range
        if atom in range_1:
            threshold = calculate_average_threshold(range_1)
            if node_sd[step, atom] > threshold:
                if atom not in all_occurrences_range_1:
                    all_occurrences_range_1[atom] = []
                all_occurrences_range_1[atom].append(step)
        elif atom in range_2:
            threshold = calculate_average_threshold(range_2)
            if node_sd[step, atom] > threshold:
                if atom not in all_occurrences_range_2:
                    all_occurrences_range_2[atom] = []
                all_occurrences_range_2[atom].append(step)
        elif atom in range_3:
            threshold = calculate_average_threshold(range_3)
            if node_sd[step, atom] > threshold:
                if atom not in all_occurrences_range_3:
                    all_occurrences_range_3[atom] = []
                all_occurrences_range_3[atom].append(step)
        elif atom in range_4:
            threshold = calculate_average_threshold(range_4)
            if node_sd[step, atom] > threshold:
                if atom not in all_occurrences_range_4:
                    all_occurrences_range_4[atom] = []
                all_occurrences_range_4[atom].append(step)

# Collect all unique steps from all occurrences
all_steps = set()
for atom_occurrences in [all_occurrences_range_1, all_occurrences_range_2, all_occurrences_range_3, all_occurrences_range_4]:
    for steps in atom_occurrences.values():
        all_steps.update(steps)

# Shuffle all the steps to randomly sample
all_steps = list(all_steps)  # Convert to list for shuffling
random.shuffle(all_steps)

# Select the first 100 unique steps
sampled_steps = all_steps[:100]
print("\n100 randomly sampled unique steps:")
print(sampled_steps)