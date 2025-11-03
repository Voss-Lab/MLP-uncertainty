"""
This script reads MD simulation data from an HDF5 file and calculates the average threshold for local uncertainty exceedance (spikes) 
for atoms in different bonding regimes: surface Pt, intermediate Pt, surface H, and gas phase H. 
It identifies the first occurrence of spikes for each atom in these ranges and prints the results.
"""
import h5py
import numpy as np
from ase.io import read, write

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

# Create dictionaries to store first occurrences for each range
first_occurrence_range_1 = {}
first_occurrence_range_2 = {}
first_occurrence_range_3 = {}
first_occurrence_range_4 = {}

steps = node_sd.shape[0]
atoms = node_sd.shape[1]

# Loop through steps and atoms to check for threshold exceedance
for step in range(steps):
    for atom in range(atoms):
        # Check for each atom's threshold exceedance in the corresponding range
        if atom in range_1:
            threshold = calculate_average_threshold(range_1)
            if node_sd[step, atom] > threshold and atom not in first_occurrence_range_1:
                first_occurrence_range_1[atom] = step
        elif atom in range_2:
            threshold = calculate_average_threshold(range_2)
            if node_sd[step, atom] > threshold and atom not in first_occurrence_range_2:
                first_occurrence_range_2[atom] = step
        elif atom in range_3:
            threshold = calculate_average_threshold(range_3)
            if node_sd[step, atom] > threshold and atom not in first_occurrence_range_3:
                first_occurrence_range_3[atom] = step
        elif atom in range_4:
            threshold = calculate_average_threshold(range_4)
            if node_sd[step, atom] > threshold and atom not in first_occurrence_range_4:
                first_occurrence_range_4[atom] = step

# Calculate and print the average threshold for each range
print("\nAverage threshold for surface Pt atoms (range 0-15, 48-63):")
print(calculate_average_threshold(range_1))

print("\nAverage threshold for intermediate Pt atoms (range 16-47):")
print(calculate_average_threshold(range_2))

print("\nAverage threshold for surface H (range 64-95):")
print(calculate_average_threshold(range_3))

print("\nAverage threshold for gas phase H (range 96-143):")
print(calculate_average_threshold(range_4))

# Print the first occurrence for each range
print("\nFirst occurrence of threshold exceedance for surface Pt atoms (0-15, 48-63):")
print(sorted(first_occurrence_range_1.items()))

print("\nFirst occurrence of threshold exceedance for intermediate Pt atoms (16-47):")
print(sorted(first_occurrence_range_2.items()))

print("\nFirst occurrence of threshold exceedance for surface H (64-95):")
print(sorted(first_occurrence_range_3.items()))

print("\nFirst occurrence of threshold exceedance for gas phase H (96-143):")
print(sorted(first_occurrence_range_4.items()))