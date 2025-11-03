"""
This script creates bootstrap samples in extended XYZ format for MLP training in MACE,
by randomly selecting structures with replacement from the full training dataset.
"""

import random
from ase.io import read, write

xyz_file = 'train.xyz'
num_structures = 3000 # number of structures in the training dataset
num_samples = 10
molecules = read(xyz_file, index=':')
if len(molecules) != num_structures:
    raise ValueError(f"The XYZ file contains {len(molecules)} structures, expected {num_structures}.")

# Create specified number of bootstrap samples
for i in range(num_samples):
    sampled_molecules = [random.choice(molecules) for _ in range(num_structures)]
    output_file = f'bag_sample_{i+1}.xyz'
    write(output_file, sampled_molecules)