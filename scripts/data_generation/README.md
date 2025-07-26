# Data Generation Scripts

This directory contains scripts and input files used to generate datasets for training and benchmarking the machine learning potential (MLP).

# Folder Descriptions

1. **`reduced_scale_model/`**  
   Initial structure used to kickstart active learning, representing a simplified version of the benchmark system.

2. **`structured_seed/`**  
   Scripts to generate the *structured seed* dataset referenced in **Table S1**. These are curated, motif-rich structures.

3. **`unstructured_seed/`**  
   Scripts to generate the *unstructured seed* dataset also referenced in **Table S1**. These include more random and diverse atomic arrangements in slab structures.

4. **`target_interface/`**  
   structure file for the benchmark system (target interface), used as a reference to generate data and evaluate extrapolation performance of the MLP.

5. **`utilities/`**  
   Miscellaneous helper scripts and tools used across the data generation process.

