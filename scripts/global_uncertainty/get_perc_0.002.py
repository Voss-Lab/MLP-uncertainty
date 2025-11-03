import numpy as np

# Load σE
sigmaE = np.loadtxt('std')

# Given threshold
threshold = 0.002

# Compute the percentile of this threshold
percentile = np.sum(sigmaE <= threshold) / len(sigmaE) * 100
print(f"Threshold {threshold:.6f} eV/atom corresponds to the {percentile:.1f}th percentile")

