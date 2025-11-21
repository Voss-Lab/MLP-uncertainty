import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
from sklearn.metrics import r2_score
import os

plt.rcParams['font.sans-serif'] = ['Arial', 'Liberation Sans', 'DejaVu Sans']
plt.rcParams['font.family'] = 'sans-serif'

# Ensemble directories
ensemble_dirs = ['2', '3', '4', '5']

# Create figure with 4 subplots in a row
fig, axes = plt.subplots(1, 4, figsize=(18, 5), sharex=True, sharey=True)

for i, d in enumerate(ensemble_dirs):
    # Load data
    sigmaE = np.loadtxt(os.path.join(d, 'std'))
    error = np.loadtxt(os.path.join(d, 'err'))

    # Compute metrics
    r2 = r2_score(error, sigmaE)
    calibration_error = np.mean(np.abs(sigmaE - error))

    # Hexbin plot
    hb = axes[i].hexbin(
        sigmaE, error,
        gridsize=50,
        cmap='Blues',
        mincnt=1,
        norm=LogNorm()
    )

    # Diagonal line
    axes[i].plot([0, 0.02], [0, 0.02], 'k--', linewidth=0.8)

    # Axis limits and ticks
    axes[i].set_xlim(0, 0.02)
    axes[i].set_ylim(0, 0.02)

    # Set xticks with reduced overlap
    axes[i].set_xticks(np.arange(0, 0.021, 0.01))  # larger spacing
    axes[i].set_yticks(np.arange(0, 0.021, 0.01))
    axes[i].tick_params(labelsize=16)

    # Add title inside plot (top-left)
    axes[i].text(0.001, 0.019, f'{d} members', fontsize=18, ha='left', va='top')

    # Add metrics text just below the title, left-aligned
    axes[i].text(0.001, 0.016, f'R²={r2:.2f}\nMAE={calibration_error:.3f}',
                 fontsize=14, ha='left', va='top',
                 bbox=dict(facecolor='white', alpha=0.5, edgecolor='none'))

# Common axis labels
fig.text(0.5, 0.04, r'$\sigma_E$ (eV/atom)', ha='center', fontsize=18)
fig.text(0.04, 0.5, 'Error (eV/atom)', va='center', rotation='vertical', fontsize=18)

# Colorbar (shared)
cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])
cbar = fig.colorbar(hb, cax=cbar_ax)
cbar.set_label('Counts (log scale)', fontsize=18)
cbar.ax.tick_params(labelsize=16)

# Adjust spacing to prevent overlaps
plt.subplots_adjust(left=0.1, right=0.9, bottom=0.15, top=0.9, wspace=0.20)

# Save and show
plt.savefig('hexbin.png', dpi=300)
plt.show()

