"""
This script reads MD data from an HDF5 file and plots the evolution of local uncertainty for each atom across MD frames 
"""
import h5py
import matplotlib.pyplot as plt

with h5py.File('md_data.h5', 'r') as f:
    node_sd = f['node_sd'][()]
    epot = f['epot'][()]
    sd = f['sd'][()]
steps = node_sd.shape[0]
atoms = node_sd.shape[1]

#adjust the range to correspond to atoms in the respective bonding regimes
for i in range(atoms):
    plt.figure()
    plt.plot(range(steps), node_sd[:, i],color="black")
    plt.title(f'Pt: atom #{i}',fontsize=20)
    #plt.title(f'gas phase $H_2$: atom #{i}',fontsize=20)
    #plt.title(f'H*: atom #{i}',fontsize=20)
    plt.xlabel('MD steps',fontsize=20)
    plt.ylabel('node energy SD, eV',fontsize=20)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.tight_layout()
    plt.savefig(f'plot_atom_{i}.png')  # Save each plot as a PNG file
    plt.close()