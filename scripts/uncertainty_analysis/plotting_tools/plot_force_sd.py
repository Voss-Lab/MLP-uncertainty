import h5py
import matplotlib.pyplot as plt
import numpy as np
from ase import neighborlist
from ase.io import Trajectory

# ------------------------- User Inputs -------------------------
traj_file = 'run2.traj'
force_h5  = 'force_sd.h5'
coord_h5  = 'coord.h5'
atoms_to_plot = [1206, 1047]
dt_fs = 0.5  # timestep in fs

site_colors = {
    'gas H': '#1f77b4',
    'gas H\u2099': 'red',
    'H*': 'teal',
    'H\u2099*': 'darkred'
}
# ---------------------------------------------------------------

# Load force SD and coordination
with h5py.File(force_h5, 'r') as f:
    force_sd = f['forcenormstd'][()]

with h5py.File(coord_h5, 'r') as f:
    coord_all = f['coord'][()]

nsteps, natoms = force_sd.shape
time_ps = np.arange(nsteps) * dt_fs / 1000

# Load trajectory for adsorption type
traj = Trajectory(traj_file)
slab_indices = [i for i, atom in enumerate(traj[-1]) if atom.symbol == 'Pt']
ads_indices  = [i for i, atom in enumerate(traj[-1]) if atom.symbol == 'H']

def classify_adsorption(atoms, Ai, slab_indices, all_ads_indices):
    cutoffs = neighborlist.natural_cutoffs(atoms)
    nl = neighborlist.NeighborList(cutoffs, self_interaction=False, bothways=True)
    nl.update(atoms)
    conn = nl.get_connectivity_matrix(sparse=False)
    neighbors = np.nonzero(conn[Ai])[0]
    in_ads = set(neighbors).intersection(all_ads_indices)
    in_slab = set(neighbors).intersection(slab_indices)
    if len(in_ads) == 0:
        return 'gas H' if len(in_slab)==0 else 'H*'
    else:
        return 'gas H\u2099' if len(in_slab)==0 else 'H\u2099*'

# Compute adsorption series
ads_series = {i: [] for i in atoms_to_plot}
for step, atoms in enumerate(traj):
    for i in atoms_to_plot:
        ads_series[i].append(classify_adsorption(atoms, i, slab_indices, ads_indices))

# Plotting function for force SD
def plot_force_sd(time_ps, sd_array, coord_array, ads_array, sd_label, filename):
    fig, ax = plt.subplots(figsize=(6,3))
    line, = ax.plot(time_ps, sd_array, color='black', label=sd_label, linewidth=1.5, zorder=2)
    ax.set_xlabel('Time (ps)', fontsize=16)
    ax.set_ylabel('$\sigma F$ (eV/Å)', fontsize=16)
    ax.tick_params(axis='both', labelsize=12)

    # Coordination on right-hand axis
    ax2 = ax.twinx()
    ax2.set_ylabel('Coordination', fontsize=16)
    ax2.set_ylim(-0.5, 12)  # fixed for visibility
    ax2.tick_params(axis='both', labelsize=12)

    scatter_handles = []
    for site_type, color in site_colors.items():
        mask = np.array(ads_array) == site_type
        if np.any(mask):
            sc = ax2.scatter(time_ps[mask], coord_array[mask],
                             c=color, label=site_type, alpha=0.8,
                             s=25, edgecolors='none', zorder=3)
            scatter_handles.append(sc)

    # Combined legend
    handles, labels = [line], [sd_label]
    scat_handles, scat_labels = ax2.get_legend_handles_labels()
    handles += scat_handles
    labels += scat_labels
    ax.legend(handles, labels, fontsize=12, loc='upper left', scatterpoints=1, markerscale=1.5)

    plt.tight_layout()
    plt.savefig(filename, dpi=300)
    plt.close()

# Generate force SD plots
for i in atoms_to_plot:
    plot_force_sd(time_ps, force_sd[:, i], coord_all[:, i], ads_series[i],
                  '$\sigma F$', f'plotf_atom_{i}.png')

print("Force SD plots with coordination overlay saved!")

