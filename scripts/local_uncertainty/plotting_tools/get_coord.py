import h5py
import numpy as np
from ase import neighborlist
from ase.io import Trajectory

traj = Trajectory('run2.traj')
natoms = len(traj[0])
nsteps = len(traj)

coord_all = np.zeros((nsteps, natoms), dtype=int)

for step, atoms in enumerate(traj):
    cutoffs = neighborlist.natural_cutoffs(atoms)
    nl = neighborlist.NeighborList(cutoffs, self_interaction=False, bothways=True)
    nl.update(atoms)
    conn = nl.get_connectivity_matrix(sparse=False)
    # count neighbors for each atom (or only slab neighbors if desired)
    for i in range(natoms):
        coord_all[step, i] = np.sum(conn[i])  # total coordination

# save to HDF5
with h5py.File('coord.h5', 'w') as f:
    f.create_dataset('coord', data=coord_all)

