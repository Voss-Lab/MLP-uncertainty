import numpy as np
from ase.io import Trajectory, write
from ase import Atoms

# ---------------- User parameters ----------------
traj_file = "run2.traj"                  # ASE trajectory file
first_cross_file = "first_crossings.txt" # text file with crossings

#1206:399, 1047:595

central_atom = 1206   # atom index of node_sd crossing
central_step = 399    # step at which node_sd crossing happens
delta = 10            # +/- step window for nearby force_sd crossings
# -------------------------------------------------

# --- Load trajectory ---
traj = Trajectory(traj_file)
nframes = len(traj)
print(f"Loaded trajectory with {nframes} frames.")

# --- Parse the first_crossings.txt file ---
first_node_cross = {}
first_force_cross = {}

with open(first_cross_file, "r") as f:
    section = None
    for line in f:
        line = line.strip()
        if not line:
            continue
        if line.startswith("First crossing for node_energy_sd"):
            section = "node"
            continue
        elif line.startswith("First crossing for forcenormstd"):
            section = "force"
            continue
        elif line.startswith("Atom"):
            # Defensive parsing: handle "Atom X : Step Y" safely
            try:
                parts = line.replace("Atom", "").replace("Step", "").replace(":", "").split()
                if len(parts) >= 2:
                    atom, step = int(parts[0]), int(parts[1])
                    if section == "node":
                        first_node_cross[atom] = step
                    elif section == "force":
                        first_force_cross[atom] = step
            except ValueError:
                continue  # skip malformed lines
        else:
            continue

# --- Identify cluster atoms within ±delta steps of central_step ---
cluster_atoms = [
    atom for atom, step in first_force_cross.items()
    if abs(step - central_step) <= delta and step > 0
]

print(f"\nCentral atom: {central_atom}")
print(f"Central step: {central_step}")
print(f"Considering ±{delta} steps window: [{central_step - delta}, {central_step + delta}]")
print(f"Number of cluster atoms: {len(cluster_atoms)}")
print(f"Cluster atom indices: {cluster_atoms}")

# --- Extract frame at central_step ---
frame = traj[central_step]
positions = frame.get_positions()
cell = frame.get_cell()                # ASE cell matrix
center_pos = positions[central_atom]

# --- Compute minimum-image (periodic) distances ---
if cluster_atoms:
    cluster_positions = positions[cluster_atoms]

    # Convert to fractional coordinates
    frac_all = positions @ np.linalg.inv(cell)
    frac_center = frac_all[central_atom]
    delta_frac = frac_all[cluster_atoms] - frac_center
    delta_frac -= np.round(delta_frac)          # wrap to [-0.5,0.5] box

    # Convert back to Cartesian
    delta_cart = delta_frac @ cell
    distances = np.linalg.norm(delta_cart, axis=1)

    print("\nMinimum-image distances from central atom to cluster atoms:")
    for a, d in zip(cluster_atoms, distances):
        print(f"Atom {a:4d}  distance = {d:8.4f} Å")

    print(f"\nMean distance: {np.mean(distances):.4f} Å")
    print(f"Max distance:  {np.max(distances):.4f} Å")
else:
    print("\n No force_SD cluster atoms found within ±10 steps of central step.")
    distances = []

# --- Save a visualization snapshot with PBC-corrected positions ---
# Add central atom last for easy highlighting
all_indices = cluster_atoms + [central_atom]
# Positions relative to central atom, wrapped by PBC
cluster_positions_wrapped = np.vstack([delta_cart, np.zeros((1,3))])  # central atom at origin
symbols = [frame.get_chemical_symbols()[i] for i in all_indices]

subset = Atoms(symbols=symbols, positions=cluster_positions_wrapped)
subset.center(vacuum=5.0)

xyz_name = f"cluster_central{central_atom}_step{central_step}_pbc.xyz"
write(xyz_name, subset)
print(f"\n Cluster snapshot saved to: {xyz_name}")

