import h5py
import numpy as np
from ase.io import Trajectory
from mace.calculators import MACECalculator_node_force
import logging
import os

# ------------------------ CONFIG ------------------------
h5_file = "md_data.h5"             # input file (contains node_sd)
traj_file = "run2.traj"            # trajectory
force_h5_file = "force_sd.h5"      # output file for forcenormstd
model_paths = ["01_swa.model","02_swa.model","03_swa.model","04_swa.model","05_swa.model"]
output_txt = "first_crossings.txt"
device = "cuda"
dtype = "float64"
# --------------------------------------------------------

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

# --- Step 1. Load node_sd from HDF5 ---
with h5py.File(h5_file, "r") as f:
    node_sd = f["node_sd"][()]   # shape: (nsteps, natoms)
nsteps, natoms = node_sd.shape
logging.info(f"Loaded node_sd from {h5_file}, shape = {node_sd.shape}")

# --- Step 2. Load or compute forcenormstd ---
if os.path.exists(force_h5_file):
    try:
        with h5py.File(force_h5_file, "r") as f:
            force_sd_all = f["forcenormstd"][()]
        logging.info(f"Loaded forcenormstd from {force_h5_file}, shape = {force_sd_all.shape}")
    except KeyError:
        force_sd_all = None
else:
    force_sd_all = None

if force_sd_all is None:
    logging.info("forcenormstd not found, computing with MACE calculator...")
    traj = Trajectory(traj_file)
    calc = MACECalculator_node_force(model_paths=model_paths, device=device, default_dtype=dtype)

    force_list = []
    for step, atoms in enumerate(traj):
        atoms.calc = calc
        _ = atoms.get_forces()
        forcenormstd = atoms.calc.results.get("forcenormstd", None)
        if forcenormstd is not None:
            force_list.append(forcenormstd)
        if step % 50 == 0:
            logging.info(f"Processed step {step}/{len(traj)}")

    force_sd_all = np.vstack(force_list)
    with h5py.File(force_h5_file, "w") as f:
        f.create_dataset("forcenormstd", data=force_sd_all)
    logging.info(f"Saved forcenormstd to {force_h5_file}")

# --- Step 3. Compute 3σ thresholds ---
mean_node, std_node = np.mean(node_sd), np.std(node_sd)
node_thresh = mean_node + 3 * std_node

mean_force, std_force = np.mean(force_sd_all), np.std(force_sd_all)
force_thresh = mean_force + 3 * std_force

# --- Step 4. Identify first crossings ---
def find_first_crossings(values, threshold):
    """Return {atom_idx: first_step_idx} where value > threshold."""
    nsteps, natoms = values.shape
    first_cross = {}
    for step in range(nsteps):
        above = np.where(values[step] > threshold)[0]
        for atom_idx in above:
            if atom_idx not in first_cross:
                first_cross[atom_idx] = step
        if len(first_cross) == natoms:
            break
    return first_cross

node_first = find_first_crossings(node_sd, node_thresh)
force_first = find_first_crossings(force_sd_all, force_thresh)

# --- Step 5. Compare sets ---
node_set = set(node_first.keys())
force_set = set(force_first.keys())

only_node = node_set - force_set
only_force = force_set - node_set
both = node_set & force_set

# --- Step 6. Write summary to file ---
with open(output_txt, "w") as f:
    f.write("3σ thresholds:\n")
    f.write(f"Node energy SD: {node_thresh:.4f} (mean={mean_node:.4f}, std={std_node:.4f})\n")
    f.write(f"Force norm std: {force_thresh:.4f} (mean={mean_force:.4f}, std={std_force:.4f})\n\n")

    f.write("First crossing for node_energy_sd (atom : step):\n")
    for atom, step in sorted(node_first.items()):
        f.write(f"Atom {atom} : Step {step}\n")

    f.write("\nFirst crossing for forcenormstd (atom : step):\n")
    for atom, step in sorted(force_first.items()):
        f.write(f"Atom {atom} : Step {step}\n")

    f.write(f"\nAtoms flagged only by node_energy_sd: {len(only_node)}\n")
    f.write(f"Atoms flagged only by forcenormstd: {len(only_force)}\n")
    f.write(f"Atoms flagged by both: {len(both)}\n")

logging.info(f"Results written to {output_txt}")
logging.info(f"Node-only: {len(only_node)}, Force-only: {len(only_force)}, Both: {len(both)}")

