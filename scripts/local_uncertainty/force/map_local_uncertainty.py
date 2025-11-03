import h5py
import numpy as np
from ase.io import Trajectory
from mace.calculators import MACECalculator_node_force
import logging
import os
import matplotlib.pyplot as plt

# ------------------------ CONFIG ------------------------
h5_file = "md_data.h5"             
traj_file = "run2.traj"            
force_h5_file = "force_sd.h5"      
model_paths = ["01_swa.model","02_swa.model","03_swa.model","04_swa.model","05_swa.model"]
output_txt = "first_crossings.txt"
plot_file = "cluster_species.png"
device = "cuda"
dtype = "float64"
time_per_step_ps = 0.0005  # 0.5 fs per step
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
with open(output_txt, "w", encoding="utf-8") as f:
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

# --- Step 7. Generate plot ---
to_ps = lambda step: step * time_per_step_ps
x_node, y_node = [to_ps(node_first[a]) for a in only_node], list(only_node)
x_force, y_force = [to_ps(force_first[a]) for a in only_force], list(only_force)
x_both, y_both = [to_ps(node_first[a]) for a in both], list(both)

plt.figure(figsize=(8,5))
plt.scatter(x_node, y_node, c="#D55E00", label=r'$\sigma E_{\mathrm{node}}$ only', alpha=0.7)
plt.scatter(x_force, y_force, c="#0072B2", label=r'$\sigma F$ only', alpha=0.7)
plt.scatter(x_both, y_both, c="#009E73", label="Both", alpha=0.9, marker="*", s=160)

ymax = max(y_node + y_force + y_both) if (y_node or y_force or y_both) else 0
plt.axhspan(0, 863, color="#999999", alpha=0.2)
plt.axhspan(864, 1151, color="#E69F00", alpha=0.1)
plt.axhspan(1151, ymax, color="#56B4E9", alpha=0.12)

plt.xticks(fontsize=18)
plt.yticks([], [])
plt.xlabel("Time (ps)", fontsize=20)
plt.ylabel("Bonding regime", fontsize=20)
plt.legend(fontsize=14, frameon=False)
plt.grid(axis="y", alpha=0.3)
plt.tight_layout()
plt.savefig(plot_file, dpi=300)

logging.info(f"Plot saved to {plot_file}")
logging.info(f"Node-only: {len(only_node)}, Force-only: {len(only_force)}, Both: {len(both)}")

