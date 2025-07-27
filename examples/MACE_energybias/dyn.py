import h5py
import numpy as np
import random
import time
import sys
from ase import units
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution, Stationary, ZeroRotation
from ase.md.langevin import Langevin
from ase.io import read, write
from mace.calculators import MACECalculator

np.random.seed(20)
model_paths=['01_swa.model','02_swa.model','03_swa.model','04_swa.model','05_swa.model']
mace_calc=MACECalculator(model_paths=model_paths, device='cuda', default_dtype="float64", mode="energy_bias", bias_amplitude=0.05, bias_width=0.001)
at = read('inp.xyz').copy()
at.set_calculator(mace_calc)

MaxwellBoltzmannDistribution(at, temperature_K=1000)
Stationary(at)
ZeroRotation(at)

dyn = Langevin(at, 0.5*units.fs, temperature_K=1000, friction=0.01 / units.fs, trajectory='run.traj')
print(" {:5s} {:10s} {:10s} {:5s} {:10s}".format(
    "time(fs)", "epot", "ekin", "temp(K)", "sd"))

epot=[]
sd=[]
node_sd=[]

def print_energy(atoms=at):
        time_fs=dyn.get_time()/units.fs
        pot=atoms.get_potential_energy()/len(atoms)
        ekin=atoms.get_kinetic_energy()/len(atoms)
        var=dyn.atoms.calc.results["energy_var"]
        std=np.sqrt(var)/len(atoms)
        print("{:5.1f} {:10.5e} {:10.5e} {:5.0f} {:5.7f}".format(
        time_fs, pot, ekin, ekin/(1.5*units.kB), std))
        epot.append(atoms.get_potential_energy()/len(atoms))
        sd.append(np.sqrt(var)/len(atoms))        
        node_sd.append(np.sqrt(dyn.atoms.calc.results["node_energy_var"]))
        sys.stdout.flush()

dyn.attach(print_energy, interval=1)
t0 = time.time()
dyn.run(1000)
t1 = time.time()
print("MD completed in {0:.2f} minutes!".format((t1-t0)/60))

# Convert lists to numpy arrays for storage
epot=np.array(epot)
sd=np.array(sd)
node_sd=np.array(node_sd)

# Save collected data to an HDF5 file
with h5py.File('md_data.h5', 'w') as f:
    f.create_dataset('epot', data=epot)
    f.create_dataset('sd', data=sd)
    f.create_dataset('node_sd',data=node_sd)
