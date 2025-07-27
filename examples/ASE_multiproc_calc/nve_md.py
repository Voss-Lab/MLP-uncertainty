from ase.io import read, write
from ase import Atoms
from ase.calculators.mixing import MultiProcessCalculator, AverageCalculator
from aenet.ase_calculator import ANNCalculator
import numpy as np
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md.verlet import VelocityVerlet
from ase import units
import sys

def printenergy(istep, a):
    """ print the potential, kinetic and total energy """
    epot = a.get_potential_energy() / len(a)
    ekin = a.get_kinetic_energy() / len(a)
    print("{:5d} {:15.8e} {:15.8e} {:7.1f} {:15.8e} {:15.8e} {:15.8e}".format(
        istep, epot, ekin, ekin/(1.5*units.kB), epot+ekin, a.calc.get_std_energy()/len(a), a.calc.get_std_forces()/len(a)))
    sys.stdout.flush()

at =  read('POSCAR')
calcs = []
for i in range(5):
    calcs.append(MultiProcessCalculator(ANNCalculator, {'Pt':'%d_Pt.ann' % (i+1), 'H':'%d_H.ann' % (i+1)}, atoms=at))
avecalc = AverageCalculator(calcs, at)
#avecalc = AverageCalculator(calcs, at, mode='energy_bias', bias_amplitude=2.5, bias_width=0.05)

temperature = 298
md_steps = 1000
print_steps = 1
dt = 0.5
MaxwellBoltzmannDistribution(at, temperature_K=temperature)

md = VelocityVerlet(at, dt*units.fs, trajectory='run.traj')
print(" #{:5s} {:15s} {:15s} {:7s} {:15s} {:15s} {:15s}".format(
    "step", "E_pot", "E_kin", "T", "E_tot", "std_E", "std_F"))
printenergy(0, at)
istep = 0
for i in range(int(md_steps/print_steps)):
    md.run(steps=print_steps)
    istep += print_steps
    printenergy(istep, at)

for c in calcs:
    c.stop()
