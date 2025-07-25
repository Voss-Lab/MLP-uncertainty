"""
This script processes Quantum Espresso output files to extract atomic positions, cell parameters, total energy, and forces 
Parsed data is output in the XSF format 
Usage: python <script_name.py> <pwfilename.out> <filename.xsf>
"""

from sys import argv,exit,stderr
from ase import Atoms,Atom

if len(argv)!=3 or len(argv)==3 and argv[1][-4:]=='.inp':
    print('usage: '+argv[0]+' pwfilename.out filename.xsf\nalso expects to find pwfilename.inp corresponding to pwfilename.out', file=stderr)
    exit(1)

rydberg = 0.5 * 27.211386245988
bohr = 0.529177210903

def atomic_reference(species):
    atoms = {         'H': -16.77118161, # adsorbates referenced to Pt ontop
                      'Pt': -2889.448072}  # Pt referenced to Pt bulk atom
    value = atoms[species]
    return value

out = open(argv[1], 'r')
s = out.readline()
while s!='':
    if s.find('number of atoms/cell')>-1:
        natoms = int(s.split()[-1])
        break
    s = out.readline()

at = Atoms()
cell = []
syms = []
coords = []
inp = open(argv[1].replace('.out','.inp'), 'r')
s = inp.readline()
while s!='':
    if s.find('CELL_PARAMETERS')>-1:
        for i in range(3):
            cell.append([float(x) for x in inp.readline().replace('d','e').split()])
    if s.find('ATOMIC_POSITIONS')>-1:
        if s.find('crystal')>-1:
            for i in range(natoms):
                a,b,c,d = inp.readline().split()
                while a[-1].isdigit():
                    a = a[:-1]
                syms.append(a)
                coords.append((float(b.replace('d','e')),float(c.replace('d','e')),float(d.replace('d','e'))))
        else:
            print("Don't know how to parse non-crystal unit coordinates", file=stderr)
            exit(2)
    s = inp.readline()
inp.close()

for a,b in zip(syms,coords):
    at.append(Atom(a, b))
at.set_cell(cell, scale_atoms=True)

forces = []
#find energy and forces in output
output = out.readlines()
out.close()
j = 0
for i,s in enumerate(output):
    if s.find('!    total energy')>-1:
        j = i
energy = float(output[j].split()[-2])
while output[j].find('smearing contrib. (-TS)')<0:
    j += 1
energy -= 0.5*float(output[j].split()[-2])
#convert energy Ry->eV
energy *= rydberg

for s in syms:
    energy -= atomic_reference(s)

while output[j].find('force =')<0:
    j += 1
for i in range(j,j+natoms):
    forces.append([float(x)*rydberg/bohr for x in output[i].split()[-3:]])

xsf = open(argv[2], 'w')
print('# total energy = '+str(energy)+ ' eV', file=xsf)
print('\nSLAB\nPRIMVEC', file=xsf)
for c in cell:
    print('%.14f %.14f %.14f' % tuple(c), file=xsf)
print('PRIMCOORD', file=xsf)
print('%d 1' % natoms, file=xsf)
for i,a in enumerate(at):
    print('%s  %15.10f %15.10f %15.10f %15.10f %15.10f %15.10f' % (a.symbol, a.position[0], a.position[1], a.position[2], forces[i][0], forces[i][1], forces[i][2]), file=xsf)
xsf.close()
