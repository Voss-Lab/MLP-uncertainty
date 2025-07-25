"""
By systematically combining strains and top-layer displacements, 
this script generates a dataset of clean Pt slab structures and saves them to an ASE database file.
"""

import numpy as np
from ase import build
from ase.io import read, write
from ase.db import connect

db=connect('inp.db')
strains = (0.93,0.96,0.97,0.98,0.99,1.0,1.01,1.02,1.03,1.04,1.07)
fewerstrains = (0.95,0.975,1.0,1.025,1.05)
toplayermoves = (-0.1,-0.05,0.05,0.1,0.3,1.0)

#Applies a relative strain in the xy-plane to the atomic simulation cell, 
#adjusting both the cell dimensions and atom positions.
def strain_atoms(atoms, relstrain):
    cell = np.array(atoms.cell.tolist())
    cell[0] *= relstrain
    cell[1] *= relstrain
    atoms.set_cell(cell, scale_atoms=True)

#Move the top layer of atoms along the z-axis by a specified displacement.
def move_top_layer(atoms, disp):
    #last atom typically has highest z when constructed with ase-build
    #move all atom with same z (i.e. full top layer)
    z = atoms.positions[-1][2]
    delta = np.zeros_like(atoms.positions)
    for i,p in enumerate(atoms.positions):
        if abs(p[2]-z)<1e-2:
            delta[i][2] = disp
    atoms.positions += delta

def generate_slabs():
    at = read(Pt_Pt.xsf)
    el = at[0].symbol
    test = set([x.symbol for x in at])
    if len(test)>1:
        raise ValueError('This script is only intended to be run for monometallic systems.')
    #make sure there is more than one atom in bulk fcc unit cell
    #to get distances
    pos = at.repeat((2,2,2)).positions
    #find minimum spacing
    mindist = 1e80
    for i in range(len(pos)):
        for j in range(i+1,len(pos)):
            delta = pos[i]-pos[j]
            l = np.dot(delta,delta)
            mindist = min(l, mindist)
    mindist = mindist**0.5
    #fcc lattice constant
    fcclatt = 2.**0.5 * mindist

    inp_structs=[]
    for s in strains:
        #go through many strain values not moving top layer
        #set size=(1,1,3) to generate smaller cell
        top = 0.0
        surface = build.fcc111(at[0].symbol, a=fcclatt, size=(2,1,3), vacuum=10.0)
        strain_atoms(surface, s)
        inp_structs.append(surface)

    for s in fewerstrains:
        #go through fewer strain values combined with top layer
        #displacements
        for top in toplayermoves:
            surface = build.fcc111(at[0].symbol, a=fcclatt, size=(2,1,3), vacuum=10.0)
            strain_atoms(surface, s)
            move_top_layer(surface, top)
            inp_structs.append(surface)
            
    for inp in inp_structs:
        db.write(inp)

def main():
    generate_slabs()


if __name__ == '__main__':
    main()
