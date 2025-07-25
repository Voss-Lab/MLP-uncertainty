"""
The script takes an fcc Pt bulk structure with 4 atoms as input and generates multiple 8-atom bulk structures. 
Structural diversity is introduced through straining and rattling, with the resulting structures saved in an ASE database file.
"""

import numpy as np
import random
from ase import build
from ase.io import read, write
from ase.db import connect

db=connect('inp.db')
strains = (0.93,0.96,0.97,0.98,0.99,1.0,1.01,1.02,1.03,1.04,1.07)
fewerstrains = (0.95,0.975,1.0,1.025,1.05)

#Applies a relative strain in the xy-plane to the atomic simulation cell, 
#adjusting both the cell dimensions and atom positions.
def strain_atoms_xy(atoms, relstrain):
    cell = np.array(atoms.cell.tolist())
    cell[0] *= relstrain
    cell[1] *= relstrain
    atoms.set_cell(cell, scale_atoms=True)

#Applies a relative strain to the entire atomic simulation cell,
#adjusting all three cell dimensions and atom positions.    
def strain_atoms(atoms, relstrain):
    cell = np.array(atoms.cell.tolist())
    cell[0] *= relstrain
    cell[1] *= relstrain
    cell[2] *= relstrain
    atoms.set_cell(cell, scale_atoms=True)

def generate_bulk():
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
        #go through many strain values 
        #scaling lattice vectors a and b only
        bulk = build.bulk(at[0].symbol, 'fcc', a=fcclatt)
        supercell = bulk.repeat(2)
        strain_atoms_xy(supercell, s)
        inp_structs.append(supercell)

    for s in strains:
        #go through many strain values 
        #scaling all the lattice vectors
        bulk = build.bulk(at[0].symbol, 'fcc', a=fcclatt)
        supercell = bulk.repeat(2)
        strain_atoms(supercell, s)
        inp_structs.append(supercell)
        
    for s in fewerstrains:
        for n in random.sample(list(range(50)),k=7):
            #go through many strain values 
            #scaling lattice vectors a and b only
            #combined with rattle at many seed values
            bulk = build.bulk(at[0].symbol, 'fcc', a=fcclatt)
            supercell = bulk.repeat(2)
            strain_atoms_xy(supercell, s)
            supercell.rattle(stdev=0.1, seed=n)
            inp_structs.append(supercell)    
            
    for s in fewerstrains:
        for n in random.sample(list(range(50)),k=7):
            #go through many strain values 
            #scaling all lattice vectors 
            #combined with rattle at many seed values
            bulk = build.bulk(at[0].symbol, 'fcc', a=fcclatt)
            supercell = bulk.repeat(2)
            strain_atoms(supercell, s)
            supercell.rattle(stdev=0.1, seed=n)
            inp_structs.append(supercell)  
            
    for inp in inp_structs:
        db.write(inp)        
      
def main():
    generate_bulk()


if __name__ == '__main__':
    main()
