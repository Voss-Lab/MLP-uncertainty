"""
By systematically combining strains and adsorbate placements, 
this script generates a dataset of Pt slab structures with H atom and saves them to an ASE database file.
"""

import numpy as np
from ase import build
from ase.build import add_adsorbate
from ase.io import read, write
from ase.db import connect

db=connect('inp.db')
strains = (0.93,0.96,0.98,0.99,1.0,1.01,1.02,1.04,1.07)
fcc111sites = ('ontop', 'bridge', 'fcc', 'hcp')
ontoph = (1.4, 1.6, 1.8, 2.0, 2.5, 5.0)
bridgeh = (0.6, 0.8, 1.0, 1.2, 1.5, 3.0, 5.0)
hollowh = (0.4, 0.6, 0.8, 1.0, 1.5, 3.0, 5.0)
ads = 'H'

#Applies a relative strain in the xy-plane to the atomic simulation cell, 
#adjusting both the cell dimensions and atom positions.
def strain_atoms(atoms, relstrain):
    cell = np.array(atoms.cell.tolist())
    cell[0] *= relstrain
    cell[1] *= relstrain
    atoms.set_cell(cell, scale_atoms=True)

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
        for sites in fcc111sites:        
            if sites == 'ontop':
                for oh in ontoph:
                    #set size=(2,1,3) for larger u.c. i.e., 0.5 ML coverage
                    surface = build.fcc111(at[0].symbol, a=fcclatt, size=(1,1,3), vacuum=10.0)
                    add_adsorbate(surface, ads, oh, sites)
                    strain_atoms(surface, s)
                    inp_structs.append(surface)
            elif sites == 'bridge':
                for bh in bridgeh:
                    surface = build.fcc111(at[0].symbol, a=fcclatt, size=(1,1,3), vacuum=10.0)
                    add_adsorbate(surface, ads, bh, sites)
                    strain_atoms(surface, s)
                    inp_structs.append(surface)
            else:
                for hh in hollowh:
                    surface = build.fcc111(at[0].symbol, a=fcclatt, size=(1,1,3), vacuum=10.0)
                    add_adsorbate(surface, ads, hh, sites)
                    strain_atoms(surface, s)
                    inp_structs.append(surface)
                    
    for inp in inp_structs:
        db.write(inp)
                                              
def main():
    generate_slabs()


if __name__ == '__main__':
    main()

