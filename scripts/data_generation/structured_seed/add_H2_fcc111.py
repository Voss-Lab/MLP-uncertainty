"""
By systematically combining strains and adsorbate placements, 
this script generates a dataset of Pt slab structures with H2 molecule and saves them to an ASE database file.
"""

import numpy as np
from ase import build
from ase.build import add_adsorbate
from ase.build import molecule
from ase.io import read, write
from ase.ga.utilities import closest_distances_generator
from ase.ga.utilities import get_all_atom_types
from ase.ga.utilities import atoms_too_close
from ase.ga.utilities import gather_atoms_by_tag
from ase.db import connect

db=connect('inp.db')
fewerstrains = (0.93, 0.96, 0.99, 1.0, 1.01, 1.04, 1.07)
fcc111sites = ('ontop', 'bridge', 'fcc', 'hcp')
ontoph = (1.4, 1.6, 1.8, 2.0, 2.5, 3.0)
bridgeh = (0.6, 0.8, 1.0, 1.2, 1.5, 3.0)
hollowh = (0.4, 0.6, 0.8, 1.0, 1.5, 3.0)
mol=molecule('H2')

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

    inp_structs = []
    for s in fewerstrains:
        for sites in fcc111sites:        
            if sites == 'ontop':
                for oh in ontoph:
                    # generate n samples at each height
                    for j in range(5):
                        #set size=(2,1,3) for larger u.c. i.e., 0.5 ML coverage
                        surface = build.fcc111(at[0].symbol, a=fcclatt, size=(1,1,3), vacuum=10.0)
                        rndmol = mol.copy()
                        # stretch H-H bond between 0.5 and 2.5 Ang
                        rndmol[1].position[2] += (2.5-0.5)*np.random.random() - 1.75
                        axis = np.random.random(3)
                        axis /= np.linalg.norm(axis)
                        rndmol.rotate(np.random.random()*90,axis)
                        gather_atoms_by_tag(rndmol)
                        #add the stretched, rotated adsorbate                            
                        add_adsorbate(surface, rndmol, oh, sites)
                        strain_atoms(surface, s)
                        z=surface.positions[-3][2]
                        z1=surface.positions[-2][2]
                        z2=surface.positions[-1][2]
                        unique_atom_types = get_all_atom_types(surface, surface.get_atomic_numbers())
                        blmin = closest_distances_generator(atom_numbers=unique_atom_types,ratio_of_covalent_radii=0.7)
                        # make sure H2 is atleast 0.2 Ang above the surface
                        # make sure the atoms are not too close
                        if (z1-z)>0.2:
                           if (z2-z)>0.2:  
                              if (atoms_too_close(surface, blmin))== False:
                                  inp_structs.append(surface)
                                  j+=1                                                         
            elif sites == 'bridge':
                for bh in bridgeh:
                    for j in range(5):
                        surface = build.fcc111(at[0].symbol, a=fcclatt, size=(1,1,3), vacuum=10.0)
                        rndmol = mol.copy()
                        rndmol[1].position[2] += (2.5-0.5)*np.random.random() - 1.75
                        axis = np.random.random(3)
                        axis /= np.linalg.norm(axis)
                        rndmol.rotate(np.random.random()*90,axis)
                        gather_atoms_by_tag(rndmol)                           
                        add_adsorbate(surface, rndmol, bh, sites)
                        strain_atoms(surface, s)
                        z=surface.positions[-3][2]
                        z1=surface.positions[-2][2]
                        z2=surface.positions[-1][2]
                        unique_atom_types = get_all_atom_types(surface, surface.get_atomic_numbers())
                        blmin = closest_distances_generator(atom_numbers=unique_atom_types,ratio_of_covalent_radii=0.7)
                        if (z1-z)>0.2:
                           if (z2-z)>0.2:  
                              if (atoms_too_close(surface, blmin))== False:
                                  inp_structs.append(surface)
                                  j+=1                                             
            else:
                for hh in hollowh:
                    for j in range(5):
                        surface = build.fcc111(at[0].symbol, a=fcclatt, size=(1,1,3), vacuum=10.0)
                        rndmol = mol.copy()
                        rndmol[1].position[2] += (2.5-0.5)*np.random.random() - 1.75
                        axis = np.random.random(3)
                        axis /= np.linalg.norm(axis)
                        rndmol.rotate(np.random.random()*90,axis)
                        gather_atoms_by_tag(rndmol)                          
                        add_adsorbate(surface, rndmol, hh, sites)
                        strain_atoms(surface, s)
                        z=surface.positions[-3][2]
                        z1=surface.positions[-2][2]
                        z2=surface.positions[-1][2]
                        unique_atom_types = get_all_atom_types(surface, surface.get_atomic_numbers())
                        blmin = closest_distances_generator(atom_numbers=unique_atom_types,ratio_of_covalent_radii=0.7)
                        if (z1-z)>0.2:
                           if (z2-z)>0.2:  
                              if (atoms_too_close(surface, blmin))== False:
                                  inp_structs.append(surface)
                                  j+=1      

        for inp in inp_structs:
            db.write(inp)   
              
def main():
    generate_slabs()

if __name__ == '__main__':
    main()

