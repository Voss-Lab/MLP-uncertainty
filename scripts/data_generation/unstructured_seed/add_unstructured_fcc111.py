"""
By combining strains, adsorbate placements, followed by random displacements of top layer with adsorbate, 
this script generates a dataset of Pt slab structures with H atom and saves them to an ASE database file.
"""

import numpy as np
from ase import build
from ase.build import add_adsorbate
from ase.io import read, write
from ase.ga.utilities import closest_distances_generator
from ase.ga.utilities import get_all_atom_types
from ase.ga.utilities import atoms_too_close
from ase.db import connect

db=connect('inp.db')

fewerstrains = (0.95,0.975,1.0,1.025,1.05)
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

#Move the top layer and adsorbate along the z-axis by a specified displacement.
def move_top_layer(atoms):
   #only one adsorbate atom:
   #but last index is metal surface atom
   z = atoms.positions[-2][2]
   delta = np.zeros_like(atoms.positions)
   if True:
       #+/- 0.2 AA in each direction
       #np.random.random(3) generates values between 0 and 1
       #'-5' between +/-0.5; '*2' between +/-1; '*0.2' between +/- 0.2       
       surfvec = (np.random.random(3) - 0.5)*2  * 0.2
       surfvec[2] = max(surfvec[2], -0.08)
       #+/- 1.5
       adsvec = (np.random.random(3) - 0.5)*2 * 1.5
       #at most move 0.2AA downwards in z-direction
       adsvec[2] = max(adsvec[2], -0.2)
       for i,p in enumerate(atoms.positions):
           if abs(p[2]-z)<1e-2:
              delta[i] += surfvec
           elif (p[2]-z)>0.2:
              delta[i] += adsvec
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
    i=1         
    for s in fewerstrains:
        for j in range(10):
            for sites in fcc111sites:
                if sites == 'ontop':
                    for oh in ontoph:
                        surface = build.fcc111(at[0].symbol, a=fcclatt, size=(1,1,3), vacuum=10.0)
                        add_adsorbate(surface, a, oh, sites)
                        strain_atoms(surface, s)
                        move_top_layer(surface)
                        unique_atom_types = get_all_atom_types(surface, surface.get_atomic_numbers())
                        blmin = closest_distances_generator(atom_numbers=unique_atom_types,ratio_of_covalent_radii=0.7)
                        if (atoms_too_close(surface, blmin))== False:
                             inp_structs.append(surface)
                             i+=1               
                             
                elif sites =='bridge':
                    for bh in bridgeh:
                        surface = build.fcc111(at[0].symbol, a=fcclatt, size=(1,1,3), vacuum=10.0)
                        add_adsorbate(surface, a, bh, sites)
                        strain_atoms(surface, s)
                        move_top_layer(surface)
                        unique_atom_types = get_all_atom_types(surface, surface.get_atomic_numbers())
                        blmin = closest_distances_generator(atom_numbers=unique_atom_types,ratio_of_covalent_radii=0.7)
                        if (atoms_too_close(surface, blmin))== False:
                             inp_structs.append(surface)
                             i+=1        
                else:
                    for hh in hollowh:
                        surface = build.fcc111(at[0].symbol, a=fcclatt, size=(1,1,3), vacuum=10.0)
                        add_adsorbate(surface, a, hh, sites)
                        strain_atoms(surface, s)
                        move_top_layer(surface)
                        unique_atom_types = get_all_atom_types(surface, surface.get_atomic_numbers())
                        blmin = closest_distances_generator(atom_numbers=unique_atom_types,ratio_of_covalent_radii=0.7)
                        if (atoms_too_close(surface, blmin))== False:
                             inp_structs.append(surface)
                             i+=1
                                         
        for inp in inp_structs:
            db.write(inp)

def main():
    generate_slabs()


if __name__ == '__main__':
    main()

