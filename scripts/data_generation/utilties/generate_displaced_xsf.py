"""
This script generates displaced atomic structures by applying small displacements (delta) to atomic positions 
approximate energy changes are computed using a first-order force projection
small in-plane forces below a defined threshold (force_threshold) are ignored to avoid redundant training
"""

import glob
import numpy as np
from ase import io

delta = 0.05 #small displacement in Angstroms to teach aenet the forces
force_threshold = 0.001 #force below which in-plane forces are not trained against

xsfs = glob.glob('*.xsf') #training structures in xsf format
for xsf in xsfs:
    if xsf.find('_FFF_')>-1: #skip xsf files generated from forces,
                           #if they already exist (which have '_FFF_' in their filenames
        continue

    at = io.read(xsf)
    cell = at.cell.tolist()
    pos = at.positions
    forces = at.get_forces() / 27.211386245988 #ase thinks the energies are in Hartree
    fi = open(xsf, 'r')
    energy = float(fi.readline().split()[-2])
    fi.close()
    natoms = len(at)
    for iatom in range(natoms):  #for every atom:
        for j in range(3):   #x,y,z-direction
            if j<2 and abs(forces[iatom][j]) < force_threshold:
                continue #for tiny in-plane forces move on to next loop iteration
                         #without creating displaced XSF file
            for k in (-1,1): #plus/minus displacement
                displace = np.zeros_like(pos)
                displace[iatom][j] = delta*k
                newpos = pos + displace
                newenergy = energy - np.dot(displace.flat,forces.flat) #project forces onto displacement to get approximate energy change

                newxsf = open(xsf.replace('.xsf','_FFF_%d_%d_%d.xsf' % (iatom,j,k)), 'w')
                print('# total energy = '+str(newenergy)+ ' eV', file=newxsf)
                print('\nSLAB\nPRIMVEC', file=newxsf)
                for c in cell:
                    print('%.12f %.12f %.12f' % tuple(c), file=newxsf)
                print('PRIMCOORD', file=newxsf)
                print('%d 1' % natoms, file=newxsf)
                for i,a in enumerate(at):
                    print('%s  %.12f %.12f %.12f %.12f %.12f %.12f' % (a.symbol, newpos[i][0], newpos[i][1], newpos[i][2], 0., 0., 0.), file=newxsf) #simply putting zero for the forces on the displaced atoms: we don't know those without doing extra DFT
                newxsf.close()
