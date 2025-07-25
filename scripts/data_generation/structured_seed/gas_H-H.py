"""
Script for generating H2 molecule structures with varying bond lengths.

This script creates a set of H2 molecules with bond lengths ranging from 
0.5 Å to 1.0 Å (or customizable range), stores them in a cubic simulation 
box of dimensions 10 Å x 10 Å x 10 Å, and organizes them into directories 
named based on the bond length.

"""

import os
import numpy as np
from ase import Atoms

dr='/path_to_dir'
os.chdir(dr)
inp_structs=[]

# Loop over bond lengths from 0.5 Å to 1.0 Å, with 30 evenly spaced values
for bl in np.linspace(0.5,1,num=30):
# fewer samples for bond lengths from 1 Å to 5 Å
# for bl in np.linspace(1, 5, num=20): 
    d = '%0.2f' %bl+'_'+'H2'
    os.makedirs(d, exist_ok=True)
    # vary cell to get box with sides 5,10, and 15A
    h2=Atoms('H2',cell=[10.0,10.0,10.0],positions=[(0,0,0),(0,0,bl)])
    h2.center()
    inp_structs.append(h2)