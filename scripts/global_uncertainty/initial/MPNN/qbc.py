"""
this script picks structures in the test set with high prediction uncertainty, based on a set threshold
the file 'std' is generated from execution of error_uncertainty.sh 
"""

from ase.io.trajectory import Trajectory
from ase.io import read, write

filename = 'std'
sd_idx = 0 # 0-based index
traj = Trajectory('test.traj')

image_num = []
with open(filename, 'r') as file:
    for line_number, line in enumerate(file, start=1):
        columns = line.split()
        sd = float(columns[sd_idx])
        if sd > 0.002:
           print('sd=',sd, 'struct=',line_number)
           image_num.append(line_number)

for j,atoms in enumerate(traj):
    for im_n in image_num:
        if j==im_n:
           write('%d.xyz' %j,atoms)
