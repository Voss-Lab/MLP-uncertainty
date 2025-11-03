"""
This script generates multiple input files (generate.in) to execute using generate.x in aenet
It uses bootstrapping to sample structures for training 
"""

import numpy
import glob

def write_input(data,metals,species,n):
    ffile = open(str(n)+'_'+'generate.in','w')
    ffile.write("OUTPUT  ref.train\n")
    ffile.write("\n")
    ffile.write("TYPES\n")
    ffile.write(f"{len(metals)+len(species)}\n")
    for ads in species:
        ffile.write(f'{ads} 0.0 ! eV\n')
    for metal in metals:
        ffile.write(f'{metal} 0.0 ! eV\n')
    ffile.write("\n")
    ffile.write("SETUPS\n")
    for ads in species:
        ffile.write(f'{ads}   {ads}.stp\n')
    for metal in metals:
        ffile.write(f'{metal}   {metal}.stp\n')
    ffile.write("\n")
    ffile.write("FILES\n")
    ffile.write(f"{len(data)}\n")
    for entry in data:
        ffile.write(f"{entry}\n")
    ffile.close()

def main():
    metals = ['Pt']
    adsorbates = ['H']
    all_structs = []
    all_structs.extend(glob.glob('../*.xsf')) #all xsf files in the training set
    rng = numpy.random.RandomState(1)
    for n in range(5):
        bootstrap_indices = rng.choice(numpy.arange(len(all_structs)), size=len(all_structs), replace=True,)
        structs = [all_structs[i] for i in bootstrap_indices]
        write_input(structs, metals, adsorbates, n)    

if __name__ == '__main__':
    main()
