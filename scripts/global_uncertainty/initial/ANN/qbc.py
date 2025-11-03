""" 
In addition to the functionality provided by error_uncertainty.py, 
this script picks structures in the test set with high prediction uncertainty, based on a set threshold 
This script works when the test dataset is the same in all predict.out files
"""

import numpy
import glob
from ase.io import read
import matplotlib.pyplot as plt
import shutil

def get_predictions():
    natoms=[]
    struct=[]
    predicted=[]
    predict_out = glob.glob('*predict.out')
    for outs in predict_out:
        ffile = open(outs,'r')
        lines = ffile.readlines()
        ffile.close()
        for a,line in enumerate(lines):
            splitt = line.split()
            if len(splitt) > 1:
                if splitt[0] == 'Number' and splitt[1] == 'of' and splitt[2] == 'atoms':
                    nat = splitt[4]
                    natoms.append(nat)
                if splitt[0] == 'File' and splitt[1] == 'name' and splitt[3].find('.ann') == -1: # Not the potentials, just the structures
                    current_struct = splitt[3]
                    struct.append(current_struct)
                
                    abort = False
                    i = 0
                    while abort == False:
                        i += 1
                        splitt2 = lines[a+i].split()
                        if len(splitt2) > 2:
                            if splitt2[0] == 'Total':
                                pred_val = float(splitt2[3])
                                abort = True   
                    predicted.append(pred_val)                    
    return numpy.array(predicted), numpy.array(natoms)
    
def get_reference_values():
    natoms=[]
    struct=[]
    reference=[]
    ffile = open('01_predict.out','r')
    lines = ffile.readlines()
    ffile.close()
    for a,line in enumerate(lines):
        splitt = line.split()
        if len(splitt) > 1:
            if splitt[0] == 'Number' and splitt[1] == 'of' and splitt[2] == 'atoms':
                nat = splitt[4]
                natoms.append(nat)
            if splitt[0] == 'File' and splitt[1] == 'name' and splitt[3].find('.ann') == -1: # Not the potentials, just the structures
                current_struct = splitt[3]
                struct.append(current_struct)
                # Read reference value from this file! Path is the same as in predict.out
                ref = open(current_struct,'r')
                lin = ref.readlines()
                ref.close()
                ref_val = float(lin[0].split()[-2])  
                reference.append(ref_val)                    
    return numpy.array(reference), numpy.array(natoms), numpy.array(struct)

def main():
    predict_out = glob.glob('*predict.out')
    predicted, nat = get_predictions()
    #get the total # of test images
    nat = nat.astype('float')
    count=int(len(nat)/len(predict_out))
    #evaluate energy/atom
    predicted=predicted/nat
    #split the array based on number of predict.out files
    predicted=predicted.reshape(-1,count)
    average=numpy.average(predicted,axis=0)
    sd=numpy.std(predicted,axis=0)
    reference, nat2, sys = get_reference_values()
    nat2 = nat2.astype('float')
    reference=reference/nat2
    err=reference-average
    
    syssel=sys[sd>0.01]
    print(len(syssel))
    destination = ('../sd_gt_0.01')
    for ss in syssel:
        shutil.copy(ss, destination)
    
if __name__ == '__main__':
    main()
