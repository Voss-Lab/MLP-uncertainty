"""
This script analyzes the accuracy of predicted energies from ANN ensemble potentials  
by comparing them to reference values from DFT.
A scatter plot of DFT error v/s prediction standard deviation (uncertainty) is generated
The script works when the test dataset is the same in all predict.out files
"""

import numpy
import scipy
import glob
from ase.io import read
import matplotlib.pyplot as plt

def get_predictions():
    natoms=[]
    struct=[]
    predicted=[]
    # output files from execution of predict.x on the test structures by ensemble members 
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
    return numpy.array(reference), numpy.array(natoms)

def main():
    predict_out = glob.glob('*predict.out')
    predicted, nat = get_predictions()
    #get the total # of test images
    nat = nat.astype('float')
    #evaluate energy/atom
    count=int(len(nat)/len(predict_out))
    #split the array based on number of predict.out files
    predicted=predicted/nat
    predicted=predicted.reshape(-1,count)
    average=numpy.average(predicted,axis=0)
    sd=numpy.std(predicted,axis=0)
    reference, nat2 = get_reference_values()
    nat2 = nat2.astype('float')
    reference=reference/nat2
    err=reference-average
    
    plt.plot(sd,err,'.',color='orange')
    plt.plot(err,err,'-.',color='grey')
    plt.xlabel('Standard deviation, eV/atom',fontsize=12)
    plt.ylabel('Error relative to DFT, eV/atom',fontsize=12)
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    plt.savefig('DFT_err_vs_SD.png')  

    
if __name__ == '__main__':
    main()
