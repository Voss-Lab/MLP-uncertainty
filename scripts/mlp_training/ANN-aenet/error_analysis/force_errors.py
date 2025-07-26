"""
This script reads the 'predict.out' file produced from execution of predict.x in aenet 
to extract predicted and reference energies, forces of test dataset
and computes force error metrics.

    Functions:
        get_all_values()       - Extracts predicted and reference total energies and atomic forces from the 'predict.out' file.
        get_force_errors()     - Calculates maximum force error, RMS force error, and MAE force error for each structure.
        main()                 - Main function that processes the data, calculates force errors, and prints the average force errors.

    Output:
        Printed values for average "Max force error", "RMS force error", and "MA force error" (in eV/Angstrom).
""" 

import numpy
from ase.io import read
import matplotlib.pyplot as plt

def get_all_values():
    '''
        get_all_values
        Read the predict.out file
    '''
    ffile = open('predict.out','r')
    lines = ffile.readlines()
    ffile.close()
    #
    predicted = []
    reference = []
    predicted_forces = []
    reference_forces = []
    # Total energy and forces
    for a,line in enumerate(lines):
        splitt = line.split()
        if len(splitt) > 1:
            if splitt[0] == 'File' and splitt[1] == 'name' and splitt[3].find('.ann') == -1: # Not the potentials, just the structures
                current_struct = splitt[3]
                # Read reference value from this file! Path is the same as in predict.out
                ref = open(current_struct,'r')
                lin = ref.readlines()
                ref.close()
                ref_val = float(lin[0].split()[-2])
                ref_forces = read(current_struct).get_forces()/27.211386245988
                pred_forces = numpy.zeros_like(ref_forces)

                i = 0
                while lines[a+i].find('corresponding atomic forces')<0:
                    i += 1
                while lines[a+i].find('-----')<0:
                    i += 1
                i += 1
                for j in range(len(pred_forces)):
                    pred_forces[j][:] = [float(x) for x in lines[a+i].split()[-3:]]
                    i += 1
                while lines[a+i].find('Total energy')<0:
                    i += 1
                pred_val = float(lines[a+i].split()[3])

                predicted.append(pred_val)
                reference.append(ref_val)
                predicted_forces.append(pred_forces)
                reference_forces.append(ref_forces)
    return numpy.array(predicted), numpy.array(reference), predicted_forces, reference_forces

def get_force_errors(pforces,rforces):
    #max force error on individual atom
    max_forces = numpy.empty(len(pforces))
    #RMS error
    rms_forces = numpy.empty(len(pforces))
    #MAE error
    mae_forces = numpy.empty(len(pforces))
    for i,(p,r) in enumerate(zip(pforces,rforces)):
        diff = p - r
        #length of difference force vector on each atom:
        d = numpy.linalg.norm(diff,axis=1)
        max_forces[i] = numpy.max(d)
        rms_forces[i] = numpy.average(d**2)**0.5
        mae_forces[i] = numpy.average(d)
    return max_forces, rms_forces, mae_forces

def main():
    predicted, reference, predicted_forces, reference_forces = get_all_values()
    max_forces, rms_forces, mae_forces = get_force_errors(predicted_forces, reference_forces)
    #max_...,rms_..., and mae_forces are numpy array with length equaling to number of structures in predict.out
    print(f'Average "Max force error": {numpy.average(max_forces):10.6f} eV/Ang')
    print(f'Average "RMS force error": {numpy.average(rms_forces):10.6f} eV/Ang')
    print(f'Average "MA force error": {numpy.average(mae_forces):10.6f} eV/Ang')


if __name__ == '__main__':
    main()
