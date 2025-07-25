"""
This script reads the 'predict.out' file produced from the execution of predict.x in aenet 
to extract predicted and reference total energies of test dataset 
and compute error metrics in per-atom energies
It also generates a plot comparing the predicted energies against the reference DFT values

    Functions:
        get_all_values()       - Extracts predicted energies, reference energies, and the number of atoms for each structure from the 'predict.out' file.
        calc_errors()          - Calculates ME, MAE, and RMSE between predicted and reference energies.
        plot_stuff()           - Plots a parity plot comparing predicted vs reference energies and includes error metrics in the plot title.
        main()                 - Main function that ties everything together, calculates the per-atom energies, errors, and generates the plot.

    Output:
        A plot ('energy.png') and printed error metrics (ME, MAE, RMSE) are produced.
""" 
 
import numpy
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
    natoms = []
    # Total energy
    for a,line in enumerate(lines):
        splitt = line.split()
        if len(splitt) > 1:
            if splitt[0] == 'Number' and splitt[1] == 'of' and splitt[2] == 'atoms':
                nat = splitt[4]
                natoms.append(nat)
            if splitt[0] == 'File' and splitt[1] == 'name' and splitt[3].find('.ann') == -1: # Not the potentials, just the structures
                current_struct = splitt[3]
                # Read reference value from this file! Path is the same as in predict.out
                ref = open(current_struct,'r')
                lin = ref.readlines()
                ref.close()
                ref_val = float(lin[0].split()[-2])
            
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
                reference.append(ref_val)                
    return numpy.array(predicted), numpy.array(reference), numpy.array(natoms)

def calc_errors(predicted,reference):
    me = numpy.sum(predicted - reference)/len(reference)
    mae = numpy.sum(abs(predicted - reference))/len(reference)
    rmse = numpy.sqrt((sum((predicted - reference)**2))/len(reference))
    return me, mae, rmse

def plot_stuff(predicted,reference):
    # Plot all values, including global ME, MAE
    me, mae, rmse = calc_errors(predicted,reference)
    plt.title(f'ME = {me:5.4f} eV ; MAE = {mae:5.4f} eV ; RMSE = {rmse:5.4f} eV',fontsize=12)
    plt.plot(reference,predicted,'o',color='orange')
    plt.plot(reference,reference,'--',color='black')
    plt.xlabel('DFT reference E [eV]',fontsize=15)
    plt.ylabel('ML predicted E [eV]',fontsize=15)
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    plt.savefig('energy.png')

def main():
    predicted, reference, nat = get_all_values()
    nat = nat.astype('float')
    predicted, reference = (predicted/nat), (reference/nat)  
    me_tot, mae_tot, rmse_tot = calc_errors(predicted,reference)
    print(f'ME  (all): {me_tot:10.4f} eV')
    print(f'MAE (all): {mae_tot:10.4f} eV')
    print(f'RMSE (all): {rmse_tot:10.4f} eV')
    plot_stuff(predicted,reference)

if __name__ == '__main__':
    main()
