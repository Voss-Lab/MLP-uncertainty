'''
This script reads the 'predict.out' file produced from the execution of predict.x in aenet 
to extract predicted and reference total energies of clean and combined systems (relaxed, unrelaxed) in the test dataset, 
and compute binding energies and associated error metrics 
It also generates a parity plot of the binding energies 

    Functions:
        get_all_unrelaxed_values()    - Extracts unrelaxed predicted and reference values for combined systems and slabs.
        get_all_relaxed_values()      - Extracts relaxed predicted and reference values for combined systems and slabs.
        get_adsorption_energies()     - Computes adsorption energies for the combined and clean slabs.
        calc_errors()                 - Calculates ME, MAE, and RMSE between predicted and reference adsorption energies.
        plot_stuff()                  - Plots a parity plot comparing predicted vs reference adsorption energies and includes error metrics.
        main()                        - Main function that ties everything together, calculates adsorption energies, errors, and generates the plot.

    Output:
        A plot ('Eads_parity_all.png') and printed error metrics (ME, MAE, RMSE) are produced.

'''

import numpy
import matplotlib.pyplot as plt

def get_all_unrelaxed_values():
    ffile = open('predict.out','r')
    lines = ffile.readlines()
    ffile.close()
    #
    Ecomb_predicted_unrelax = []
    Ecomb_reference_unrelax = []
    sys_comb        = []
    Eslab_predicted_unrelax = []
    Eslab_reference_unrelax = []
    sys_slab        = []
    # Total energy
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

                abort = False
                i = 0
                while abort == False:
                    i += 1
                    splitt2 = lines[a+i].split()
                    if len(splitt2) > 2:
                        if splitt2[0] == 'Total':
                            # Predicted value!
                            pred_val = float(splitt2[3])
                            abort = True

                # store structure -> to be able to compute adsorption energies later on
                if current_struct.find('unrelaxed') != -1: 
                    if  current_struct.find('combined') != -1:
                        sys_comb.append(current_struct.split('/')[-1])
                        Ecomb_predicted_unrelax.append(pred_val)
                        Ecomb_reference_unrelax.append(ref_val)
                    elif current_struct.find('clean') != -1:
                        sys_slab.append(current_struct.split('/')[-1])
                        Eslab_predicted_unrelax.append(pred_val)
                        Eslab_reference_unrelax.append(ref_val)
    return numpy.array(Ecomb_predicted_unrelax), numpy.array(Ecomb_reference_unrelax), numpy.array(sys_comb), numpy.array(Eslab_predicted_unrelax), numpy.array(Eslab_reference_unrelax), numpy.array(sys_slab)

def get_all_relaxed_values():
    ffile = open('predict.out','r')
    lines = ffile.readlines()
    ffile.close()
    #
    Ecomb_predicted_relax = []
    Ecomb_reference_relax = []
    sys_comb2        = []
    Eslab_predicted_relax = []
    Eslab_reference_relax = []
    sys_slab2        = []
    # Total energy
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

                abort = False
                i = 0
                while abort == False:
                    i += 1
                    splitt2 = lines[a+i].split()
                    if len(splitt2) > 2:
                        if splitt2[0] == 'Total':
                            # Predicted value!
                            pred_val = float(splitt2[3])
                            abort = True

                # store structure -> to be able to compute adsorption energies later on
                if current_struct.find('unrelaxed') == -1: 
                    if current_struct.find('relaxed') != -1:  
                        if current_struct.find('combined') != -1:
                            sys_comb2.append(current_struct.split('/')[-1])
                            Ecomb_predicted_relax.append(pred_val)
                            Ecomb_reference_relax.append(ref_val)
                        elif current_struct.find('clean') != -1:
                            sys_slab2.append(current_struct.split('/')[-1])
                            Eslab_predicted_relax.append(pred_val)
                            Eslab_reference_relax.append(ref_val)
    return numpy.array(Ecomb_predicted_relax), numpy.array(Ecomb_reference_relax), numpy.array(sys_comb2), numpy.array(Eslab_predicted_relax), numpy.array(Eslab_reference_relax), numpy.array(sys_slab2)    

def get_adsorption_energies(Ecombpred,Ecombref,sys_comb,Eslabpred,Eslabref,sys_slab):
    '''
        get_adsorption_energies
        Compute adsorption energies based on the combined systems and the clean slabs
    '''
    Eads_pred = []
    Eads_ref  = []
    
    # get the first split('_') from the name of .xsf files > same for both clean and combined slabs
    for s_comb,system_comb in enumerate(sys_comb):
        # number which corresponds to comb system
        metal1_comb = system_comb.split('_')[0]
        metal2_comb = system_comb.split('_')[1]
        idx_comb = s_comb

        # go through slab systems -> find corresponding slab
        for s_slab,system_slab in enumerate(sys_slab):
            metal1_slab = system_slab.split('_')[0]
            metal2_slab = system_slab.split('_')[1]
            if metal1_slab == metal1_comb and metal2_slab == metal2_comb:
                idx_slab = s_slab
                break

        # Compute adsorption energies
        Eads_pred.append(Ecombpred[idx_comb] - Eslabpred[idx_slab])
        Eads_ref.append(Ecombref[idx_comb] - Eslabref[idx_slab])
    return numpy.array(Eads_pred), numpy.array(Eads_ref)

def calc_errors(predicted,reference):
    '''
        Calculate ME, MAE, and RMSE between predicted and reference values
    '''
    me = numpy.sum(predicted - reference)/len(reference)
    mae = numpy.sum(abs(predicted - reference))/len(reference)
    rmse = numpy.sqrt((sum((predicted - reference)**2))/len(reference))
    return me, mae, rmse
    
def plot_stuff(predicted, reference):
    # Plot all values, including global ME, MAE
    me, mae, rmse = calc_errors(predicted,reference)
    plt.title(f'{me:5.4f} eV ME; {mae:5.4f} eV MAE; {rmse:5.4f} eV RMSE',fontsize=12, color='black')
    plt.plot(reference,predicted,'o',color='orange')
    plt.plot(reference,reference,'--',color='grey')
    plt.xlabel("DFT calculated $ΔE_{ads}$ [eV]",fontsize=15, labelpad=10)
    plt.ylabel("ML predicted $ΔE_{ads}$ [eV]",fontsize=15, labelpad=10)
    plt.tight_layout()
    plt.savefig('Eads_parity_all.png')
    
def main():
    
    # calculate adsorption energies
    Ecomb_predicted, Ecomb_reference, sys_comb, Eslab_predicted, Eslab_reference, sys_slab = get_all_unrelaxed_values() 
    Ecomb_predicted2, Ecomb_reference2, sys_comb2, Eslab_predicted2, Eslab_reference2, sys_slab2 = get_all_relaxed_values()
    Eads_pred_u, Eads_ref_u = get_adsorption_energies(Ecomb_predicted, Ecomb_reference, sys_comb, Eslab_predicted, Eslab_reference, sys_slab)
    Eads_pred_r, Eads_ref_r = get_adsorption_energies(Ecomb_predicted2, Ecomb_reference2, sys_comb2, Eslab_predicted2, Eslab_reference2, sys_slab2)
    Eads_pred = numpy.concatenate((Eads_pred_u, Eads_pred_r))
    Eads_ref = numpy.concatenate((Eads_ref_u, Eads_ref_r))
    me_tot, mae_tot, rmse_tot = calc_errors(Eads_pred, Eads_ref)
    print(f'ME: {me_tot:10.4f} eV')
    print(f'MAE: {mae_tot:10.4f} eV')
    print(f'RMSE: {rmse_tot:10.4f} eV')
    plot_stuff(Eads_pred, Eads_ref)


if __name__ == '__main__':
    main()
