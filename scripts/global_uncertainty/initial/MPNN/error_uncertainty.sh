"""
This script predicts the energies of test structures using each MPNN ensemble member, 
saves the outputs as individual files (e.g., 01_out, 02_out, ...), 
and processes the results to visualize the error-uncertainty correlation 
"""

mkdir energies
tr '\n' '%' < test.xyz |sed 's/%Lattice/ /g'|tr '%' '\n'|grep pbc|sed s/=.*energy=//g|awk '{printf("%.8f\n", $2/$1)}' > energies/ref_per_atom
tr '\n' '%' < 01_out.xyz |sed 's/%Lattice/ /g'|tr '%' '\n'|grep pbc|sed s/=.*energy=//g|awk '{printf("%.8f\n", $2/$1)}' > energies/01
tr '\n' '%' < 02_out.xyz |sed 's/%Lattice/ /g'|tr '%' '\n'|grep pbc|sed s/=.*energy=//g|awk '{printf("%.8f\n", $2/$1)}' > energies/02
tr '\n' '%' < 03_out.xyz |sed 's/%Lattice/ /g'|tr '%' '\n'|grep pbc|sed s/=.*energy=//g|awk '{printf("%.8f\n", $2/$1)}' > energies/03
tr '\n' '%' < 04_out.xyz |sed 's/%Lattice/ /g'|tr '%' '\n'|grep pbc|sed s/=.*energy=//g|awk '{printf("%.8f\n", $2/$1)}' > energies/04
tr '\n' '%' < 05_out.xyz |sed 's/%Lattice/ /g'|tr '%' '\n'|grep pbc|sed s/=.*energy=//g|awk '{printf("%.8f\n", $2/$1)}' > energies/05
cd energies

paste 0* > data
awk '{
    sum = 0; sumsq = 0; n = 0;
    for (i = 1; i <= NF; i++) {
        sum += $i;
        sumsq += $i*$i;
        n++;
    }
    mean = sum / n;
    stdev = sqrt((sumsq - sum*sum/n) / (n-1));
    printf "%.4f %.4f\n", mean, stdev;
}' data > mean_std
awk '{print $1}' mean_std > mean
awk '{print $2}' mean_std > std


awk '{
    getline line < "mean";
    split(line, second_values);
    for (i = 1; i <= NF; i++) {
        diff = $i - second_values[i];
        printf "%.4f ", (diff >= 0 ? diff : -diff);
    }
    printf "\n";
}' ref_per_atom > err


"""
#use the below script to generate the Error-uncertainty correlation plot 

import matplotlib.pyplot as plt
import numpy as np

sd = np.loadtxt("std")
err = np.loadtxt("err")
plt.plot(sd,err,'x',color='blue',markersize=3)
plt.plot(err,err,'-.',color='black')
plt.xlabel('Ensemble SD, eV/atom',fontsize=12)
plt.ylabel('Error relative to DFT, eV/atom',fontsize=12)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.xlim(None,0.02)
plt.ylim(None,0.02)
plt.tight_layout()
plt.savefig('err_sd.png')
"""