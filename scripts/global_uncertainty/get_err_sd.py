
[sumanbr@sdfiana006 MPNN]$ cat error_uncertainty.sh
"""
This script predicts the energies of test structures using each MPNN ensemble member,
saves the outputs as individual files (e.g., 01_out, 02_out, ...),
and stores the computed error and uncertainty for further analysis

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

