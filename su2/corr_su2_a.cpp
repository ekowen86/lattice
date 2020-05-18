//
//  corr_su2_a.cpp
//
//  Created by Evan Owen on 9/19/18.
//  Copyright © 2018 Evan Owen. All rights reserved.
//
//  usage: corr_su2 output1 output2 N T D beta precision gauge n_therm n_sweeps n_data
//         output1: output file name for four-point data
//         output2: output file name for correlator data
//         N: lattice size (spacial)
//         T: lattice size (time)
//         D: number of dimensions
//         beta: beta
//         precision: target precision for relaxation (0 for double-precision)
//         gauge: "C" for coulomb, "L" for landau
//         n_therm: minimum number of thermalization sweeps
//         n_sweeps: number of sweeps per data collection
//         n_data: number of data values to collect (zero to collect forever)

#include <iostream>
#include <string>
#include <limits>
#include <fstream>
#include "su2_a.hpp"

using namespace std;

int main(int argc, const char * argv[]) {
    
    const char* output1 = argv[1];
    const char* output2 = argv[2];
    int N = atoi(argv[3]);
    int T = atoi(argv[4]);
    int D = atoi(argv[5]);
    double beta = stod(argv[6]);
    long double precision = stod(argv[7]);
    precision = max(precision, numeric_limits<long double>::epsilon());
    bool coulomb = strcmp(argv[8], "C") ? false : true;
    int n_therm = atoi(argv[9]);
    int n_sweeps = atoi(argv[10]);
    int n_data = atoi(argv[11]);
    bool forever = false;
    if (n_data == 0) forever = true;
        
    // create a lattice and thermalize it
    su2_a lattice = su2_a(N, T, D, beta);
    n_therm = lattice.thermalize(n_therm);
    
    for (int n = 0; n < n_data || forever; n++) {
        
        // copy the lattice and relax it
        su2_a relax_lattice = su2_a(&lattice);
        long double dU = 1;
        long double U1 = relax_lattice.relax(coulomb);
        long double U2;
        while (dU > precision) {
            U2 = U1;
            for (int i = 0; i < 10; i++) {
                // do 10 relaxation sweeps
                U1 = relax_lattice.relax(coulomb);
            }
            dU = abs(U1 - U2) / abs(U2);
            cout << "U = " << U1 << " / ";
            cout << "dU = " << dU << "\n";
        }
        
        relax_lattice.write_four_point(output1, coulomb);
        relax_lattice.write_correlator(output2, coulomb);
        cout << "n = " << n << "\n\n";
        
        lattice.sweep(n_sweeps);
    }
    
    return 0;
}
