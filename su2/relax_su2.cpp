//
//  relax_su2.cpp
//
//  Created by Evan Owen on 4/3/18.
//  Copyright © 2018 Evan Owen. All rights reserved.
//
//  usage: relax_su2 output N D beta precision gauge n_therm n_sweeps n_data
//         output: output file name
//         N: lattice size
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
#include "su2.hpp"

using namespace std;

int main(int argc, const char * argv[]) {

    const char* output = argv[1];
    int N = atoi(argv[2]);
    int D = atoi(argv[3]);
    double beta = stod(argv[4]);
    double precision = stod(argv[5]);
    precision = max(precision, numeric_limits<double>::epsilon());
    bool coulomb = strcmp(argv[6], "C") ? false : true;
    int n_therm = atoi(argv[7]);
    int n_sweeps = atoi(argv[8]);
    int n_data = atoi(argv[9]);
    bool forever = false;
    if (n_data == 0) forever = true;
    
    int n_prop = N / 2; // maximum propagator length/separation

    // create a lattice and thermalize it
    su2 lattice = su2(N, D, beta);
    n_therm = lattice.thermalize(n_therm);

    for (int n = 0; n < n_data || forever; n++) {
        
        // copy the lattice and relax it
        su2 relax_lattice = su2(&lattice);
        double dU = 1;
        double U1 = relax_lattice.relax(coulomb);
        double U2;
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
        
        relax_lattice.write_four_point(output, n_prop, coulomb);
        cout << "n = " << n << "\n\n";
        
        lattice.sweep(n_sweeps);
    }

    return 0;
}
