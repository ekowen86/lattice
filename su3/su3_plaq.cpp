//
//  su3_plaq.cpp
//
//  Created by Evan Owen on 10/19/18.
//  Copyright © 2018 Evan Owen. All rights reserved.
//
//  usage: su3_plaq N T D beta_min beta_max beta_step precision n_therm n_sweeps n_data parallel
//         N: lattice size (space dimension)
//         T: lattice size (time dimension)
//         D: number of dimensions
//         beta_min: minimum beta value
//         beta_max: maximum beta value
//         beta_step: beta increment
//         precision: target precision for relaxation (0 for double-precision)
//         n_therm: minimum number of thermalization sweeps
//         n_sweeps: number of sweeps per data collection
//         n_data: number of data values to collect (zero to collect forever)
//         parellel: use parallel processors for gauge fixing

#include <iostream>
#include <iomanip>
#include "su3.hpp"

using namespace std;

int main(int argc, const char * argv[]) {
    
    int N = atoi(argv[1]);
    int T = atoi(argv[2]);
    int D = atoi(argv[3]);
    double beta_min = stod(argv[4]);
    double beta_max = stod(argv[5]);
    double beta_step = stod(argv[6]);
    double precision = stod(argv[7]);
    precision = max(precision, numeric_limits<double>::epsilon());
    int n_therm = atoi(argv[8]);
    int n_sweeps = atoi(argv[9]);
    int n_data = atoi(argv[10]);
    bool parallel = strcmp(argv[11], "true") == 0;

//    double U[60];

    for (double beta = beta_min; beta <= beta_max; beta += beta_step) {
        
        su3_lattice lattice = su3_lattice(N, T, D, beta);
        lattice.parallel = parallel;
        lattice.thermalize(n_therm);
        
        cout << setprecision(3);
        cout << "beta = " << beta << "\n\n";
        cout << setprecision(12);

        cout << lattice.plaq() * 3 << "\n";

        for (int n = 0; n < n_data; n++) {
            
            lattice.sweep(n_sweeps);
            cout << lattice.plaq() * 3 << "\n";
//            // copy the lattice and relax it
//            su3_lattice relax_lattice = su3_lattice(&lattice);
//            relax_lattice.relax(precision, true);
//            relax_lattice.quark_rand(U);
//            for (int q = 0; q < 60; q++){
//                cout << U[q] << "\t";
//            }
//            cout << "\n";
//            lattice.sweep(n_sweeps);
        }
        
        cout << "\n\n";
    }
    return 0;
}
