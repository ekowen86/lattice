//
//  su3_relax.cpp
//
//  Created by Evan Owen on 11/1/18.
//  Copyright © 2018 Evan Owen. All rights reserved.
//
//  usage: su3_relax output N T D beta precision n_therm n_sweeps n_data parallel
//         output: output file name
//         N: lattice size (space dimension)
//         T: lattice size (time dimension)
//         D: number of dimensions
//         beta: beta
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

//    const char* output = argv[1];
    int N = atoi(argv[2]);
    int T = atoi(argv[3]);
    int D = atoi(argv[4]);
    double beta = stod(argv[5]);
    double precision = stod(argv[6]);
    precision = max(precision, numeric_limits<double>::epsilon());
    int n_therm = atoi(argv[7]);
    int n_sweeps = atoi(argv[8]);
    int n_data = atoi(argv[9]);
    bool parallel = strcmp(argv[10], "true") == 0;
    
    su3_lattice lattice = su3_lattice(N, T, D, beta);
    lattice.parallel = parallel;
//    cout << "n_threads: " << lattice.n_threads << "\n";
    lattice.thermalize(n_therm);
    
    double U[QUARK_MAX];
    double U_sum[QUARK_MAX];
    for (int q = 0; q < QUARK_MAX; q++) U_sum[q] = 0.0;

    for (int n = 0; n < n_data; n++) {
        
        // copy the lattice and relax it
        su3_lattice relax_lattice = su3_lattice(&lattice);
        long double avg_link = relax_lattice.relax(precision, true);
        
//        cout << setprecision(12);
//        for (int s = 0; s < relax_lattice.n_sites; s++) {
//            cout << relax_lattice.site[s].contract_rand(3, 4, 4, 48, 1) << "\n";
//        }

        relax_lattice.quark(U);
        for (int q = 0; q < QUARK_MAX; q++) U_sum[q] += U[q];

        cout << setprecision(12);
        cout << n;
        for (int q = 0; q < QUARK_MAX; q++) {
            cout << "\t" << U[q];
        }
        cout << "\t" << avg_link;
        cout << "\t" << relax_lattice.plaq();
        cout << "\n";

        lattice.sweep(n_sweeps);
    }

    cout << "average:";
    for (int q = 0; q < QUARK_MAX; q++) {
        cout << "\t" << U_sum[q] / n_data;
    }
    cout << "\n";

    return 0;
}
