//
//  su2_x_test.cpp
//
//  Created by Evan Owen on 4/3/21.
//  Copyright © 2021 Evan Owen. All rights reserved.
//
//  usage: su2_test N T N5 D beta_start beta_stop beta_inc eps5 n_sweeps n_data
//         N: lattice size (spacial dimensions)
//         T: lattice size (time dimension)
//         N5: number of coupled parallel lattices in extra dimension
//         D: number of dimensions (excluding extra dimension)
//         beta: gauge coupling
//         eps5: ratio of coupling in extra dimension to normal coupling
//         n_therm: minimum number of thermalization sweeps
//         n_sweeps: number of sweeps per data collection
//         n_data: number of data values to collect

#include <iostream>
#include <iomanip>
#include <thread>

#include "su2_x.hpp"

using namespace std;

int main(int argc, const char* argv[]) {

    int N = atoi(argv[1]); cout << "N: " << N << " ";
    int T = atoi(argv[2]); cout << "T: " << T << " ";
    int N5 = atoi(argv[3]); cout << "N5: " << N5 << " ";
    int D = atoi(argv[4]); cout << "D: " << D << endl;

    cout << setprecision(4) << fixed;
    double beta_start = stod(argv[5]); cout << "beta_start: " << beta_start << " ";
    double beta_stop = stod(argv[6]); cout << "beta_stop: " << beta_stop << " ";
    double beta_inc = stod(argv[7]); cout << "beta_inc: " << beta_inc << " ";
    double eps5 = stod(argv[8]); cout << "eps5: " << eps5 << endl;

    int n_sweeps = atoi(argv[9]); cout << "n_sweeps: " << n_sweeps << " ";
    int n_data = atoi(argv[10]); cout << "n_data: " << n_data << endl;

    cout << thread::hardware_concurrency() << " parallel cores available" << endl;

    // initialize the lattice
    su2_x_lattice lattice = su2_x_lattice(N, T, N5, D, beta_start, eps5);
    lattice.parallel = true;

    double plaq, plaq_sum;
    for (double beta = beta_start; beta <= beta_stop; beta += beta_inc) {
        cout << setprecision(4) << fixed;
        cout << endl << "beta: " << beta << endl;
        lattice.verbose = 1;
        lattice.beta = beta;
        lattice.thermalize();
        lattice.verbose = 0;
        plaq_sum = 0.0;
        cout << setprecision(6) << fixed;
        for (int n = 0; n < n_data; n++) {
            printf("%d", n);
            for (int n5 = 0; n5 < N5; n5++) {
                printf(" %e", lattice.plaq(n5));
            }
            printf("\n");
            plaq = lattice.plaq(lattice.n5_center);
            plaq_sum += plaq;
            lattice.heat_bath(n_sweeps);
            // lattice.hmc(n_sweeps);
        }
        cout << setprecision(6) << fixed;
        double plaq_avg = plaq_sum / double(n_data);
        cout << "plaq: " << plaq_avg;
        cout << " E: " << (1.0 - plaq_avg) << endl;
    }

    return 0;
}
