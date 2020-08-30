//
//  su3_x_test.cpp
//  su3_x_test
//
//  Created by Evan Owen on 8/21/20.
//

#include <iostream>
#include <iomanip>

#include "su3_x.hpp"

using namespace std;

int main() {

    // print parameters
    int N = 6; cout << "N = " << N << ", ";
    int T = 6; cout << "T = " << T << ", ";
    int N5 = 3; cout << "N5 = " << N5 << ", ";
    int D = 4; cout << "D = " << D << ", ";
    double beta = 5.7; cout << "beta = " << setprecision(2) << fixed << beta << ", ";
    double eps5 = 0.5; cout << "eps5 = " << setprecision(2) << fixed << eps5 << ", ";
    int n_data = 20; cout << "n_data = " << n_data << endl;
    
    // initialize and thermalize the lattice
    su3_x_lattice lattice = su3_x_lattice(N, T, N5, D, beta, eps5);
    lattice.parallel = true;

//    lattice.verbose = 1;
//    lattice.thermalize();
//    lattice.verbose = 0;
//    for (int n = 0; n < n_data; n++) {
//
//        // update the lattice
//        if (n != 0) lattice.hmc(20);
//
//        // copy the lattice and smear it
//        su3_x_lattice smear_lattice(lattice);
////        smear_lattice.stout_smear(smear_lattice.n5_center, 0.2, 5);
//
//        // print results
//        cout << setprecision(6);
//        cout << smear_lattice.polyakov_loop(1, smear_lattice.n5_center);
//        for (int r = 2; r <= (N / 2); r++) {
//            cout << ", " << smear_lattice.polyakov_loop(r, smear_lattice.n5_center);
//        }
//        cout << endl;
//    }

    double plaq, plaq_sum;
    for (double beta = 0.5; beta <= 10.0; beta += 0.5) {
        lattice.verbose = 1;
        lattice.beta = beta;
        lattice.thermalize();
        lattice.verbose = 0;
        plaq_sum = 0.0;
        for (int n = 0; n < n_data; n++) {
            plaq = lattice.plaq(lattice.n5_center);
            plaq_sum += plaq;
            cout << setprecision(6) << fixed;
            cout << plaq << endl;
//            lattice.heat_bath(20);
            lattice.hmc(20);
        }
        cout << setprecision(2) << fixed;
        cout << "beta = " << beta;
        cout << setprecision(6) << fixed;
        cout << ", plaq = " << plaq_sum / double(n_data);
        cout << ", E = " << (1.0 - plaq_sum / double(n_data) / 3.0) << endl;
    }
    
    return 0;
}
