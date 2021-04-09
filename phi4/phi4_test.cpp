//
//  phi4_test.cpp
//  phi4_test
//
//  Created by Evan Owen on 11/16/20.
//

#include <iostream>
#include "phi4.hpp"

int main(int argc, const char * argv[]) {
    
    int L = 32;
    int T = 32;
    int D = 2;
    double mu2 = 0.2;
    double lambda = 0.2;
    
    int n_data = 2000;
    int n_sweep = 10;
    
    phi4_lattice lattice(L, T, D, mu2, lambda);
    lattice.parallel = false;
    lattice.n_wolff = n_sweep;
    lattice.thermalize();

    double m1 = 0.0;
    double m2 = 0.0;
    double m4 = 0.0;
    for (int n = 0; n < n_data; n++) {
        if (n) lattice.heat_bath(n_sweep);
//        if (n) lattice.wolff_update();
        double m = lattice.phi_ave();
        m1 += m;
        m2 += m * m;
        m4 += m * m * m * m;
    }
    m1 /= double(n_data);
    m2 /= double(n_data);
    m4 /= double(n_data);
    printf("m1 = %6e\n", m1);
    printf("m2 = %6e\n", m2);
    printf("m4 = %6e\n", m4);
    printf("U4 = %6e\n", (3.0 / 2.0) * (1.0 - (m4 / (m2 * m2 * 3.0))));

    return 0;
}
