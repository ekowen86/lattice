//
//  phi4_mass.cpp
//
//  Created by Evan Owen on 11/7/19.
//  Copyright © 2019 Evan Owen. All rights reserved.
//

#include <iostream>
#include "phi4.hpp"

using namespace std;

int main(int argc, const char * argv[]) {
    
    int N = 32;
    int T = 32;
    int D = 2;
    double mu = 0.2;
    double lambda = 0.2;
    
    int n_data = 200;
    int n_sweep = 10;
    
    phi4_lattice lattice(N, T, D, mu, lambda);
    lattice.parallel = true;
    cout << "phiAve = " << lattice.phi_ave() << endl;
    lattice.thermalize();
    
    double m1 = 0.0;
    double m2 = 0.0;
    double m4 = 0.0;
    for (int n = 0; n < n_data; n++) {
        if (n) lattice.heat_bath(n_sweep);
        double m = lattice.phi_ave();
        m1 += m / double(n_data);
        m2 += m * m / double(n_data);
        m4 += m * m * m * m / double(n_data);
    }
    cout << "m1 = " << m1 << endl;
    cout << "m2 = " << m2 << endl;
    cout << "m4 = " << m4 << endl;
    cout << "U4 = " << (3.0 / 2.0) * (1.0 - (m4 / (m2 * m2 * 3.0))) << endl;

    return 0;
}
