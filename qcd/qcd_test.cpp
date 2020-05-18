//
//  qcd_test.cpp
//  qcd_test
//
//  Created by Evan Owen on 2/28/20.
//  Copyright © 2020 Evan Owen. All rights reserved.
//

#include <iostream>
#include "qcd.hpp"

int main(int argc, const char * argv[]) {
    int N = 8; cout << "N = " << N;
    int T = 4; cout << ", T = " << T;
    int D = 4; cout << ", D = " << D;
    int F = 2; cout << ", F = " << F;
    double beta = 4.0; cout << ", beta = " << beta;
    double kappa = 0.22; cout << ", kappa = " << kappa << "\n";
    double m[] = {1, 1, 5};
    qcd_lattice lattice(N, T, D, F, kappa, beta, m);
    lattice.parallel = true;
//    lattice.thermalize(1000);
    
    for (; beta < 4.501; beta += 0.05) {
        cout << "beta = " << beta << "\n";
        lattice.beta = beta;
        lattice.thermalize(1000000000);
        
        // wait until the average plaquette comes to equilibrium
        double oldP = lattice.plaquette();
        double oldF = lattice.fermion_action();
        double dP = 1.0;
        double dF = -1.0;
        int count = 0;
        while (count < 3) {
            double P = 0.0;
            double F = 0.0;
            for (int i = 0; i < 100; i++) {
                lattice.sweep();
                P += lattice.plaquette();
                F += lattice.fermion_action();
            }
            P /= 100;
            dP = P - oldP;
            oldP = P;
            F /= 100;
            dF = F - oldF;
            oldF = F;
            cout << "P = " << P;
            cout << ", dP = " << dP;
            cout << ", F = " << F;
            cout << ", dF = " << dF;
            cout << ", xi = " << lattice.xi_avg();
            cout << ", psi = " << lattice.psi_bar_psi(0) << "\n";
            if (dP < 0 && real(dF) > 0) {
                count++;
            } else {
                count = 0;
            }
        }
    }
    
    return 0;
}
