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
    
    phi4_lattice lattice(12, 12, 4, 1.0, 0.186, 1.145);
    cout << "phiAve = " << lattice.phi_ave() << "\n";
    lattice.thermalize(1000);
    lattice.thermalize(1000);
    lattice.thermalize(1000);
    lattice.thermalize(1000);
    lattice.thermalize(1000);
    lattice.thermalize(1000);
    lattice.thermalize(1000);
    lattice.thermalize(1000);
    lattice.thermalize(1000);
    lattice.thermalize(1000);

    return 0;
}
