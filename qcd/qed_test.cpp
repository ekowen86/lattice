//
//  qed_test.cpp
//  qed_test
//
//  Created by Evan Owen on 2/16/20.
//  Copyright © 2020 Evan Owen. All rights reserved.
//

#include <iostream>
#include "qed.hpp"

int main(int argc, const char * argv[]) {
    int N = 8;
    int T = 8;
    int D = 4;
    int F = 1;
    double a = 0.5;
    double beta = 1.5;
    double m[] = {1.7, 1.7, 7.5};
    qed_lattice lattice(N, T, D, F, a, beta, m);
    lattice.parallel = true;
    lattice.thermalize(100);
    return 0;
}
