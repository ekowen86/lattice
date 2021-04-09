//
//  test_exp_su2.cpp
//  test the function su2_x_site::exp_su2
//
//  Created by Evan Owen on 4/5/21.
//  Copyright © 2021 Evan Owen. All rights reserved.
//

#include <iostream>
#include <unsupported/Eigen/MatrixFunctions> // for (slow) matrix exponential

#include "su2_x.hpp"

using namespace std;

int main(int argc, const char* argv[]) {

    int N = 1;
    int T = 1;
    int N5 = 1;
    int D = 1;
    double beta = 1.0;
    double eps5 = 1.0;

    // initialize a lattice with one site and one link
    su2_x_lattice lattice = su2_x_lattice(N, T, N5, D, beta, eps5);
    su2_x_site site = lattice.site[0];

    // generate a random su(2) algebra element
    site.init_momenta();

    // get an arbitrary traceless anti-hermitian matrix
    su2_link s = site.p_link[0] * I;

    // exponentiate it via exp_su2
    su2_link U = site.exp_su2(s);

    // print results
    printf("U:\n");
    cout << U << endl << endl;

    printf("U^dag:\n");
    cout << U.adjoint() << endl << endl;

    printf("U^dag U:\n");
    cout << U.adjoint() * U << endl << endl;

    printf("det(U):\n");
    cout << U.determinant() << endl << endl;

    // repeat using the Eigen matrix exponential and compare results
    su2_link U_Eigen = s.exp();

    printf("Eigen result:\n");
    cout << U_Eigen << endl << endl;

    printf("Error:\n");
    cout << U_Eigen - U << endl << endl;

    printf("norm(Error):\n");
    double error_norm = (U_Eigen - U).norm();
    cout << error_norm << endl << endl;

    printf("Test: %s\n", (error_norm < 1e-15) ? "*** Passed ***" : "*** Failed ***");

    return 0;
}
