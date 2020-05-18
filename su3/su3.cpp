//
//  su3.cpp
//
//  Created by Evan Owen on 10/19/18.
//  Copyright © 2018 Evan Owen. All rights reserved.
//

#include <iostream>
#include <iomanip>
#include <fstream>
#include <random>
#include <cmath>
#include <future>
#include "Eigen/Dense"
#include "su3.hpp"

using namespace std;

bool su3_site::lock() {
    if (is_locked || forward[0]->is_locked || backward[0]->is_locked) {
        this_thread::sleep_for(chrono::milliseconds(10));
        return false;
    }
    is_locked = true;
    forward[0]->is_locked = true;
    backward[0]->is_locked = true;
    return true;
}

void su3_site::unlock() {
    is_locked = false;
    forward[0]->is_locked = false;
    backward[0]->is_locked = false;
}

void su3_site::reset_links(bool cold) {
    for (int d = 0; d < lattice->D; d++) set_link(d, cold ? su3_identity : create_link(true));
}

void su3_site::copy_links(su3_site* site) {
    for (int d = 0; d < lattice->D; d++) set_link(d, su3_link(site->link[d]));
}

void su3_site::set_link(int d, su3_link value) {
    complex<double> det = value.determinant();
    if (abs(norm(det) - 1.0) > 1e-10) {
        // make the determinant 1
        double n = sqrt(norm(det));
        double a = arg(det);
        value /= cbrt(n);
        value *= polar(1.0, -a / 3.0);
    }
    link[d] = value;
    link_inverse[d] = value.transpose().conjugate();
}

su2_link su3_site::create_su2(bool random) {
    
    double a, g;
    
    if (random) {
        a = lattice->rand(-1.0, 1.0);
        g = sqrt(1 - a * a);
    } else {
        double x;
        double z = lattice->z;
        
        do {
            x = lattice->rand(exp(-2.0 * z), 1.0);
            a = 1.0 + log(x) / z;
            g = sqrt(1 - a * a);
        } while (lattice->rand() < g);
    }
    
    double theta = acos(lattice->rand(-1.0, 1.0));
    double phi = lattice->rand(0.0, 2.0 * M_PI);
    double a0 = a;
    double a1 = g * sin(theta) * cos(phi);
    double a2 = g * sin(theta) * sin(phi);
    double a3 = g * cos(theta);
    
    su2_link m;
    m(0,0) = complex<double>(a0, a3);
    m(0,1) = complex<double>(a2, a1);
    m(1,0) = complex<double>(-a2, a1);
    m(1,1) = complex<double>(a0, -a3);
    return m;
}

su3_link su3_site::create_link(bool random) {
    
    su3_link R = su3_identity;
    R.block<2,2>(0,0) = create_su2(random);
    
    su3_link S = su3_identity;
    su2_link s = create_su2(random);
    S(0,0) = s(0,0);
    S(0,2) = s(0,1);
    S(2,0) = s(1,0);
    S(2,2) = s(1,1);
    
    su3_link T = su3_identity;
    T.block<2,2>(1,1) = create_su2(random);
    
    su3_link RST = R * S * T;
    if (rand() > 0.5) return RST.transpose().conjugate();
    return RST;
}

double su3_site::plaq() {
    // compute the average plaquette at this site
    double U = 0;
    su3_link u;
    for (int d1 = 0; d1 < lattice->D; d1++) {
        for (int d2 = d1 + 1; d2 < lattice->D; d2++) {
            u = link[d1];
            u *= forward[d1]->link[d2];
            u *= forward[d2]->link_inverse[d1];
            u *= link_inverse[d2];
            U += u.trace().real() / 3;
        }
    }
    return U / lattice->D / (lattice->D - 1) * 2;
}

double su3_site::wilson_loop(int a, int b) {
    // compute the average a x b wilson loop at this site
    double U = 0;
    su3_link u;
    su3_site* s = this;
    int x;
    for (int d1 = 0; d1 < lattice->D; d1++) {
        for (int d2 = d1 + 1; d2 < lattice->D; d2++) {
            u = su3_identity;
            for (x = 0; x < a; x++) {
                u *= s->link[d1];
                s = s->forward[d1];
            }
            for (x = 0; x < b; x++) {
                u *= s->link[d2];
                s = s->forward[d2];
            }
            for (x = 0; x < a; x++) {
                s = s->backward[d1];
                u *= s->link_inverse[d1];
            }
            for (x = 0; x < b; x++) {
                s = s->backward[d2];
                u *= s->link_inverse[d2];
            }
            U += u.trace().real() / 3;
        }
    }
    return U / lattice->D / (lattice->D - 1);
}

double su3_site::polyakov_loop(int r) {
    // compute the average product of polyakov loops at this site
    su3_site* s = this;
    int x;
    
    // polyakov loop starting at this site
    su3_link u1 = su3_identity;
    for (x = 0; x < lattice->T; x++) {
        u1 *= s->link[0];
        s = s->forward[0];
    }
    
    // multiply by adjacent polyakov loops a distance r away
    double U = 0;
    su3_link u2;
    for (int d = 1; d < lattice->D; d++) {
        u2 = su3_identity;
        s = this;
        for (x = 0; x < r; x++) s = s->forward[d];
        
        for (x = 0; x < lattice->T; x++) {
            s = s->backward[0];
            u2 *= s->link_inverse[0];
        }
        U += (u1 * u2).trace().real() / 3;
    }
    return U / (lattice->D - 1);
}

double su3_site::correlator(int T) {
    // correlator of duration T
    su3_link u = su3_identity;
    su3_site* s = this;
    int x;
    for (x = 0; x < T; x++) {
        u *= s->link[0];
        s = s->forward[0];  // time direction is always 0
    }
    return u.trace().real() / 3;
}

double su3_site::mag_U() {
    return this->link[0].trace().real() / 3;
}

double su3_site::abs_U() {
    return abs(this->link[0].trace().real()) / 3;
}

double su3_site::four_point(int T, int R) {
    // four-point function of duration T and spacing R
    double U = 0;
    su3_link u;
    su3_site* s = this;
    int x;
    for (int d = 1; d < lattice->D; d++) { // d is the spacing direction
        u = su3_identity;
        for (x = 0; x < T; x++) {
            u *= s->link[0];
            s = s->forward[0];  // time direction is always 0
        }
        for (x = 0; x < R; x++) s = s->forward[d];
        for (x = 0; x < T; x++) {
            s = s->backward[0];  // time direction is always 0
            u *= s->link_inverse[0];
        }
        for (x = 0; x < R; x++) s = s->backward[d];
        U += u.trace().real() / 3;
    }
    return U / (lattice->D - 1);
}

double su3_site::contract_24(su3_link* A) {
    return contract_3(contract_8(A[0], A[1], A[2], A[3], A[4], A[5], A[6], A[7]),
                      contract_8(A[8], A[9], A[10], A[11], A[12], A[13], A[14], A[15]),
                      contract_8(A[16], A[17], A[18], A[19], A[20], A[21], A[22], A[23]));
}

double su3_site::contract_12(su3_link* A) {
    return contract_3(contract_4(A[0], A[1], A[2], A[3]), contract_4(A[4], A[5], A[6], A[7]), contract_4(A[8], A[9], A[10], A[11]));
}

su3_link su3_site::contract_8(su3_link A, su3_link B, su3_link C, su3_link D, su3_link E, su3_link F, su3_link G, su3_link H) {
    return contract_2(contract_4(A, B, C, D), contract_4(E, F, G, H));
}

double su3_site::contract_6(su3_link A, su3_link B, su3_link C, su3_link D, su3_link E, su3_link F) {
    return contract_3(contract_2(A, B), contract_2(C, D), contract_2(E, F));
}

su3_link su3_site::contract_4(su3_link A, su3_link B, su3_link C, su3_link D) {
    return contract_2(contract_2(A, B), contract_2(C, D));
}

double su3_site::contract_3(su3_link A, su3_link B, su3_link C) {
    complex<double> U = 0.0;
    int a1, a2, a3, b1, b2, b3;
    for (a1 = 0; a1 < 3; a1++) {
        a2 = (a1+1) % 3;
        a3 = (a1+2) % 3;
        for (b1 = 0; b1 < 3; b1++) {
            b2 = (b1+1) % 3;
            b3 = (b1+2) % 3;
            U += A(a1,b1) * B(a2,b2) * C(a3,b3);
            U -= A(a1,b1) * B(a2,b3) * C(a3,b2);
            U -= A(a1,b1) * B(a3,b2) * C(a2,b3);
            U += A(a1,b1) * B(a3,b3) * C(a2,b2);
        }
    }
    return U.real() / 6.0;
}

su3_link su3_site::contract_2(su3_link B, su3_link C) {
    su3_link A = su3_zero;
    int a1, a2, a3, b1, b2, b3;
    for (a1 = 0; a1 < 3; a1++) {
        a2 = (a1+1) % 3;
        a3 = (a1+2) % 3;
        for (b1 = 0; b1 < 3; b1++) {
            b2 = (b1+1) % 3;
            b3 = (b1+2) % 3;
            A(a1,b1) += B(a2,b2) * C(a3,b3);
            A(a1,b1) -= B(a2,b3) * C(a3,b2);
            A(a1,b1) -= B(a3,b2) * C(a2,b3);
            A(a1,b1) += B(a3,b3) * C(a2,b2);
        }
    }
    return A / 2.0;
}

double su3_site::contract_1(su3_link A, su3_link B) {
    complex<double> U = 0.0;
    for (int a = 0; a < 3; a++) {
        for (int b = 0; b < 3; b++) {
            U += A(a,b) * B(a,b);
        }
    }
    return U.real() / 3;
}

// A-----B-----C  +-----> d1

void su3_site::quark_0(double* U) {
    double U0 = 0.0;
    for (int d1 = 1; d1 < lattice->D; d1++) {
        su3_site* A = this;
        su3_site* B = A->forward[d1];
        su3_site* C = B->forward[d1];
        
        su3_link a = Q(A);
        su3_link b = Q(B);
        su3_link c = Q(C);

        U0 += contract_3(a, b, c);
    }
    U[0] += U0 / (lattice->D - 1);
}

// A2----B2----C2   ^ d2
//  |     |     |   |
//  |     |     |   |
// A1----B1----C1   +-----> d1

void su3_site::quark_1(double* U) {
    double U1 = 0.0;
    double U2 = 0.0;
    double U3 = 0.0;
    double U4 = 0.0;
    su3_link A, B, C;
    for (int d1 = 1; d1 < lattice->D; d1++) {
        for (int d2 = 1; d2 < lattice->D; d2++) {
            if (d1 == d2) continue;
            su3_link A1 = link[0];
            su3_link B1 = forward[d1]->link[0];
            su3_link C1 = forward[d1]->forward[d1]->link[0];
            su3_link A2 = forward[d2]->link[0];
            su3_link B2 = forward[d2]->forward[d1]->link[0];
            su3_link C2 = forward[d2]->forward[d1]->forward[d1]->link[0];
            
            A = contract_2(A1, A2);
            B = contract_2(B1, B2);
            C = contract_2(C1, C2);
            U1 += contract_3(A, B, C);
            
            A = contract_2(A1, A2);
            B = contract_2(B1, C1);
            C = contract_2(B2, C2);
            U2 += contract_3(A, B, C);
            A = contract_2(A1, B1);
            B = contract_2(A2, B2);
            C = contract_2(C1, C2);
            U2 += contract_3(A, B, C);

            U3 += contract_3(A1, B1, C1) * contract_3(A2, B2, C2);
            
            U4 += contract_3(A1, A2, B1) * contract_3(B2, C1, C2);
            U4 += contract_3(A1, A2, B2) * contract_3(B1, C1, C2);
        }
    }
    U[0] += U1 / (lattice->D - 1) / (lattice->D - 2);
    U[1] += U2 / (lattice->D - 1) / (lattice->D - 2) / 2;
    U[2] += U3 / (lattice->D - 1) / (lattice->D - 2);
    U[3] += U4 / (lattice->D - 1) / (lattice->D - 2) / 2;
}

void su3_site::quark_2(double* U) {
    double U5 = 0.0;
    double U6 = 0.0;
    double U7 = 0.0;
    double U8 = 0.0;
    double U9 = 0.0;
    double U10 = 0.0;
    double U11 = 0.0;
    double U12 = 0.0;
    double U13 = 0.0;
    su3_link A, B, C;
    for (int d1 = 1; d1 < lattice->D; d1++) {

        su3_site* A = this;
        su3_site* B = A->forward[d1];
        su3_site* C = B->forward[d1];

        for (int d2 = 1; d2 < lattice->D; d2++) {
            if (d1 == d2) continue;
            
            su3_site* A1 = A;
            su3_site* B1 = B;
            su3_site* C1 = C;

            for (int d3 = 1; d3 < lattice->D; d3++) {
                if (d2 == d3) continue;
                if (d1 == d3) continue;
                
                su3_site* A11 = A1; su3_link a11 = Q(A11);
                su3_site* A12 = A11->forward[d2]; su3_link a12 = Q(A12);
                su3_site* A21 = A11->forward[d3]; su3_link a21 = Q(A21);
                su3_site* A22 = A21->forward[d2]; su3_link a22 = Q(A22);
                
                su3_site* B11 = B1; su3_link b11 = Q(B11);
                su3_site* B12 = B11->forward[d2]; su3_link b12 = Q(B12);
                su3_site* B21 = B11->forward[d3]; su3_link b21 = Q(B21);
                su3_site* B22 = B21->forward[d2]; su3_link b22 = Q(B22);
                
                su3_site* C11 = C1; su3_link c11 = Q(C11);
                su3_site* C12 = C11->forward[d2]; su3_link c12 = Q(C12);
                su3_site* C21 = C11->forward[d3]; su3_link c21 = Q(C21);
                su3_site* C22 = C21->forward[d2]; su3_link c22 = Q(C22);

                // vertically contracted
                su3_link A1 = contract_2(a11, a12);
                su3_link A2 = contract_2(a21, a22);
                su3_link B1 = contract_2(b11, b12);
                su3_link B2 = contract_2(b21, b22);
                su3_link C1 = contract_2(c11, c12);
                su3_link C2 = contract_2(c21, c22);

                // horizontally contracted
                su3_link AA1 = contract_2(a11, a21);
                su3_link AA2 = contract_2(a12, a22);
//                su3_link BB1 = contract_2(B11, B21);
//                su3_link BB2 = contract_2(B12, B22);
                su3_link CC1 = contract_2(c11, c21);
                su3_link CC2 = contract_2(c12, c22);
                su3_link AB11 = contract_2(a11, b11);
                su3_link AB12 = contract_2(a12, b12);
                su3_link AB21 = contract_2(a21, b21);
                su3_link AB22 = contract_2(a22, b22);
                su3_link BC11 = contract_2(b11, c11);
                su3_link BC12 = contract_2(b12, c12);
                su3_link BC21 = contract_2(b21, c21);
                su3_link BC22 = contract_2(b22, c22);

                //    A12---B12---C12   ^ d2
                //   / |   / |   / |    |
                // A22---B22---C22 |    |
                //  | A11-|-B11-|-C11   +-----> d1
                //  | /   | /   | /    /
                // A21---B21---C21   L d3

                U5 += contract_6(A1, A2, B1, B2, C1, C2);

                U6 += contract_6(A1, B1, A2, B2, C1, C2);
                U6 += contract_6(A1, A2, B1, C1, B2, C2);

                U7 += contract_3(A1, B1, C1) * contract_3(A2, B2, C2);

                U8 += contract_3(A1, A2, B1) * contract_3(B2, C1, C2);
                U8 += contract_3(A1, A2, B2) * contract_3(B1, C1, C2);
                
                double U9_1 = contract_3(A1, BC11, BC12);
                double U9_2 = contract_3(A2, BC21, BC22);
                double U9_3 = contract_3(C1, AB11, AB12);
                double U9_4 = contract_3(C2, AB21, AB22);
                U9 += U9_1 * U9_2;
                U9 += U9_3 * U9_4;
                U9 += U9_1 * U9_4;
                U9 += U9_2 * U9_3;

                double U10_1 = contract_3(A1, AB21, AB22);
                double U10_2 = contract_3(C2, BC12, BC11);
                double U10_3 = contract_3(B2, AA1, AA2);
                double U10_4 = contract_3(B1, CC1, CC2);
                U10 += U10_1 * U10_2;
                U10 += U10_3 * U10_4;
                U10 += U10_1 * U10_4;
                U10 += U10_2 * U10_3;
                U10_1 = contract_3(A2, AB11, AB12);
                U10_2 = contract_3(C1, BC22, BC21);
                U10_3 = contract_3(B1, AA1, AA2);
                U10_4 = contract_3(B2, CC1, CC2);
                U10 += U10_1 * U10_2;
                U10 += U10_3 * U10_4;
                U10 += U10_1 * U10_4;
                U10 += U10_2 * U10_3;

                U11 +=
                contract_3(a11, b11, c11) * contract_3(a12, b12, c12) *
                contract_3(a21, b21, c21) * contract_3(a22, b22, c22);
                
                double U12_1 = contract_3(a11, a12, b12);
                double U12_2 = contract_3(b11, c11, c12);
                double U12_3 = contract_3(a21, a22, b22);
                double U12_4 = contract_3(b21, c21, c22);
                double U12_5 = contract_3(a11, a12, b11);
                double U12_6 = contract_3(b12, c11, c12);
                double U12_7 = contract_3(a21, a22, b21);
                double U12_8 = contract_3(b22, c21, c22);
                U12 += U12_1 * U12_2 * U12_3 * U12_4;
                U12 += U12_5 * U12_6 * U12_7 * U12_8;
                U12 += U12_1 * U12_2 * U12_7 * U12_8;
                U12 += U12_5 * U12_6 * U12_3 * U12_4;
                
                U13 +=
                contract_3(a11, a21, a22) * contract_3(a12, b12, b11) *
                contract_3(b22, b21, c21) * contract_3(c22, c12, c11);
                U13 +=
                contract_3(a21, a11, a12) * contract_3(a22, b22, b21) *
                contract_3(b12, b11, c11) * contract_3(c12, c22, c21);
                U13 +=
                contract_3(c21, c11, c12) * contract_3(c22, b22, b21) *
                contract_3(b12, b11, a11) * contract_3(a12, a22, a21);
                U13 +=
                contract_3(c21, c11, c12) * contract_3(c22, b22, b21) *
                contract_3(b12, b11, a11) * contract_3(a12, a22, a21);
            }
        }
    }
    U[0] += U5 / (lattice->D - 1) / (lattice->D - 2) / (lattice->D - 3);
    U[1] += U6 / (lattice->D - 1) / (lattice->D - 2) / (lattice->D - 3) / 2;
    U[2] += U7 / (lattice->D - 1) / (lattice->D - 2) / (lattice->D - 3);
    U[3] += U8 / (lattice->D - 1) / (lattice->D - 2) / (lattice->D - 3) / 2;
    U[4] += U9 / (lattice->D - 1) / (lattice->D - 2) / (lattice->D - 3) / 4;
    U[5] += U10 / (lattice->D - 1) / (lattice->D - 2) / (lattice->D - 3) / 8;
    U[6] += U11 / (lattice->D - 1) / (lattice->D - 2) / (lattice->D - 3);
    U[7] += U12 / (lattice->D - 1) / (lattice->D - 2) / (lattice->D - 3) / 4;
    U[8] += U13 / (lattice->D - 1) / (lattice->D - 2) / (lattice->D - 3) / 4;
}

//    A14---B14---C14
//   / |   / |   / |
// A24---B24---C24 |
//  | A13-|-B13-|-C13
//  |/ |  |/ |  |/ |
// A23---B23---C23 |
//  | A12-|-B12-|-C12   ^ d2
//  |/ |  |/ |  |/ |    |
// A22---B22---C22 |    |
//  | A11-|-B11-|-C11   +-----> d1
//  |/    |/    |/     /
// A21---B21---C21    L d3

void su3_site::quark_3(double* U) {
    double U14 = 0.0;
    double U15 = 0.0;
    double U16 = 0.0;
    double U17 = 0.0;
    double U18 = 0.0;
    double U19 = 0.0;
    double U20 = 0.0;
    double U21 = 0.0;
    double U22 = 0.0;
    double U23 = 0.0;
    double U24 = 0.0;
    double U25 = 0.0;
    double U26 = 0.0;
    for (int d1 = 1; d1 < lattice->D; d1++) {
        for (int d2 = 1; d2 < lattice->D; d2++) {
            if (d1 == d2) continue;
            for (int d3 = 1; d3 < lattice->D; d3++) {
                if (d2 == d3) continue;
                if (d1 == d3) continue;
                
                su3_site* A11 = this; su3_link a11 = Q(A11);
                su3_site* A12 = A11->forward[d2]; su3_link a12 = Q(A12);
                su3_site* A13 = A12->forward[d2]; su3_link a13 = Q(A13);
                su3_site* A14 = A13->forward[d2]; su3_link a14 = Q(A14);
                su3_site* A21 = A11->forward[d3]; su3_link a21 = Q(A21);
                su3_site* A22 = A21->forward[d2]; su3_link a22 = Q(A22);
                su3_site* A23 = A22->forward[d2]; su3_link a23 = Q(A23);
                su3_site* A24 = A23->forward[d2]; su3_link a24 = Q(A24);
                
                su3_site* B11 = A11->forward[d1]; su3_link b11 = Q(B11);
                su3_site* B12 = B11->forward[d2]; su3_link b12 = Q(B12);
                su3_site* B13 = B12->forward[d2]; su3_link b13 = Q(B13);
                su3_site* B14 = B13->forward[d2]; su3_link b14 = Q(B14);
                su3_site* B21 = B11->forward[d3]; su3_link b21 = Q(B21);
                su3_site* B22 = B21->forward[d2]; su3_link b22 = Q(B22);
                su3_site* B23 = B22->forward[d2]; su3_link b23 = Q(B23);
                su3_site* B24 = B23->forward[d2]; su3_link b24 = Q(B24);
                
                su3_site* C11 = B11->forward[d1]; su3_link c11 = Q(C11);
                su3_site* C12 = C11->forward[d2]; su3_link c12 = Q(C12);
                su3_site* C13 = C12->forward[d2]; su3_link c13 = Q(C13);
                su3_site* C14 = C13->forward[d2]; su3_link c14 = Q(C14);
                su3_site* C21 = C11->forward[d3]; su3_link c21 = Q(C21);
                su3_site* C22 = C21->forward[d2]; su3_link c22 = Q(C22);
                su3_site* C23 = C22->forward[d2]; su3_link c23 = Q(C23);
                su3_site* C24 = C23->forward[d2]; su3_link c24 = Q(C24);

                U14 += contract_3(contract_8(a11, a12, a21, a22, b11, b21, c11, c21),
                                  contract_8(a13, b13, a23, b23, b12, c12, b22, c22),
                                  contract_8(a14, b14, a24, b24, c13, c14, c23, c24));
                U14 += contract_3(contract_8(c11, c12, c21, c22, b11, b21, a11, a21),
                                  contract_8(c13, b13, c23, b23, b12, a12, b22, a22),
                                  contract_8(c14, b14, c24, b24, a13, a14, a23, a24));

                U15 += contract_3(contract_8(a21, a22, b21, b22, c21, c22, c11, c12),
                                  contract_8(a11, a12, b11, b12, a23, a13, b23, b13),
                                  contract_8(a24, a14, b24, b14, c23, c24, c13, c14));
                U15 += contract_3(contract_8(a11, a12, b11, b12, c11, c12, c21, c22),
                                  contract_8(a21, a22, b21, b22, a13, a23, b13, b23),
                                  contract_8(a14, a24, b14, b24, c13, c14, c23, c24));
                U15 += contract_3(contract_8(c21, c22, b21, b22, a21, a22, a11, a12),
                                  contract_8(c11, c12, b11, b12, c23, c13, b23, b13),
                                  contract_8(c24, c14, b24, b14, a23, a24, a13, a14));
                U15 += contract_3(contract_8(c11, c12, b11, b12, a11, a12, a21, a22),
                                  contract_8(c21, c22, b21, b22, c13, c23, b13, b23),
                                  contract_8(c14, c24, b14, b24, a13, a14, a23, a24));

                U16 +=
                contract_3(a11, b11, c11) * contract_3(a21, b21, c21) *
                contract_3(a12, b12, c12) * contract_3(a22, b22, c22) *
                contract_3(a13, b13, c13) * contract_3(a23, b23, c23) *
                contract_3(a14, b14, c14) * contract_3(a24, b24, c24);
                
                U17 +=
                contract_3(a11, a21, a22) * contract_3(b11, b21, c21) *
                contract_3(c22, c12, c11) * contract_3(b22, b23, c23) *
                contract_3(a12, b12, b13) * contract_3(a23, a13, a14) *
                contract_3(a24, b24, b14) * contract_3(c24, c14, c13);
                U17 +=
                contract_3(a21, a11, a12) * contract_3(b21, b11, c11) *
                contract_3(c12, c22, c21) * contract_3(b12, b13, c13) *
                contract_3(a22, b22, b23) * contract_3(a13, a23, a24) *
                contract_3(a14, b14, b24) * contract_3(c14, c24, c23);
                U17 +=
                contract_3(c11, c21, c22) * contract_3(b11, b21, a21) *
                contract_3(a22, a12, a11) * contract_3(b22, b23, a23) *
                contract_3(c12, b12, b13) * contract_3(c23, c13, c14) *
                contract_3(c24, b24, b14) * contract_3(a24, a14, a13);
                U17 +=
                contract_3(c21, c11, c12) * contract_3(b21, b11, a11) *
                contract_3(a12, a22, a21) * contract_3(b12, b13, a13) *
                contract_3(c22, b22, b23) * contract_3(c13, c23, c24) *
                contract_3(c14, b14, b24) * contract_3(a14, a24, a23);
                
                U18 +=
                contract_3(a21, b21, b22) * contract_3(b11, c11, c21) *
                contract_3(a11, a12, b12) * contract_3(a22, a23, a13) *
                contract_3(c22, c12, c13) * contract_3(b23, c23, c24) *
                contract_3(a14, a24, b24) * contract_3(b13, b14, c14);
                U18 +=
                contract_3(a11, b11, b12) * contract_3(b21, c21, c11) *
                contract_3(a21, a22, b22) * contract_3(a12, a13, a23) *
                contract_3(c12, c22, c23) * contract_3(b13, c13, c14) *
                contract_3(a24, a14, b14) * contract_3(b23, b24, c24);
                U18 +=
                contract_3(c21, b21, b22) * contract_3(b11, a11, a21) *
                contract_3(c11, c12, b12) * contract_3(c22, c23, c13) *
                contract_3(a22, a12, a13) * contract_3(b23, a23, a24) *
                contract_3(c14, c24, b24) * contract_3(b13, b14, a14);
                U18 +=
                contract_3(c11, b11, b12) * contract_3(b21, a21, a11) *
                contract_3(c21, c22, b22) * contract_3(c12, c13, c23) *
                contract_3(a12, a22, a23) * contract_3(b13, a13, a14) *
                contract_3(c24, c14, b14) * contract_3(b23, b24, a24);
                
                U19 += contract_3(contract_8(a11, a12, a21, a22, a13, a14, a23, a24),
                                  contract_8(b11, b12, b21, b22, b13, b14, b23, b24),
                                  contract_8(c11, c12, c21, c22, c13, c14, c23, c24));

                U20 +=
                contract_3(a11, a21, b21) * contract_3(b11, c11, c21) *
                contract_3(a12, a22, b22) * contract_3(b12, c12, c22) *
                contract_3(a13, a23, b23) * contract_3(b13, c13, c23) *
                contract_3(a14, a24, b24) * contract_3(b14, c14, c24);
                U20 +=
                contract_3(a11, a21, b11) * contract_3(b21, c11, c21) *
                contract_3(a12, a22, b12) * contract_3(b22, c12, c22) *
                contract_3(a13, a23, b13) * contract_3(b23, c13, c23) *
                contract_3(a14, a24, b14) * contract_3(b24, c14, c24);

                U21 +=
                contract_3(a11, a21, b21) * contract_3(b11, c11, c21) *
                contract_3(a12, a22, b12) * contract_3(b22, c12, c22) *
                contract_3(a13, a23, b23) * contract_3(b13, c13, c23) *
                contract_3(a14, a24, b14) * contract_3(b24, c14, c24);
                U21 +=
                contract_3(a11, a21, b11) * contract_3(b21, c11, c21) *
                contract_3(a12, a22, b22) * contract_3(b12, c12, c22) *
                contract_3(a13, a23, b13) * contract_3(b23, c13, c23) *
                contract_3(a14, a24, b24) * contract_3(b14, c14, c24);
                
                U22 += contract_3(contract_8(a11, b11, a21, b21, c11, c12, c21, c22),
                                  contract_8(a12, a22, b12, b22, a13, a23, b13, b23),
                                  contract_8(a14, a24, b14, b24, c13, c23, c14, c24));
                U22 += contract_3(contract_8(c11, b11, c21, b21, a11, a12, a21, a22),
                                  contract_8(c12, c22, b12, b22, c13, c23, b13, b23),
                                  contract_8(c14, c24, b14, b24, a13, a23, a14, a24));

                U23 += contract_3(contract_8(a11, b11, a21, b21, a12, b12, a22, b22),
                                  contract_8(a13, b13, a23, b23, a14, b14, a24, b24),
                                  contract_8(c11, c21, c12, c22, c13, c23, c14, c24));
                U23 += contract_3(contract_8(c11, b11, c21, b21, c12, b12, c22, b22),
                                  contract_8(c13, b13, c23, b23, c14, b14, c24, b24),
                                  contract_8(a11, a21, a12, a22, a13, a23, a14, a24));
                
                U24 += contract_3(contract_8(a21, a22, a23, a24, b22, b23, a12, a13),
                                  contract_8(a11, b11, b12, c12, b21, c21, c11, c22),
                                  contract_8(a14, b14, b24, c24, b13, c13, c14, c23));
                U24 += contract_3(contract_8(a11, a12, a13, a14, b12, b13, a22, a23),
                                  contract_8(a21, b21, b22, c22, b11, c11, c21, c12),
                                  contract_8(a24, b24, b14, c14, b23, c23, c24, c13));
                U24 += contract_3(contract_8(c21, c22, c23, c24, b22, b23, c12, c13),
                                  contract_8(c11, b11, b12, a12, b21, a21, a11, a22),
                                  contract_8(c14, b14, b24, a24, b13, a13, a14, a23));
                U24 += contract_3(contract_8(c11, c12, c13, c14, b12, b13, c22, c23),
                                  contract_8(c21, b21, b22, a22, b11, a11, a21, a12),
                                  contract_8(c24, b24, b14, a14, b23, a23, a24, a13));

                U25 += contract_3(contract_8(a11, a21, b11, b21, a12, a22, b12, b22),
                                  contract_8(c11, c21, c12, c22, c13, c23, b13, b23),
                                  contract_8(a13, a14, a23, a24, b14, b24, c14, c24));
                U25 += contract_3(contract_8(c11, c21, b11, b21, c12, c22, b12, b22),
                                  contract_8(a11, a21, a12, a22, a13, a23, b13, b23),
                                  contract_8(c13, c14, c23, c24, b14, b24, a14, a24));
                U25 += contract_3(contract_8(a14, a24, b14, b24, a13, a23, b13, b23),
                                  contract_8(c14, c24, c13, c23, c12, c22, b12, b22),
                                  contract_8(a12, a11, a22, a21, b11, b21, c11, c21));
                U25 += contract_3(contract_8(c14, c24, b14, b24, c13, c23, b13, b23),
                                  contract_8(a14, a24, a13, a23, a12, a22, b12, b22),
                                  contract_8(c12, c11, c22, c21, b11, b21, a11, a21));
                
                U26 += contract_3(contract_8(a11, a21, b11, b21, a12, a22, b12, b22),
                                  contract_8(c11, c21, c12, c22, c13, c23, b13, b23),
                                  contract_8(a13, a14, a23, a24, b14, b24, c14, c24));
                U26 += contract_3(contract_8(c11, c21, b11, b21, c12, c22, b12, b22),
                                  contract_8(a11, a21, a12, a22, a13, a23, b13, b23),
                                  contract_8(c13, c14, c23, c24, b14, b24, a14, a24));
                U26 += contract_3(contract_8(a14, a24, b14, b24, a13, a23, b13, b23),
                                  contract_8(c14, c24, c13, c23, c12, c22, b12, b22),
                                  contract_8(a12, a11, a22, a21, b11, b21, c11, c21));
                U26 += contract_3(contract_8(c14, c24, b14, b24, c13, c23, b13, b23),
                                  contract_8(a14, a24, a13, a23, a12, a22, b12, b22),
                                  contract_8(c12, c11, c22, c21, b11, b21, a11, a21));
            }
        }
    }
    U[0] += U14 / (lattice->D - 1) / (lattice->D - 2) / (lattice->D - 3) / 2;
    U[1] += U15 / (lattice->D - 1) / (lattice->D - 2) / (lattice->D - 3) / 4;
    U[2] += U16 / (lattice->D - 1) / (lattice->D - 2) / (lattice->D - 3);
    U[3] += U17 / (lattice->D - 1) / (lattice->D - 2) / (lattice->D - 3) / 4;
    U[4] += U18 / (lattice->D - 1) / (lattice->D - 2) / (lattice->D - 3) / 4;
    U[5] += U19 / (lattice->D - 1) / (lattice->D - 2) / (lattice->D - 3);
    U[6] += U20 / (lattice->D - 1) / (lattice->D - 2) / (lattice->D - 3) / 2;
    U[7] += U21 / (lattice->D - 1) / (lattice->D - 2) / (lattice->D - 3) / 2;
    U[8] += U22 / (lattice->D - 1) / (lattice->D - 2) / (lattice->D - 3) / 2;
    U[9] += U23 / (lattice->D - 1) / (lattice->D - 2) / (lattice->D - 3) / 2;
    U[10] += U24 / (lattice->D - 1) / (lattice->D - 2) / (lattice->D - 3) / 4;
    U[11] += U25 / (lattice->D - 1) / (lattice->D - 2) / (lattice->D - 3) / 4;
    U[12] += U26 / (lattice->D - 1) / (lattice->D - 2) / (lattice->D - 3) / 4;
}

//          A14---B14---C14
//         / |   / |   / |
//       A24---B24---C24 |
//      / |   / |   / | C13
//    A34---B34---C34 |/ |
//   / |   / |   / | C23 |
// A44---B44---C44 |/ | C12   ^ d2
//  | A23-|-B33-|-C33 |/ |    |
//  |/ |  |/ |  |/ | C22 |    |
// A43---B43---C43 |/ | C11   +-----> d1
//  | A22-|-B32-|-C32 |/     /
//  |/ |  |/ |  |/ | C21    /
// A42---B42---C42 |/      L d3
//  | A21-|-B31-|-C31
//  |/    |/    |/
// A41---B41---C41

void su3_site::quark_4(double* U) {
    double U30 = 0.0;
    double U31 = 0.0;
    double U32 = 0.0;
    double U33 = 0.0;
    for (int d1 = 1; d1 < lattice->D; d1++) {
        for (int d2 = 1; d2 < lattice->D; d2++) {
            if (d1 == d2) continue;
            for (int d3 = 1; d3 < lattice->D; d3++) {
                if (d2 == d3) continue;
                if (d1 == d3) continue;
                
                su3_site* A11 = this; su3_link a11 = Q(A11);
                su3_site* A12 = A11->forward[d2]; su3_link a12 = Q(A12);
                su3_site* A13 = A12->forward[d2]; su3_link a13 = Q(A13);
                su3_site* A14 = A13->forward[d2]; su3_link a14 = Q(A14);
                su3_site* A21 = A11->forward[d3]; su3_link a21 = Q(A21);
                su3_site* A22 = A21->forward[d2]; su3_link a22 = Q(A22);
                su3_site* A23 = A22->forward[d2]; su3_link a23 = Q(A23);
                su3_site* A24 = A23->forward[d2]; su3_link a24 = Q(A24);
                su3_site* A31 = A21->forward[d3]; su3_link a31 = Q(A31);
                su3_site* A32 = A31->forward[d2]; su3_link a32 = Q(A32);
                su3_site* A33 = A32->forward[d2]; su3_link a33 = Q(A33);
                su3_site* A34 = A33->forward[d2]; su3_link a34 = Q(A34);
                su3_site* A41 = A31->forward[d3]; su3_link a41 = Q(A41);
                su3_site* A42 = A41->forward[d2]; su3_link a42 = Q(A42);
                su3_site* A43 = A42->forward[d2]; su3_link a43 = Q(A43);
                su3_site* A44 = A43->forward[d2]; su3_link a44 = Q(A44);

                su3_site* B11 = A11->forward[d1]; su3_link b11 = Q(B11);
                su3_site* B12 = B11->forward[d2]; su3_link b12 = Q(B12);
                su3_site* B13 = B12->forward[d2]; su3_link b13 = Q(B13);
                su3_site* B14 = B13->forward[d2]; su3_link b14 = Q(B14);
                su3_site* B21 = B11->forward[d3]; su3_link b21 = Q(B21);
                su3_site* B22 = B21->forward[d2]; su3_link b22 = Q(B22);
                su3_site* B23 = B22->forward[d2]; su3_link b23 = Q(B23);
                su3_site* B24 = B23->forward[d2]; su3_link b24 = Q(B24);
                su3_site* B31 = B21->forward[d3]; su3_link b31 = Q(B31);
                su3_site* B32 = B31->forward[d2]; su3_link b32 = Q(B32);
                su3_site* B33 = B32->forward[d2]; su3_link b33 = Q(B33);
                su3_site* B34 = B33->forward[d2]; su3_link b34 = Q(B34);
                su3_site* B41 = B31->forward[d3]; su3_link b41 = Q(B41);
                su3_site* B42 = B41->forward[d2]; su3_link b42 = Q(B42);
                su3_site* B43 = B42->forward[d2]; su3_link b43 = Q(B43);
                su3_site* B44 = B43->forward[d2]; su3_link b44 = Q(B44);

                su3_site* C11 = B11->forward[d1]; su3_link c11 = Q(C11);
                su3_site* C12 = C11->forward[d2]; su3_link c12 = Q(C12);
                su3_site* C13 = C12->forward[d2]; su3_link c13 = Q(C13);
                su3_site* C14 = C13->forward[d2]; su3_link c14 = Q(C14);
                su3_site* C21 = C11->forward[d3]; su3_link c21 = Q(C21);
                su3_site* C22 = C21->forward[d2]; su3_link c22 = Q(C22);
                su3_site* C23 = C22->forward[d2]; su3_link c23 = Q(C23);
                su3_site* C24 = C23->forward[d2]; su3_link c24 = Q(C24);
                su3_site* C31 = C21->forward[d3]; su3_link c31 = Q(C31);
                su3_site* C32 = C31->forward[d2]; su3_link c32 = Q(C32);
                su3_site* C33 = C32->forward[d2]; su3_link c33 = Q(C33);
                su3_site* C34 = C33->forward[d2]; su3_link c34 = Q(C34);
                su3_site* C41 = C31->forward[d3]; su3_link c41 = Q(C41);
                su3_site* C42 = C41->forward[d2]; su3_link c42 = Q(C42);
                su3_site* C43 = C42->forward[d2]; su3_link c43 = Q(C43);
                su3_site* C44 = C43->forward[d2]; su3_link c44 = Q(C44);

                U30 += contract_6(contract_8(a11, a12, a21, a22, a13, a14, a23, a24),
                                  contract_8(a31, a32, a41, a42, a33, a34, a43, a44),
                                  contract_8(b11, b12, b21, b22, b13, b14, b23, b24),
                                  contract_8(b31, b32, b41, b42, b33, b34, b43, b44),
                                  contract_8(c11, c12, c21, c22, c13, c14, c23, c24),
                                  contract_8(c31, c32, c41, c42, c33, c34, c43, c44));
                
                U31 += contract_6(contract_8(a11, a12, a21, a22, b11, b12, b21, b22),
                                  contract_8(a31, a32, a41, a42, b31, b32, b41, b42),
                                  contract_8(a13, a14, a23, a24, b13, b14, b23, b24),
                                  contract_8(a33, a34, a43, a44, b33, b34, b44, b44),
                                  contract_8(c11, c12, c21, c22, c13, c14, c23, c24),
                                  contract_8(c31, c32, c41, c42, c33, c34, c43, c44));
                U31 += contract_6(contract_8(c11, c12, c21, c22, b11, b12, b21, b22),
                                  contract_8(c31, c32, c41, c42, b31, b32, b41, b42),
                                  contract_8(c13, c14, c23, c24, b13, b14, b23, b24),
                                  contract_8(c33, c34, c43, c44, b33, b34, b44, b44),
                                  contract_8(a11, a12, a21, a22, a13, a14, a23, a24),
                                  contract_8(a31, a32, a41, a42, a33, a34, a43, a44));

                U32 +=
                contract_3(a11, b11, c11) * contract_3(a12, b12, c12) * contract_3(a13, b13, c13) * contract_3(a14, b14, b14) *
                contract_3(a21, b21, c21) * contract_3(a22, b22, c22) * contract_3(a23, b23, c23) * contract_3(a24, b24, c24) *
                contract_3(a31, b31, c31) * contract_3(a32, b32, c32) * contract_3(a33, b33, c33) * contract_3(a34, b34, c34) *
                contract_3(a41, b41, c41) * contract_3(a42, b42, c42) * contract_3(a43, b43, c43) * contract_3(a44, b44, c44);
                
                U33 +=
                contract_3(a11, a21, b21) * contract_3(c21, c11, b11) * contract_3(a31, a41, b41) * contract_3(c41, c31, b31) *
                contract_3(a12, a22, b22) * contract_3(c22, c12, b12) * contract_3(a32, a42, b42) * contract_3(c42, c32, b32) *
                contract_3(a13, a23, b23) * contract_3(c23, c13, b13) * contract_3(a33, a43, b43) * contract_3(c43, c33, b33) *
                contract_3(a14, a24, b24) * contract_3(c24, c14, b14) * contract_3(a34, a44, b44) * contract_3(c44, c34, b34);
                U33 +=
                contract_3(a11, a21, b11) * contract_3(c21, c11, b21) * contract_3(a31, a41, b31) * contract_3(c41, c31, b41) *
                contract_3(a12, a22, b12) * contract_3(c22, c12, b22) * contract_3(a32, a42, b32) * contract_3(c42, c32, b42) *
                contract_3(a13, a23, b13) * contract_3(c23, c13, b23) * contract_3(a33, a43, b33) * contract_3(c43, c33, b43) *
                contract_3(a14, a24, b14) * contract_3(c24, c14, b24) * contract_3(a34, a44, b34) * contract_3(c44, c34, b44);
            }
        }
    }
    U[0] += U30 / (lattice->D - 1) / (lattice->D - 2) / (lattice->D - 3);
    U[1] += U31 / (lattice->D - 1) / (lattice->D - 2) / (lattice->D - 3) / 2;
    U[2] += U32 / (lattice->D - 1) / (lattice->D - 2) / (lattice->D - 3);
    U[3] += U33 / (lattice->D - 1) / (lattice->D - 2) / (lattice->D - 3) / 2;
}

void su3_site::quark_5(double* U) {
    double U40 = 0.0;
    double U41 = 0.0;
    double U42 = 0.0;
    double U43 = 0.0;
    for (int d1 = 1; d1 < lattice->D; d1++) {
        for (int d2 = 1; d2 < lattice->D; d2++) {
            if (d1 == d2) continue;
            for (int d3 = 1; d3 < lattice->D; d3++) {
                if (d2 == d3) continue;
                if (d1 == d3) continue;
                
                su3_site* A11 = this; su3_link a11 = Q(A11);
                su3_site* A12 = A11->forward[d2]; su3_link a12 = Q(A12);
                su3_site* A13 = A12->forward[d2]; su3_link a13 = Q(A13);
                su3_site* A14 = A13->forward[d2]; su3_link a14 = Q(A14);
                su3_site* A21 = A11->forward[d3]; su3_link a21 = Q(A21);
                su3_site* A22 = A21->forward[d2]; su3_link a22 = Q(A22);
                su3_site* A23 = A22->forward[d2]; su3_link a23 = Q(A23);
                su3_site* A24 = A23->forward[d2]; su3_link a24 = Q(A24);
                su3_site* A31 = A21->forward[d3]; su3_link a31 = Q(A31);
                su3_site* A32 = A31->forward[d2]; su3_link a32 = Q(A32);
                su3_site* A33 = A32->forward[d2]; su3_link a33 = Q(A33);
                su3_site* A34 = A33->forward[d2]; su3_link a34 = Q(A34);
                su3_site* A41 = A31->forward[d3]; su3_link a41 = Q(A41);
                su3_site* A42 = A41->forward[d2]; su3_link a42 = Q(A42);
                su3_site* A43 = A42->forward[d2]; su3_link a43 = Q(A43);
                su3_site* A44 = A43->forward[d2]; su3_link a44 = Q(A44);
                
                su3_site* B11 = A11->forward[d1]; su3_link b11 = Q(B11);
                su3_site* B12 = B11->forward[d2]; su3_link b12 = Q(B12);
                su3_site* B13 = B12->forward[d2]; su3_link b13 = Q(B13);
                su3_site* B14 = B13->forward[d2]; su3_link b14 = Q(B14);
                su3_site* B21 = B11->forward[d3]; su3_link b21 = Q(B21);
                su3_site* B22 = B21->forward[d2]; su3_link b22 = Q(B22);
                su3_site* B23 = B22->forward[d2]; su3_link b23 = Q(B23);
                su3_site* B24 = B23->forward[d2]; su3_link b24 = Q(B24);
                su3_site* B31 = B21->forward[d3]; su3_link b31 = Q(B31);
                su3_site* B32 = B31->forward[d2]; su3_link b32 = Q(B32);
                su3_site* B33 = B32->forward[d2]; su3_link b33 = Q(B33);
                su3_site* B34 = B33->forward[d2]; su3_link b34 = Q(B34);
                su3_site* B41 = B31->forward[d3]; su3_link b41 = Q(B41);
                su3_site* B42 = B41->forward[d2]; su3_link b42 = Q(B42);
                su3_site* B43 = B42->forward[d2]; su3_link b43 = Q(B43);
                su3_site* B44 = B43->forward[d2]; su3_link b44 = Q(B44);
                
                su3_site* C11 = B11->forward[d1]; su3_link c11 = Q(C11);
                su3_site* C12 = C11->forward[d2]; su3_link c12 = Q(C12);
                su3_site* C13 = C12->forward[d2]; su3_link c13 = Q(C13);
                su3_site* C14 = C13->forward[d2]; su3_link c14 = Q(C14);
                su3_site* C21 = C11->forward[d3]; su3_link c21 = Q(C21);
                su3_site* C22 = C21->forward[d2]; su3_link c22 = Q(C22);
                su3_site* C23 = C22->forward[d2]; su3_link c23 = Q(C23);
                su3_site* C24 = C23->forward[d2]; su3_link c24 = Q(C24);
                su3_site* C31 = C21->forward[d3]; su3_link c31 = Q(C31);
                su3_site* C32 = C31->forward[d2]; su3_link c32 = Q(C32);
                su3_site* C33 = C32->forward[d2]; su3_link c33 = Q(C33);
                su3_site* C34 = C33->forward[d2]; su3_link c34 = Q(C34);
                su3_site* C41 = C31->forward[d3]; su3_link c41 = Q(C41);
                su3_site* C42 = C41->forward[d2]; su3_link c42 = Q(C42);
                su3_site* C43 = C42->forward[d2]; su3_link c43 = Q(C43);
                su3_site* C44 = C43->forward[d2]; su3_link c44 = Q(C44);
                
                su3_site* D11 = C11->forward[d1]; su3_link d11 = Q(D11);
                su3_site* D12 = D11->forward[d2]; su3_link d12 = Q(D12);
                su3_site* D13 = D12->forward[d2]; su3_link d13 = Q(D13);
                su3_site* D14 = D13->forward[d2]; su3_link d14 = Q(D14);
                su3_site* D21 = D11->forward[d3]; su3_link d21 = Q(D21);
                su3_site* D22 = D21->forward[d2]; su3_link d22 = Q(D22);
                su3_site* D23 = D22->forward[d2]; su3_link d23 = Q(D23);
                su3_site* D24 = D23->forward[d2]; su3_link d24 = Q(D24);
                su3_site* D31 = D21->forward[d3]; su3_link d31 = Q(D31);
                su3_site* D32 = D31->forward[d2]; su3_link d32 = Q(D32);
                su3_site* D33 = D32->forward[d2]; su3_link d33 = Q(D33);
                su3_site* D34 = D33->forward[d2]; su3_link d34 = Q(D34);
                su3_site* D41 = D31->forward[d3]; su3_link d41 = Q(D41);
                su3_site* D42 = D41->forward[d2]; su3_link d42 = Q(D42);
                su3_site* D43 = D42->forward[d2]; su3_link d43 = Q(D43);
                su3_site* D44 = D43->forward[d2]; su3_link d44 = Q(D44);

                su3_site* E11 = D11->forward[d1]; su3_link e11 = Q(E11);
                su3_site* E12 = E11->forward[d2]; su3_link e12 = Q(E12);
                su3_site* E13 = E12->forward[d2]; su3_link e13 = Q(E13);
                su3_site* E14 = E13->forward[d2]; su3_link e14 = Q(E14);
                su3_site* E21 = E11->forward[d3]; su3_link e21 = Q(E21);
                su3_site* E22 = E21->forward[d2]; su3_link e22 = Q(E22);
                su3_site* E23 = E22->forward[d2]; su3_link e23 = Q(E23);
                su3_site* E24 = E23->forward[d2]; su3_link e24 = Q(E24);
                su3_site* E31 = E21->forward[d3]; su3_link e31 = Q(E31);
                su3_site* E32 = E31->forward[d2]; su3_link e32 = Q(E32);
                su3_site* E33 = E32->forward[d2]; su3_link e33 = Q(E33);
                su3_site* E34 = E33->forward[d2]; su3_link e34 = Q(E34);
                su3_site* E41 = E31->forward[d3]; su3_link e41 = Q(E41);
                su3_site* E42 = E41->forward[d2]; su3_link e42 = Q(E42);
                su3_site* E43 = E42->forward[d2]; su3_link e43 = Q(E43);
                su3_site* E44 = E43->forward[d2]; su3_link e44 = Q(E44);

                su3_site* F11 = E11->forward[d1]; su3_link f11 = Q(F11);
                su3_site* F12 = F11->forward[d2]; su3_link f12 = Q(F12);
                su3_site* F13 = F12->forward[d2]; su3_link f13 = Q(F13);
                su3_site* F14 = F13->forward[d2]; su3_link f14 = Q(F14);
                su3_site* F21 = F11->forward[d3]; su3_link f21 = Q(F21);
                su3_site* F22 = F21->forward[d2]; su3_link f22 = Q(F22);
                su3_site* F23 = F22->forward[d2]; su3_link f23 = Q(F23);
                su3_site* F24 = F23->forward[d2]; su3_link f24 = Q(F24);
                su3_site* F31 = F21->forward[d3]; su3_link f31 = Q(F31);
                su3_site* F32 = F31->forward[d2]; su3_link f32 = Q(F32);
                su3_site* F33 = F32->forward[d2]; su3_link f33 = Q(F33);
                su3_site* F34 = F33->forward[d2]; su3_link f34 = Q(F34);
                su3_site* F41 = F31->forward[d3]; su3_link f41 = Q(F41);
                su3_site* F42 = F41->forward[d2]; su3_link f42 = Q(F42);
                su3_site* F43 = F42->forward[d2]; su3_link f43 = Q(F43);
                su3_site* F44 = F43->forward[d2]; su3_link f44 = Q(F44);
                
                U40 += contract_3(contract_4(contract_8(a11, a21, a12, a22, b11, b21, b12, b22),
                                             contract_8(a31, a41, a32, a42, b31, b41, b32, b42),
                                             contract_8(a13, a23, a14, a24, b13, b23, b14, b24),
                                             contract_8(a33, a43, a34, a44, b33, b43, b34, b44)),
                                  contract_4(contract_8(c11, c21, c12, c22, d11, d21, d12, d22),
                                             contract_8(c31, c41, c32, c42, d31, d41, d32, d42),
                                             contract_8(c13, c23, c14, c24, d13, d23, d14, d24),
                                             contract_8(c33, c43, c34, c44, d33, d43, d34, d44)),
                                  contract_4(contract_8(e11, e21, e12, e22, f11, f21, f12, f22),
                                             contract_8(e31, e41, e32, e42, f31, f41, f32, f42),
                                             contract_8(e13, e23, e14, e24, f13, f23, f14, f24),
                                             contract_8(e33, e43, e34, e44, f33, f43, f34, f44)));
                
                U41 += contract_3(contract_4(contract_8(a11, a21, a12, a22, b11, b21, b12, b22),
                                             contract_8(a31, a41, a32, a42, b31, b41, b32, b42),
                                             contract_8(c11, c21, c12, c22, d11, d21, d12, d22),
                                             contract_8(c31, c41, c32, c42, d31, d41, d32, d42)),
                                  contract_4(contract_8(a13, a23, a14, a24, b13, b23, b14, b24),
                                             contract_8(a33, a43, a34, a44, b33, b43, b34, b44),
                                             contract_8(c13, c23, c14, c24, d13, d23, d14, d24),
                                             contract_8(c33, c43, c34, c44, d33, d43, d34, d44)),
                                  contract_4(contract_8(e11, e21, e12, e22, f11, f21, f12, f22),
                                             contract_8(e31, e41, e32, e42, f31, f41, f32, f42),
                                             contract_8(e13, e23, e14, e24, f13, f23, f14, f24),
                                             contract_8(e33, e43, e34, e44, f33, f43, f34, f44)));
                U41 += contract_3(contract_4(contract_8(a11, a21, a12, a22, b11, b21, b12, b22),
                                             contract_8(a31, a41, a32, a42, b31, b41, b32, b42),
                                             contract_8(a13, a23, a14, a24, b13, b23, b14, b24),
                                             contract_8(a33, a43, a34, a44, b33, b43, b34, b44)),
                                  contract_4(contract_8(c11, c21, c12, c22, d11, d21, d12, d22),
                                             contract_8(c31, c41, c32, c42, d31, d41, d32, d42),
                                             contract_8(e11, e21, e12, e22, f11, f21, f12, f22),
                                             contract_8(e31, e41, e32, e42, f31, f41, f32, f42)),
                                  contract_4(contract_8(c13, c23, c14, c24, d13, d23, d14, d24),
                                             contract_8(c33, c43, c34, c44, d33, d43, d34, d44),
                                             contract_8(e13, e23, e14, e24, f13, f23, f14, f24),
                                             contract_8(e33, e43, e34, e44, f33, f43, f34, f44)));
                
                U42 +=
                contract_3(a11, b11, c11) * contract_3(d11, e11, f11) *
                contract_3(a12, b12, c12) * contract_3(d12, e12, f12) *
                contract_3(a13, b13, c13) * contract_3(d13, e13, f13) *
                contract_3(a14, b14, c14) * contract_3(d14, e14, f14) *
                contract_3(a21, b21, c21) * contract_3(d21, e21, f21) *
                contract_3(a22, b22, c22) * contract_3(d22, e22, f22) *
                contract_3(a23, b23, c23) * contract_3(d23, e23, f23) *
                contract_3(a24, b24, c24) * contract_3(d24, e24, f24) *
                contract_3(a31, b31, c31) * contract_3(d31, e31, f21) *
                contract_3(a32, b32, c32) * contract_3(d32, e32, f32) *
                contract_3(a33, b33, c33) * contract_3(d33, e33, f33) *
                contract_3(a34, b34, c34) * contract_3(d34, e34, f34) *
                contract_3(a41, b41, c41) * contract_3(d41, e41, f41) *
                contract_3(a42, b42, c42) * contract_3(d42, e42, f42) *
                contract_3(a43, b43, c43) * contract_3(d43, e43, f43) *
                contract_3(a44, b44, c44) * contract_3(d44, e44, f44);
                
                U43 +=
                contract_3(a11, a12, b12) * contract_3(b11, c11, c12) *
                contract_3(a21, a22, b22) * contract_3(b21, c21, c22) *
                contract_3(a31, a32, b32) * contract_3(b31, c31, c32) *
                contract_3(a41, a42, b42) * contract_3(b41, c41, c42) *
                contract_3(a13, a14, b14) * contract_3(b13, c13, c14) *
                contract_3(a23, a24, b24) * contract_3(b23, c23, c24) *
                contract_3(a33, a34, b34) * contract_3(b33, c33, c34) *
                contract_3(a43, a44, b44) * contract_3(b43, c43, c44) *
                contract_3(d11, d12, e12) * contract_3(e11, f11, f12) *
                contract_3(d21, d22, e22) * contract_3(e21, f21, f22) *
                contract_3(d31, d32, e32) * contract_3(e31, f31, f32) *
                contract_3(d41, d42, e42) * contract_3(e41, f41, f42) *
                contract_3(d13, d14, e14) * contract_3(e13, f13, f14) *
                contract_3(d23, d24, e24) * contract_3(e23, f23, f24) *
                contract_3(d33, d34, e34) * contract_3(e33, f33, f34) *
                contract_3(d43, d44, e44) * contract_3(e43, f43, f44);
                U43 +=
                contract_3(a11, a12, b11) * contract_3(b12, c11, c12) *
                contract_3(a21, a22, b21) * contract_3(b22, c21, c22) *
                contract_3(a31, a32, b31) * contract_3(b32, c31, c32) *
                contract_3(a41, a42, b41) * contract_3(b42, c41, c42) *
                contract_3(a13, a14, b13) * contract_3(b14, c13, c14) *
                contract_3(a23, a24, b23) * contract_3(b24, c23, c24) *
                contract_3(a33, a34, b33) * contract_3(b34, c33, c34) *
                contract_3(a43, a44, b43) * contract_3(b44, c43, c44) *
                contract_3(d11, d12, e11) * contract_3(e12, f11, f12) *
                contract_3(d21, d22, e21) * contract_3(e22, f21, f22) *
                contract_3(d31, d32, e31) * contract_3(e32, f31, f32) *
                contract_3(d41, d42, e41) * contract_3(e42, f41, f42) *
                contract_3(d13, d14, e13) * contract_3(e14, f13, f14) *
                contract_3(d23, d24, e23) * contract_3(e24, f23, f24) *
                contract_3(d33, d34, e33) * contract_3(e34, f33, f34) *
                contract_3(d43, d44, e43) * contract_3(e44, f43, f44);
            }
        }
    }
    U[0] += U40 / (lattice->D - 1) / (lattice->D - 2) / (lattice->D - 3);
    U[1] += U41 / (lattice->D - 1) / (lattice->D - 2) / (lattice->D - 3) / 2;
    U[2] += U42 / (lattice->D - 1) / (lattice->D - 2) / (lattice->D - 3);
    U[3] += U43 / (lattice->D - 1) / (lattice->D - 2) / (lattice->D - 3) / 2;
}

//  :     :     :
//  :     :     :
// A4----B4----C4
//  |     |     |
//  |     |     |
// A3----B3----C3
//  |     |     |
//  |     |     |
// A2----B2----C2   ^ d2
//  |     |     |   |
//  |     |     |   |
// A1----B1----C1   +-----> d1

void su3_site::quark_special(double* U) {
    double U1 = 0.0;
    double U2 = 0.0;
    double U3 = 0.0;
    double U4 = 0.0;
    double U5 = 0.0;
    double U6 = 0.0;
    double U7 = 0.0;
    su3_site* A[8];
    su3_site* B[8];
    su3_site* C[8];
    su3_link a[8];
    su3_link b[8];
    su3_link c[8];
    for (int d1 = 1; d1 < lattice->D; d1++) {
        for (int d2 = 1; d2 < lattice->D; d2++) {
            if (d1 == d2) continue;
            
            A[0] = this; a[0] = A[0]->link[0];
            B[0] = forward[d1]; b[0] = B[0]->link[0];
            C[0] = forward[d1]->forward[d1]; c[0] = C[0]->link[0];
            for (int i = 1; i < 8; i++) {
                A[i] = A[i-1]->forward[d2]; a[i] = A[i]->link[0];
                B[i] = B[i-1]->forward[d2]; b[i] = B[i]->link[0];
                C[i] = C[i-1]->forward[d2]; c[i] = C[i]->link[0];
            }
            
            U1 += contract_3(contract_8(a[0], a[7], a[1], a[6], a[2], a[5], a[3], a[4]),
                             contract_8(b[0], b[7], b[1], b[6], b[2], b[5], b[3], b[4]),
                             contract_8(c[0], c[7], c[1], c[6], c[2], c[5], c[3], c[4]));
            U2 += contract_3(contract_8(a[0], a[4], a[1], a[5], a[2], a[6], a[3], a[7]),
                             contract_8(b[0], b[4], b[1], b[5], b[2], b[6], b[3], b[7]),
                             contract_8(c[0], c[4], c[1], c[5], c[2], c[6], c[3], c[7]));
            U3 += contract_3(a[0], b[7], c[7]) * contract_3(a[7], b[0], c[0]) *
                  contract_3(a[1], b[1], c[6]) * contract_3(a[6], b[6], c[1]) *
                  contract_3(a[2], b[5], c[5]) * contract_3(a[5], b[2], c[2]) *
                  contract_3(a[3], b[3], c[4]) * contract_3(a[4], b[4], c[3]);
            U4 += contract_3(a[0], b[4], c[4]) * contract_3(a[4], b[0], c[0]) *
                  contract_3(a[1], b[1], c[5]) * contract_3(a[5], b[5], c[1]) *
                  contract_3(a[2], b[6], c[6]) * contract_3(a[6], b[2], c[2]) *
                  contract_3(a[3], b[3], c[7]) * contract_3(a[7], b[7], c[3]);
            U5 += contract_3(a[0], a[1], b[0]) * contract_3(a[3], b[3], b[2]);
            U6 += contract_3(a[0], a[1], b[0]) * contract_3(a[4], b[4], b[3]);
            U7 += contract_3(a[0], a[1], b[0]) * contract_3(a[5], b[5], b[4]);
        }
    }
    U[0] += U1 / (lattice->D - 1) / (lattice->D - 2);
    U[1] += U2 / (lattice->D - 1) / (lattice->D - 2);
    U[2] += U3 / (lattice->D - 1) / (lattice->D - 2);
    U[3] += U4 / (lattice->D - 1) / (lattice->D - 2);
    U[4] += U5 / (lattice->D - 1) / (lattice->D - 2);
    U[5] += U6 / (lattice->D - 1) / (lattice->D - 2);
    U[6] += U7 / (lattice->D - 1) / (lattice->D - 2);
}

double su3_site::contract_distance(int x1, int y1, int x2, int y2) {
    double U = 0.0;
    for (int d1 = 1; d1 < lattice->D; d1++) {
        for (int d2 = 1; d2 < lattice->D; d2++) {
            if (d1 == d2) continue;
            su3_site* A = this;
            su3_site* B = this;
            su3_site* C = this;
            for (int x = 0; x < x1; x++) B = B->forward[d1];
            for (int y = 0; y < y1; y++) B = B->forward[d2];
            for (int x = 0; x < x2; x++) C = C->forward[d1];
            for (int y = 0; y < y2; y++) C = C->forward[d2];
            
            U += contract_3(Q(A), Q(B), Q(C));
        }
    }
    return U / (lattice->D - 1) / (lattice->D - 2);
}

void su3_site::quark_distance(double* U) {
    U[0] += contract_distance(1, 0, 0, 1); // L = 1.381
    U[1] += contract_distance(1, 1, 0, 2); // L = 1.609
    U[2] += contract_distance(0, 2, 2, 1); // L = 2.157
    U[3] += contract_distance(1, 3, 3, 1); // L = 3.051
    U[4] += contract_distance(0, 4, 3, 2); // L = 3.737
    U[5] += contract_distance(1, 4, 4, 1); // L = 4.163
    U[6] += contract_distance(0, 5, 4, 3); // L = 4.824
    U[7] += contract_distance(1, 5, 5, 2); // L = 5.161
    U[8] += contract_distance(3, 5, 6, 0); // L = 5.887
    U[9] += contract_distance(2, 6, 6, 1); // L = 6.270
    U[10] += contract_distance(0, 7, 6, 4); // L = 6.973
    U[11] += contract_distance(7, 2, 2, 7); // L = 7.210
    U[12] += contract_distance(3, 7, 8, 1); // L = 7.829
    U[13] += contract_distance(0, 8, 7, 4); // L = 8.042
    U[14] += contract_distance(2, 8, 8, 2); // L = 8.326
}

void su3_site::contract_hadron(su3_link* A, int size, int n) {
    
    // n had better be 3 x 2^n !!
    if (n == 3) return;
    
    int index[size];
    for (int i = 0; i < size; i++) index[i] = i;
    
    // randomly shuffle the links
    for (int i = 0; i < size; i++) {
        int r = rand() % size;
        int temp = index[r];
        index[r] = index[i];
        index[i] = temp;
    }

    // contract random pairs of links
    int n2 = n / 2;
    su3_link a[n2];
    for (int i = 0; i < n2; i++) {
        int i2 = i * 2;
        int r1 = index[i2];
        int r2 = index[i2 + 1];
        a[i] = contract_2(A[r1], A[r2]);
    }
    
    // copy the contracted values back to the array
    for (int i = 0; i < n2; i++) A[i] = a[i];
    
    // recursively contract the array down smaller
    if (n2 != 3) contract_hadron(A, n2, n2);
}

double su3_site::contract_baryon(su3_link* A, int size, int n) {
    int index[size];
    for (int i = 0; i < size; i++) index[i] = i;
    
    // randomly shuffle the array indices
    for (int i = 0; i < size; i++) {
        int r = rand() % n;
        int temp = index[r];
        index[r] = index[i];
        index[i] = temp;
    }
    
    int n3 = n / 3;
    double a[n3];
    for (int i = 0; i < n3; i++) {
        int i3 = i * 3;
        int r1 = index[i3];
        int r2 = index[i3 + 1];
        int r3 = index[i3 + 2];
        a[i] = contract_3(A[r1], A[r2], A[r3]);
    }
    
    double U = 1.0;
    for (int i = 0; i < n3; i++) U *= a[i];
    return U;
}

double su3_site::separate_baryons(su3_link* A, int size, int n) {
    int index[size];
    for (int i = 0; i < size; i++) index[i] = i;
    
    // randomly shuffle the array indices
    for (int i = 0; i < size; i++) {
        int r = rand() % n;
        int temp = index[r];
        index[r] = index[i];
        index[i] = temp;
    }
    
    int n3 = n / 3;
    double U = 0.0;
    for (int i = 0; i < n3; i++) {
        int i3 = i * 3;
        int r1 = index[i3];
        int r2 = index[i3 + 1];
        int r3 = index[i3 + 2];
        U += contract_3(A[r1], A[r2], A[r3]);
    }
    
    return U / n3;
}

double su3_site::contract_rand(int D1, int D2, int D3, int n, int type) {
    int size = D1 * D2 * D3;
    su3_link A[size];
    
    // pick random directions for x, y, and z
    int d1 = rand() % 3;
    int d2 = (d1 + 1) % 3;
    int d3 = (d2 + 1) % 3;
    d1++; d2++; d3++;
    
    int s = 0;
    su3_site* X, *Y, *Z;
    X = this;
    for (int x = 0; x < D1; x++) {
        Y = X;
        for (int y = 0; y < D2; y++) {
            Z = Y;
            for (int z = 0; z < D3; z++) {
                A[s++] = Q(Z);
                Z = Z->forward[d3];
            }
            Y = Y->forward[d2];
        }
        X = X->forward[d1];
    }
    
    if (type == 0) {
        contract_hadron(A, size, n);
        return contract_3(A[0], A[1], A[2]);
    } else if (type == 1){
        return contract_baryon(A, size, n);
    } else {
        return separate_baryons(A, size, n);
    }
}

void su3_site::quark_rand(double* U, int type) {

    // V = 3
    U[ 0] += contract_rand(3, 1, 1, 3, type);

    // V = 6
    U[ 1] += contract_rand(3, 2, 1, 6, type);
    U[ 2] += contract_rand(6, 1, 1, 6, type);

    // V = 12
    U[ 3] += contract_rand(3, 4, 1, 12, type);
    U[ 4] += contract_rand(3, 2, 2, 12, type);
    U[ 5] += contract_rand(6, 2, 1, 12, type);
    U[ 6] += contract_rand(12, 1, 1, 12, type);

    // V = 24
    U[ 7] += contract_rand(3, 8, 1, 24, type);
    U[ 8] += contract_rand(3, 4, 2, 24, type);
    U[ 9] += contract_rand(6, 4, 1, 24, type);
    U[10] += contract_rand(6, 2, 2, 24, type);
    U[11] += contract_rand(12, 2, 1, 24, type);

    // V = 48
    U[12] += contract_rand(3, 8, 2, 48, type);
    U[13] += contract_rand(3, 4, 4, 48, type);
    U[14] += contract_rand(6, 8, 1, 48, type);
    U[15] += contract_rand(6, 4, 2, 48, type);
    U[16] += contract_rand(12, 4, 1, 48, type);
    U[17] += contract_rand(12, 2, 2, 48, type);

    // V = 96
    U[18] += contract_rand(3, 8, 4, 96, type);
    U[19] += contract_rand(6, 8, 2, 96, type);
    U[20] += contract_rand(6, 4, 4, 96, type);
    U[21] += contract_rand(12, 8, 1, 96, type);
    U[22] += contract_rand(12, 4, 2, 96, type);

    // V = 192
    U[23] += contract_rand(3, 8, 8, 192, type);
    U[24] += contract_rand(6, 8, 4, 192, type);
    U[25] += contract_rand(12, 8, 2, 192, type);
    U[26] += contract_rand(12, 4, 4, 192, type);
    
    // V = 384
    U[27] += contract_rand(6, 8, 8, 384, type);
    U[28] += contract_rand(12, 8, 4, 384, type);

    // V = 768
    U[29] += contract_rand(12, 8, 8, 768, type);
}

void su3_site::quark_density(double* U, int type) {
    U[ 0] += contract_rand(4, 4, 4, 48, type);
    U[ 1] += contract_rand(4, 4, 5, 48, type);
    U[ 2] += contract_rand(4, 5, 5, 48, type);
    U[ 3] += contract_rand(5, 5, 5, 48, type);
    U[ 4] += contract_rand(5, 5, 6, 48, type);
    U[ 5] += contract_rand(5, 6, 6, 48, type);
    U[ 6] += contract_rand(6, 6, 6, 48, type);
    U[ 7] += contract_rand(6, 6, 7, 48, type);
    U[ 8] += contract_rand(6, 7, 7, 48, type);
    U[ 9] += contract_rand(7, 7, 7, 48, type);
    U[10] += contract_rand(7, 7, 8, 48, type);
    U[11] += contract_rand(7, 8, 8, 48, type);
    U[12] += contract_rand(8, 8, 8, 48, type);
    U[13] += contract_rand(8, 8, 9, 48, type);
    U[14] += contract_rand(8, 9, 9, 48, type);
    U[15] += contract_rand(9, 9, 9, 48, type);
    U[16] += contract_rand(9, 9, 10, 48, type);
    U[17] += contract_rand(9, 10, 10, 48, type);
    U[18] += contract_rand(10, 10, 10, 48, type);
    U[19] += contract_rand(10, 10, 11, 48, type);
    U[20] += contract_rand(10, 11, 11, 48, type);
    U[21] += contract_rand(11, 11, 11, 48, type);
    U[22] += contract_rand(11, 11, 12, 48, type);
    U[23] += contract_rand(11, 12, 12, 48, type);
    U[24] += contract_rand(12, 12, 12, 48, type);
}

// one-dimensional quark gas

void su3_site::quark_gas_1(double* U, bool baryon) {
    
    for (int R = 0; R < lattice->N / 2; R++) {
        double u = 0.0;
        for (int d1 = 1; d1 < lattice->D; d1++) {
            for (int d2 = 1; d2 < lattice->D; d2++) {
                if (d1 == d2) continue;
                
                su3_site* A1 = this;
                su3_site* A2 = forward[d1];
                su3_site* A3 = forward[d2];

                su3_site* B1 = A3->forward[d1]->forward[d1];
                for (int i = 0; i < R; i++) B1 = B1->forward[d1];
                su3_site* B2 = B1->backward[d1];
                su3_site* B3 = B1->backward[d2];
                
                if (baryon) {
                    u += contract_3(Q(A1), Q(A2), Q(A3)) * contract_3(Q(B1), Q(B2), Q(B3));
                } else {
                    u += contract_6(Q(A1), Q(B1), Q(A2), Q(B2), Q(A3), Q(B3));
                }
            }
        }
        U[R] += u / (lattice->D - 1) / (lattice->D - 2);
    }
}

void su3_site::quark_gas_2(double* U, bool baryon) {
    
    for (int R = 0; R < lattice->N / 2; R++) {
        double u = 0.0;
        for (int d1 = 1; d1 < lattice->D; d1++) {
            for (int d2 = 1; d2 < lattice->D; d2++) {
                if (d1 == d2) continue;
                
                su3_site* A1 = this;
                su3_site* A2 = A1->forward[d1];
                su3_site* A3 = A1->forward[d2];
                
                su3_site* B1 = A3->forward[d1]->forward[d1];
                for (int i = 0; i < R; i++) B1 = B1->forward[d1];
                su3_site* B2 = B1->backward[d1];
                su3_site* B3 = B1->backward[d2];
                
                su3_site* C1 = A3->forward[d2];
                for (int i = 0; i < R; i++) C1 = C1->forward[d2];
                su3_site* C2 = C1->forward[d1];
                su3_site* C3 = C1->forward[d2];
                
                su3_site* D1 = C3->forward[d1]->forward[d1];
                for (int i = 0; i < R; i++) D1 = D1->forward[d1];
                su3_site* D2 = D1->backward[d1];
                su3_site* D3 = D1->backward[d2];

                if (baryon) {
                    u += contract_3(Q(A1), Q(A2), Q(A3)) * contract_3(Q(B1), Q(B2), Q(B3)) * contract_3(Q(C1), Q(C2), Q(C3)) * contract_3(Q(D1), Q(D2), Q(D3));
                } else {
                    u += contract_3(contract_4(Q(A1), Q(B1), Q(C1), Q(D1)),
                                    contract_4(Q(A2), Q(B2), Q(C2), Q(D2)),
                                    contract_4(Q(A3), Q(B3), Q(C3), Q(D3)));
                }
            }
        }
        U[R] += u / (lattice->D - 1) / (lattice->D - 2);
    }
}

void su3_site::quark_gas_3(double* U, bool baryon) {
    
    for (int R = 0; R < lattice->N / 2; R++) {
        double u = 0.0;
        for (int d1 = 1; d1 < lattice->D; d1++) {
            for (int d2 = 1; d2 < lattice->D; d2++) {
                if (d1 == d2) continue;
                for (int d3 = 1; d3 < lattice->D; d3++) {
                    if (d2 == d3) continue;
                    if (d1 == d3) continue;

                    su3_site* A11 = this;
                    su3_site* A12 = A11->forward[d1];
                    su3_site* A13 = A11->forward[d2];
                    
                    su3_site* B11 = A13->forward[d1]->forward[d1];
                    for (int i = 0; i < R; i++) B11 = B11->forward[d1];
                    su3_site* B12 = B11->backward[d1];
                    su3_site* B13 = B11->backward[d2];
                    
                    su3_site* C11 = A13->forward[d2];
                    for (int i = 0; i < R; i++) C11 = C11->forward[d2];
                    su3_site* C12 = C11->forward[d1];
                    su3_site* C13 = C11->forward[d2];
                    
                    su3_site* D11 = C13->forward[d1]->forward[d1];
                    for (int i = 0; i < R; i++) D11 = D11->forward[d1];
                    su3_site* D12 = D11->backward[d1];
                    su3_site* D13 = D11->backward[d2];

                    su3_site* A21 = A11->forward[d3];
                    su3_site* A22 = A21->forward[d1];
                    su3_site* A23 = A21->forward[d2];
                    
                    su3_site* B21 = A23->forward[d1]->forward[d1];
                    for (int i = 0; i < R; i++) B21 = B21->forward[d1];
                    su3_site* B22 = B21->backward[d1];
                    su3_site* B23 = B21->backward[d2];
                    
                    su3_site* C21 = A23->forward[d2];
                    for (int i = 0; i < R; i++) C21 = C21->forward[d2];
                    su3_site* C22 = C21->forward[d1];
                    su3_site* C23 = C21->forward[d2];
                    
                    su3_site* D21 = C23->forward[d1]->forward[d1];
                    for (int i = 0; i < R; i++) D21 = D21->forward[d1];
                    su3_site* D22 = D21->backward[d1];
                    su3_site* D23 = D21->backward[d2];

                    if (baryon) {
                        u += contract_3(Q(A11), Q(A12), Q(A13)) * contract_3(Q(B11), Q(B12), Q(B13)) * contract_3(Q(C11), Q(C12), Q(C13)) * contract_3(Q(D11), Q(D12), Q(D13)) * contract_3(Q(A21), Q(A22), Q(A23)) * contract_3(Q(B21), Q(B22), Q(B23)) * contract_3(Q(C21), Q(C22), Q(C23)) * contract_3(Q(D21), Q(D22), Q(D23));
                    } else {
                        u += contract_6(contract_4(Q(A11), Q(B11), Q(C11), Q(D11)),
                                        contract_4(Q(A12), Q(B12), Q(C12), Q(D12)),
                                        contract_4(Q(A13), Q(B13), Q(C13), Q(D13)),
                                        contract_4(Q(A21), Q(B21), Q(C21), Q(D21)),
                                        contract_4(Q(A22), Q(B22), Q(C22), Q(D22)),
                                        contract_4(Q(A23), Q(B23), Q(C23), Q(D23)));
                    }
                }
            }
        }
        U[R] += u / (lattice->D - 1) / (lattice->D - 2) / (lattice->D - 3);
    }
}

su3_link su3_site::overrelax(su3_link g) {
    
    // this function doesn't work in SU(3)
    
    double w = 1.7;
    int n_max = 5;
    
    su3_link g1 = g - su3_identity; // original relaxation matrix minus identity
    su3_link gn = su3_identity; // original g raised to the nth power
    su3_link gw = su3_zero; // overrelaxation matrix
    
    for (int n = 0; n <= n_max; n++) {
        double wn = tgamma(w + 1) / tgamma(w + 1 - n) / tgamma(n + 1);
        gw += wn * gn;
        gn *= g1;
    }
    return make_unitary(gw);
}

su3_link su3_site::make_unitary(su3_link g) {
    
    // make sure the determinant is 1
    complex<double> det = g.determinant();
//    if (abs(norm(det) - 1.0) > 1e-10) {
        double n = sqrt(norm(det));
        double a = arg(det);
        g /= cbrt(n);
        g *= polar(1.0, -a / 3.0);
//    }
    
    // create 3-vectors from two random rows
    // the third row is unused
    int r1 = rand() % 3;
    int r2 = (r1 + 1) % 3;
    int r3 = (r1 + 2) % 3;
    su3_vector u = g.row(r1);
    su3_vector v = g.row(r2);
    
    // normalize u
    u /= u.norm();
    
    // make v orthonormal to u
    v = v - u * u.dot(v);
    v /= v.norm();
    
    // create a third orthonormal vector w
    su3_vector w = u.cross(v);
    
    g.row(r1) = u;
    g.row(r2) = v;
    g.row(r3) = w;
    
    return g;
}

void su3_site::relax(bool coulomb) {
    while (!lock()) {}
    su3_link g = make_unitary(su3_identity - sum_G(coulomb) * I / 4);
//    g = overrelax(g);
    
    // apply the gauge transformation (including the time direction)
    for (int d = 0; d < lattice->D; d++) {
        set_link(d, g * link[d]);
        backward[d]->set_link(d, backward[d]->link[d] * g.transpose().conjugate());
    }
    unlock();
}

double su3_site::sum_landau() {
    double U = 0.0;
    for (int d = 0; d < lattice->D; d++) U += link[d].trace().real();
    return U / 3 / lattice->D;
}

double su3_site::sum_coulomb() {
    double U = 0.0;
    for (int d = 1; d < lattice->D; d++) U += link[d].trace().real();
    return U / 3 / (lattice->D - 1);
}

su3_link su3_site::sum_G(bool coulomb) {
    su3_link G = su3_zero;
    for (int d = coulomb ? 1 : 0; d < lattice->D; d++) {
        G += (link[d] - link_inverse[d]) / 2 / I;
        G -= (backward[d]->link[d] - backward[d]->link_inverse[d]) / 2 / I;
    }
    return G;
}

double su3_site::sweep() {
    double accept = 0;
    lock();
    for (int d = 0; d < lattice->D; d++) accept += sweep_link(d);
    unlock();
    return accept / lattice->D;
}

double su3_site::sweep_link(int d1) {
    su3_link old_link = link[d1];
    su3_link new_link = old_link * create_link(false);
    
    su3_link u;
    su3_link A = su3_zero;
    for (int d2 = 0; d2 < lattice->D; d2++) {
        if (d1 == d2) continue;
        
        // forward staple
        u = forward[d1]->link[d2];
        u *= forward[d2]->link_inverse[d1];
        u *= link_inverse[d2];
        A += u;

        // backward staple
        u = backward[d2]->forward[d1]->link_inverse[d2];
        u *= backward[d2]->link_inverse[d1];
        u *= backward[d2]->link[d2];
        A += u;
    }
    
    double dE = ((new_link - old_link) * A).trace().real() / 3;
    
    // metropolis algorithm
    if (dE >= 0 || lattice->rand() < exp(lattice->beta * dE)) {
        set_link(d1, new_link);
        return 1;
    }
    return 0;
}

su3_lattice::su3_lattice(int N, int T, int D, double beta, bool cold) {
    this->N = N;
    this->T = T;
    this->D = D;
    this->beta = beta;
    this->parallel = false;
    this->z = 50;
    this->sigma[0] = su2_identity;
    this->sigma[1] << 0, 1, 1, 0;
    this->sigma[2] << 0, -I, I, 0;
    this->sigma[3] << 1, 0, 0, -1;

    random_device rd;  // seed for the random number engine
    gen = mt19937(rd()); // mersenne_twister_engine seeded with rd()
    
    n_sites = pow(N, D - 1) * T;
    n_slice = n_sites / T; // number of sites in each time slice
    site = new su3_site[n_sites];
    
    for (int s = 0; s < n_sites; s++) {
        site[s].lattice = this;
        site[s].is_locked = NULL;
        site[s].reset_links(cold);
        
        for (int d = 0; d < D; d++) {
            update_i(s);
            move(d, +1);
            site[s].forward[d] = &site[update_s()];
            move(d, -2);
            site[s].backward[d] = &site[update_s()];
        }
    }
}

su3_lattice::su3_lattice(su3_lattice* lattice) {
    this->N = lattice->N;
    this->T = lattice->T;
    this->D = lattice->D;
    this->beta = lattice->beta;
    this->parallel = lattice->parallel;
    this->z = 50;
    this->sigma[0] = su2_identity;
    this->sigma[1] << 0, 1, 1, 0;
    this->sigma[2] << 0, -I, I, 0;
    this->sigma[3] << 1, 0, 0, -1;

    random_device rd;  // seed for the random number engine
    gen = mt19937(rd()); // mersenne_twister_engine seeded with rd()
    
    n_sites = pow(N, D - 1) * T;
    n_slice = n_sites / T; // number of sites in each time slice
    site = new su3_site[n_sites];
    
    for (int s = 0; s < n_sites; s++) {
        site[s].lattice = this;
        site[s].is_locked = NULL;
        site[s].copy_links(&lattice->site[s]);
        
        for (int d = 0; d < D; d++) {
            update_i(s);
            move(d, +1);
            site[s].forward[d] = &site[update_s()];
            move(d, -2);
            site[s].backward[d] = &site[update_s()];
        }
    }
}

su3_lattice::~su3_lattice() {
    delete[] site;
}

double su3_lattice::rand(double min, double max) {
    return ((double)gen()) / 4294967295.0 * (max - min) + min;
}

void su3_lattice::update_i(int s) {
    i[0] = s / (int)pow(N, D - 1);
    for (int d = 1; d < D; d++) i[d] = s / (int)pow(N, D - d - 1) % N;
}

int su3_lattice::update_s() {
    int s = 0;
    for (int d = 0; d < D; d++) s += i[d] * (int)pow(N, D - d - 1);
    return s;
}

void su3_lattice::move(int d, int n) {
    
    if (d == 0) {
        if (n > 0) {
            i[0] = (i[0] + n) % T;
        } else {
            i[0] = (i[0] + n + T) % T;
        }
    } else {
        if (n > 0) {
            i[d] = (i[d] + n) % N;
        } else {
            i[d] = (i[d] + n + N) % N;
        }
    }
}

double su3_lattice::plaq() {
    // compute the average plaquette of the lattice
    double U = 0;
    for (int s = 0; s < n_sites; s++) U += site[s].plaq();
    return U / n_sites;
}

double su3_lattice::wilson_loop(int a, int b) {
    // compute the average a x b wilson loop of the lattice
    double sum = 0;
    for (int s = 0; s < n_sites; s++) sum += site[s].wilson_loop(a, b);
    return sum / n_sites;
}

double su3_lattice::polyakov_loop(int r) {
    // compute the average product of polyakov loops spaced r apart
    double sum = 0;
    for (int s = 0; s < n_sites / T; s++) sum += site[s].polyakov_loop(r);
    return sum / n_sites * T;
}

double su3_lattice::correlator(int T) {
    // compute the average correlator of length T
    double sum = 0;
    for (int s = 0; s < n_sites; s++) sum += site[s].correlator(T);
    return sum / n_sites;
}

double su3_lattice::four_point(int T, int R) {
    // compute the average T x R 4-point function of the lattice
    double sum = 0;
    for (int s = 0; s < n_sites; s++) sum += site[s].four_point(T, R);
    return sum / n_sites;
}

void async_quark_distance(su3_site* site, int n, double* U) {
    for (int s = 0; s < n; s++) site[s].quark_distance(U);
}

void su3_lattice::quark_distance(double* U) {
    for (int q = 0; q < 15; q++) U[q] = 0.0;
    
    if (parallel) {
        future<void> quark_slice[T];
        double U_slice[T][15];
        for (int t = 0; t < T; t++) {
            for (int q = 0; q < 15; q++) U_slice[t][q] = 0.0;
            quark_slice[t] = async(async_quark_distance, &site[t * n_slice], n_slice, U_slice[t]);
        }
        
        for (int t = 0; t < T; t++) {
            quark_slice[t].get();
            for (int q = 0; q < 15; q++) {
                U[q] += U_slice[t][q];
            }
        }
    } else {
        for (int s = 0; s < n_sites; s++) {
            site[s].quark_distance(U);
        }
    }
    for (int q = 0; q < 15; q++) U[q] /= n_sites;
}

void async_quark_rand(su3_site* site, int n, double* U) {
    for (int s = 0; s < n; s++) {
        site[s].quark_rand(U, true);
        site[s].quark_rand(U+30, false);
    }
}

void su3_lattice::quark_rand(double* U) {
    for (int q = 0; q < 60; q++) U[q] = 0.0;
    
    if (parallel) {
        future<void> quark_slice[T];
        double U_slice[T][60];
        for (int t = 0; t < T; t++) {
            for (int q = 0; q < 60; q++) U_slice[t][q] = 0.0;
            quark_slice[t] = async(async_quark_rand, &site[t * n_slice], n_slice, U_slice[t]);
        }
        
        for (int t = 0; t < T; t++) {
            quark_slice[t].get();
            for (int q = 0; q < 60; q++) {
                U[q] += U_slice[t][q];
            }
        }
    } else {
        for (int s = 0; s < n_sites; s++) {
            site[s].quark_rand(U, true);
            site[s].quark_rand(U+30, false);
        }
    }
    for (int q = 0; q < 60; q++) U[q] /= n_sites;
}

void async_quark(su3_site* site, int n, double* U) {
    for (int s = 0; s < n; s++) {
        
        U[0] += site[s].mag_U();
        U[1] += site[s].abs_U();
//        site[s].quark_rand(U, 0);
//        site[s].quark_rand(U + 30, 1);
//        site[s].quark_density(U+30, 0);
        
//        site[s].quark_0(U);
//        site[s].quark_1(U+1);
//        site[s].quark_2(U+5);
//        site[s].quark_3(U+13);
//        site[s].quark_4(U+25);
//        site[s].quark_5(U+29);
//        site[s].quark_special(U+33);
//        site[s].quark_rand(U+40, 1);
//        site[s].quark_rand(U+70, 2);
//        site[s].quark_density(U+100, 1);
//        site[s].quark_density(U+125, 2);
//        site[s].quark_gas_1(U+150, true);
//        site[s].quark_gas_2(U+166, true);
//        site[s].quark_gas_3(U+182, true);
//        site[s].quark_gas_1(U+198, false);
//        site[s].quark_gas_2(U+214, false);
//        site[s].quark_gas_3(U+230, false);
    }
}

void su3_lattice::quark(double* U) {
    for (int q = 0; q < QUARK_MAX; q++) U[q] = 0.0;
    
    if (parallel) {
        future<void> quark_slice[T];
        double U_slice[T][QUARK_MAX];
        for (int t = 0; t < T; t++) {
            for (int q = 0; q < QUARK_MAX; q++) U_slice[t][q] = 0.0;
            quark_slice[t] = async(async_quark, &site[t * n_slice], n_slice, U_slice[t]);
        }
        
        for (int t = 0; t < T; t++) {
            quark_slice[t].get();
            for (int q = 0; q < QUARK_MAX; q++) {
                U[q] += U_slice[t][q];
            }
        }
    } else {
        for (int s = 0; s < n_sites; s++) {
            site[s].quark_0(U);
            site[s].quark_1(U+1);
            site[s].quark_2(U+5);
            site[s].quark_3(U+13);
            site[s].quark_4(U+25);
            site[s].quark_5(U+29);
            site[s].quark_special(U+33);
            site[s].quark_rand(U+40, 1);
            site[s].quark_rand(U+70, 2);
            site[s].quark_density(U+100, 1);
            site[s].quark_density(U+104, 2);
        }
    }
    for (int q = 0; q < QUARK_MAX; q++) U[q] /= n_sites;
}

double async_sweep(su3_site* site, int n) {
    double accept = 0.0;

    // start at a random spot in the time slice to mitigate parallel effects
    int start = rand() % n;
    int end = start + n;
    
    for (int s = start; s < end; s++) accept += site[s % n].sweep();

    return accept;
}

void su3_lattice::sweep(int n_sweeps) {
    
    if (n_sweeps) {
        for (int i = 0; i < n_sweeps; i++) sweep(0);
        return;
    }
    
    double accept = 0;
    
    if (parallel) {
        future<double> sweep_slice[T];
        for (int t = 0; t < T; t++) sweep_slice[t] = async(async_sweep, &site[t * n_slice], n_slice);
        for (int t = 0; t < T; t++) accept += sweep_slice[t].get();
    } else {
        for (int s = 0; s < n_sites; s++) accept += site[s].sweep();
    }
    accept /= n_sites;

    if (accept < 0.4) {
        z *= 1.2;
    } else if (accept > 0.6) {
        z *= 0.8;
    }
}

int su3_lattice::thermalize(int n_min, int n_max) {
    
    // create a cold start model
    su3_lattice coldStart = su3_lattice(N, T, D, beta, true);
    double hot_plaq = plaq();
    double cold_plaq = coldStart.plaq();
    int m;
    for (m = 0; hot_plaq < cold_plaq; m++) {
        // sweep both models until the average plaquette reaches the same value
        sweep();
        coldStart.sweep();
        hot_plaq = plaq();
        cold_plaq = coldStart.plaq();
//        cout << "dE = " << abs(cold_plaq - hot_plaq) << "\n";
        if (n_max && (m > n_max)) break;
    }
    
    // sweep a few more times
    int n_therm = max(n_min, m * 2);
    for (int i = m; i < n_therm; i++) sweep();
    
    // return value is total number of sweeps
    return n_therm;
}

double su3_lattice::ave_link_t() {
    // average trace of links in the t direction
    // this should be ~0 in Coulomb gauge
    
    double U = 0;
    
    for (int s = 0; s < n_sites; s++) U += site[s].link[0].trace().real() / 3;
    
    return U / n_sites;
}

double su3_lattice::theta(bool coulomb) {
    su3_link G1, G2;
    double th = 0.0;
    
    for (int s = 0; s < n_sites; s++) {
        G1 = site[s].sum_G(coulomb);
        G1 -= (G1.trace() * su3_identity) / 3;
        G2 = G1.transpose().conjugate();
        th += (G1 * G2).trace().real() / 3;
    }
    
    return th;
}

double async_relax(su3_site* site, int n, double error_target, bool coulomb) {
    double avg_link = 0.0;
    double th = 1.0;
    while (th > error_target) {
        
        // start at a random spot in the time slice to mitigate parallel effects
        int start = rand() % n;
        int end = start + n;
        
        for (int s = start; s < end; s++) site[s % n].relax(coulomb);
        
        // add up all the links
        su3_link G1, G2;
        th = 0.0;
        for (int s = 0; s < n; s++) {
            G1 = site[s].sum_G(coulomb);
            G1 -= (G1.trace() * su3_identity) / 3;
            G2 = G1.transpose().conjugate();
            th += (G1 * G2).trace().real() / 3;
        }
        avg_link = 0.0;
        for (int s = 0; s < n; s++) avg_link += site[s].sum_coulomb();
        avg_link /= n;
    }
    return avg_link;
}

long double su3_lattice::relax(long double error_target, bool coulomb) {
    
    double avg_link = 0.0;
    double th = 1.0;

    if (parallel) {
        future<double> avg_link_slice[T];
        for (int t = 0; t < T; t++)
            avg_link_slice[t] = async(async_relax, &site[t * n_slice], n_slice, error_target / T, coulomb);
        
        for (int t = 0; t < T; t++) avg_link += avg_link_slice[t].get();
        
        return avg_link / T;

    } else if (coulomb) {
        // for coulomb gauge, relax each time slice separately
        for (int t = 0; t < T; t++) {
            int s1 = t * n_slice;
            int s2 = s1 + n_slice;
            double avg_link_slice = 0.0;
            th = 1.0;
            while (th > (error_target / T)) {

                for (int s = s1; s < s2; s++) site[s].relax(coulomb);

                // add up all the links
                su3_link G1, G2;
                th = 0.0;
                for (int s = s1; s < s2; s++) {
                    G1 = site[s].sum_G(coulomb);
                    G1 -= (G1.trace() * su3_identity) / 3;
                    G2 = G1.transpose().conjugate();
                    th += (G1 * G2).trace().real() / 3;
                }
                avg_link_slice = 0.0;
                for (int s = s1; s < s2; s++) avg_link_slice += site[s].sum_coulomb();
                avg_link_slice /= n_slice;
//                cout << "avg_link = " << avg_link_slice << ", ";
//                cout << "theta = " << th << "\n";
            }
            avg_link += avg_link_slice;
        }
        return avg_link / T;
    } else {
        // for landau gauge, relax the entire lattice all at once
        while (th > error_target) {
            
            for (int s = 0; s < n_sites; s++) site[s].relax(coulomb);
            
            th = theta(coulomb);

            // add up all the links
            avg_link = 0.0;
            if (coulomb) {
                for (int s = 0; s < n_sites; s++) avg_link += site[s].sum_coulomb();
            } else {
                for (int s = 0; s < n_sites; s++) avg_link += site[s].sum_landau();
            }
            avg_link /= n_sites;
//            cout << "avg_link = " << avg_link << ", ";
//            cout << "theta = " << th << ", ";
//            cout << "plaq = " << plaq() << "\n";
        }
        return avg_link;
    }
}

void su3_lattice::write_four_point(const char* filename, bool coulomb) {
    
    ofstream file;
    file.open(filename, ofstream::app);
    
    int n_prop = N / 2;
    int t_prop = T / 2;
    file << setprecision(12);
    file << N << "^" << D - 1 << "*" << T << "\t";
    file << beta << "\t";
    file << (coulomb ? "C" : "L");
    for (int R = 1; R <= n_prop; R++) {
        for (int T = 1; T <= t_prop; T++) {
            file << "\t" << four_point(T, R);
        }
    }
    file << "\n";
    file.close();
}

void su3_lattice::write_correlator(const char* filename, bool coulomb) {
    
    ofstream file;
    file.open(filename, ofstream::app);
    
    int t_prop = T / 2;
    file << setprecision(12);
    file << N << "^" << D - 1 << "*" << T << "\t";
    file << beta << "\t";
    file << (coulomb ? "C*" : "L*");
    for (int T = 1; T <= t_prop; T++) {
        file << "\t" << correlator(T);
    }
    file << "\n";
    file.close();
}
