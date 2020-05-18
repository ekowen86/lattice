//
//  su3_p.cpp
//
//  Created by Evan Owen on 11/16/18.
//  Copyright © 2018 Evan Owen. All rights reserved.
//

#include <iostream>
#include <iomanip>
#include <fstream>
#include <random>
#include <cmath>
#include <future>
#include "Eigen/Dense"
#include "su3_p.hpp"

using namespace std;

void su3_site::lock() {
    
    // wait for this site and all adjacent links to be unlocked
    while (is_locked) {}
    is_locked = true;
}

void su3_site::unlock() {
    is_locked = false;
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
        // make sure the determinant is 1
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

complex<double> su3_site::contract_3(su3_link A, su3_link B, su3_link C) {
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
    return U / 6.0;
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
    return A.conjugate().transpose() / 2;
}

// *
// |
// *---*

double su3_site::contract_3x1_L() {
    complex<double> U = 0.0;
    for (int d1 = 1; d1 < lattice->D; d1++) {
        for (int d2 = d1 + 1; d2 < lattice->D; d2++) {
            su3_link A = link[0];
            su3_link B = forward[d1]->link[0];
            su3_link C = forward[d2]->link[0];
            U += contract_3(A, B, C);
        }
    }
    return U.real() / (lattice->D - 1) / (lattice->D - 2) * 2;
}

// *---*---*

double su3_site::contract_3x1() {
    complex<double> U = 0.0;
    for (int d = 1; d < lattice->D; d++) {
        su3_link A = link[0];
        su3_link B = forward[d]->link[0];
        su3_link C = forward[d]->forward[d]->link[0];
        U += contract_3(A, B, C);
    }
    return U.real() / (lattice->D - 1);
}

//   *
//  /|
// *---*
// |   |
// *---*

double su3_site::contract_3x2_L() {
    complex<double> U = 0.0;
    for (int d1 = 1; d1 < lattice->D; d1++) {
        for (int d2 = 1; d2 < lattice->D; d2++) {
            for (int d3 = 1; d3 < lattice->D; d3++) {
                if (d1 == d2) continue;
                if (d2 == d3) continue;
                if (d1 == d3) continue;
                
                su3_link A1 = link[0];
                su3_link A2 = forward[d1]->link[0];
                su3_link A = contract_2(A1, A2);
                su3_link B1 = forward[d2]->link[0];
                su3_link B2 = forward[d2]->forward[d1]->link[0];
                su3_link B = contract_2(B1, B2);
                su3_link C1 = forward[d3]->link[0];
                su3_link C2 = forward[d3]->forward[d1]->link[0];
                su3_link C = contract_2(C1, C2);
                U += contract_3(A, B, C);
            }
        }
    }
    return U.real() / (lattice->D - 1) / (lattice->D - 2) / (lattice->D - 3);
}

// *---*---*
// |   |   |
// *---*---*

double su3_site::contract_3x2() {
    complex<double> U = 0.0;
    for (int d1 = 1; d1 < lattice->D; d1++) {
        for (int d2 = 1; d2 < lattice->D; d2++) {
            if (d1 == d2) continue;
            su3_link A1 = link[0];
            su3_link A2 = forward[d1]->link[0];
            su3_link A = contract_2(A1, A2);
            su3_link B1 = forward[d2]->link[0];
            su3_link B2 = forward[d2]->forward[d1]->link[0];
            su3_link B = contract_2(B1, B2);
            su3_link C1 = forward[d2]->forward[d2]->link[0];
            su3_link C2 = forward[d2]->forward[d2]->forward[d1]->link[0];
            su3_link C = contract_2(C1, C2);
            U += contract_3(A, B, C);
        }
    }
    return U.real() / (lattice->D - 1) / (lattice->D - 2);
}

//   *---*---*   ^ d1
//  /   /   /|   |
// *---*---* *   +--> d3
// |   |   |/   /
// *---*---*   L d2

double su3_site::contract_3x2x2() {
    complex<double> U = 0.0;
    su3_site* s;
    for (int d1 = 1; d1 < lattice->D; d1++) {
        for (int d2 = 1; d2 < lattice->D; d2++) {
            for (int d3 = 1; d3 < lattice->D; d3++) {
                if (d1 == d2) continue;
                if (d2 == d3) continue;
                if (d1 == d3) continue;
                
                s = this;
                
                su3_link A11 = s->link[0];
                su3_link A12 = s->forward[d1]->link[0];
                su3_link A1 = contract_2(A11, A12);
                
                s = s->forward[d2];
                su3_link A21 = s->link[0];
                su3_link A22 = s->forward[d1]->link[0];
                su3_link A2 = contract_2(A21, A22);
                
                su3_link A = contract_2(A1, A2);
                
                s = this->forward[d3];
                su3_link B11 = s->link[0];
                su3_link B12 = s->forward[d1]->link[0];
                su3_link B1 = contract_2(B11, B12);
                
                s = s->forward[d2];
                su3_link B21 = s->link[0];
                su3_link B22 = s->forward[d1]->link[0];
                su3_link B2 = contract_2(B21, B22);
                
                su3_link B = contract_2(B1, B2);
                
                s = this->forward[d3]->forward[d3];
                su3_link C11 = s->link[0];
                su3_link C12 = s->forward[d1]->link[0];
                su3_link C1 = contract_2(B11, B12);
                
                s = s->forward[d2];
                su3_link C21 = s->link[0];
                su3_link C22 = s->forward[d1]->link[0];
                su3_link C2 = contract_2(C21, C22);
                
                su3_link C = contract_2(C1, C2);
                
                U += contract_3(A, B, C);
            }
        }
    }
    return U.real() / (lattice->D - 1) / (lattice->D - 2) / (lattice->D - 3);
}

//   *   *---*  ^ d1
//  /|   | B |  |
// *A* *---* *  +-->d3
// |/  | C |   /
// *   *---*  L d2

double su3_site::contract_3x2x2_T() {
    complex<double> U = 0.0;
    su3_site* s;
    for (int d1 = 1; d1 < lattice->D; d1++) {
        for (int d2 = 1; d2 < lattice->D; d2++) {
            for (int d3 = 1; d3 < lattice->D; d3++) {
                if (d1 == d2) continue;
                if (d2 == d3) continue;
                if (d1 == d3) continue;
                
                s = this;
                
                su3_link A11 = s->link[0];
                su3_link A12 = s->forward[d1]->link[0];
                su3_link A1 = contract_2(A11, A12);
                
                s = s->forward[d2];
                su3_link A21 = s->link[0];
                su3_link A22 = s->forward[d1]->link[0];
                su3_link A2 = contract_2(A21, A22);
                
                su3_link A = contract_2(A1, A2);
                
                s = this->forward[d3];
                su3_link B11 = s->link[0];
                su3_link B12 = s->forward[d1]->link[0];
                su3_link B1 = contract_2(B11, B12);
                
                s = s->forward[d3];
                su3_link B21 = s->link[0];
                su3_link B22 = s->forward[d1]->link[0];
                su3_link B2 = contract_2(B21, B22);
                
                su3_link B = contract_2(B1, B2);
                
                s = this->forward[d2]->forward[d3];
                su3_link C11 = s->link[0];
                su3_link C12 = s->forward[d1]->link[0];
                su3_link C1 = contract_2(B11, B12);
                
                s = s->forward[d3];
                su3_link C21 = s->link[0];
                su3_link C22 = s->forward[d1]->link[0];
                su3_link C2 = contract_2(C21, C22);
                
                su3_link C = contract_2(C1, C2);
                
                U += contract_3(A, B, C);
            }
        }
    }
    return U.real() / (lattice->D - 1) / (lattice->D - 2) / (lattice->D - 3);
}

su3_link su3_site::overrelax(su3_link g) {
    
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
    
    // make the determinant 1
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
    lock();
//    for (int d = 0; d < lattice->D; d++) backward[d]->lock();
    su3_link g = make_unitary(su3_identity - sum_G(coulomb) * I / 4);
//    g = overrelax(g);
    
    // apply the gauge transformation (including the time direction)
    for (int d = 0; d < lattice->D; d++) {
        set_link(d, g * link[d]);
        backward[d]->set_link(d, backward[d]->link[d] * g.transpose().conjugate());
    }
    unlock();
//    for (int d = 0; d < lattice->D; d++) backward[d]->unlock();
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
    lock();
    double accept = 0;
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
    this->z = 50;
    this->sigma[0] = su2_identity;
    this->sigma[1] << 0, 1, 1, 0;
    this->sigma[2] << 0, -I, I, 0;
    this->sigma[3] << 1, 0, 0, -1;
    this->n_threads = T;

    random_device rd;  // seed for the random number engine
    gen = mt19937(rd()); // mersenne_twister_engine seeded with rd()
    
    n_sites = pow(N, D - 1) * T;
    site = new su3_site[n_sites];
    
    for (int s = 0; s < n_sites; s++) {
        site[s].lattice = this;
        site[s].thread_id = gen() % n_threads;
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
    this->z = 50;
    this->sigma[0] = su2_identity;
    this->sigma[1] << 0, 1, 1, 0;
    this->sigma[2] << 0, -I, I, 0;
    this->sigma[3] << 1, 0, 0, -1;
    this->n_threads = T;

    random_device rd;  // seed for the random number engine
    gen = mt19937(rd()); // mersenne_twister_engine seeded with rd()
    
    n_sites = pow(N, D - 1) * T;
    site = new su3_site[n_sites];
    
    for (int s = 0; s < n_sites; s++) {
        site[s].lattice = this;
        site[s].thread_id = gen() % n_threads;
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

double su3_lattice::contract_3x1_L() {
    double sum = 0;
    for (int s = 0; s < n_sites; s++) sum += site[s].contract_3x1_L();
    return sum / n_sites;
}

double su3_lattice::contract_3x1() {
    double sum = 0;
    for (int s = 0; s < n_sites; s++) sum += site[s].contract_3x1();
    return sum / n_sites;
}

double su3_lattice::contract_3x2_L() {
    double sum = 0;
    for (int s = 0; s < n_sites; s++) sum += site[s].contract_3x2_L();
    return sum / n_sites;
}

double su3_lattice::contract_3x2() {
    double sum = 0;
    for (int s = 0; s < n_sites; s++) sum += site[s].contract_3x2();
    return sum / n_sites;
}

double su3_lattice::contract_3x2x2() {
    double sum = 0;
    for (int s = 0; s < n_sites; s++) sum += site[s].contract_3x2x2();
    return sum / n_sites;
}

double su3_lattice::contract_3x2x2_T() {
    double sum = 0;
    for (int s = 0; s < n_sites; s++) sum += site[s].contract_3x2x2_T();
    return sum / n_sites;
}

double async_sweep(su3_site* site, int n) {
    double sum = 0.0;
    for (int i = 0; i < n; i++) {
        sum += site->sweep();
        site++;
    }
    return sum;
}

void su3_lattice::sweep(int n_sweeps) {
    
    if (n_sweeps) {
        for (int i = 0; i < n_sweeps; i++) sweep(0);
        return;
    }
    
    future<double> accept_thread[n_threads];
    
    int n_slice = n_sites / n_threads;
    for (int n = 0; n < n_threads; n++) {
        accept_thread[n] = async(async_sweep, &site[n * n_slice], n_slice);
    }

    double accept = 0.0;
    for (int n = 0; n < n_threads; n++) accept += accept_thread[n].get();
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
        
        for (int s = 0; s < n; s++) site[s].relax(coulomb);
        
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
//    double th = 1.0;

    int n_slice = n_sites / T; // number of sites in each time slice
    future<double> avg_link_slice[n_threads];
    for (int n = 0; n < n_threads; n++) avg_link_slice[n] = async(async_relax, &site[n * n_slice], n_slice, error_target / n_threads, coulomb);
    for (int n = 0; n < n_threads; n++) avg_link += avg_link_slice[n].get();
    
    return avg_link / T;

//    if (coulomb) {
//        // for coulomb gauge, relax each time slice separately
//        int n_slice = n_sites / T; // number of sites in each time slice
//        future<double> avg_link_slice[n_threads];
//        for (int n = 0; n < n_threads; n++) avg_link_slice[n] = async(async_relax, &site[n * n_slice], n_slice, error_target / n_threads, coulomb);
//
//        for (int n = 0; n < n_threads; n++) avg_link += avg_link_slice[n].get();
//
//        return avg_link / T;
//    } else {
//        // for landau gauge, relax the entire lattice all at once
//        while (th > error_target) {
//
//            for (int s = 0; s < n_sites; s++) site[s].relax(coulomb);
//
//            th = theta(coulomb);
//
//            // add up all the links
//            avg_link = 0.0;
//            if (coulomb) {
//                for (int s = 0; s < n_sites; s++) avg_link += site[s].sum_coulomb();
//            } else {
//                for (int s = 0; s < n_sites; s++) avg_link += site[s].sum_landau();
//            }
//            avg_link /= n_sites;
////            cout << "avg_link = " << avg_link << ", ";
////            cout << "theta = " << th << ", ";
////            cout << "plaq = " << plaq() << "\n";
//        }
//        return avg_link;
//    }
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
