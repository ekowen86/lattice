//
//  su3_x.cpp
//
//  Created by Evan Owen on 10/19/18.
//  Copyright © 2018 Evan Owen. All rights reserved.
//

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <future>
#include <unsupported/Eigen/MatrixFunctions>

#include "su3_x.hpp"

using namespace std;

// MARK: su3_x_site

void su3_x_site::init(su3_x_lattice* lattice, su3_x_site* lattice_sites, int s) {
    this->lattice = lattice;
    
    gen = &lattice->gen[s / lattice->n_slice];
    forward_edge = (s % lattice->N5) == (lattice->N5 - 1);
    backward_edge = (s % lattice->N5) == 0;
    reset_links(true);
    is_locked = false;
    
    for (int d = 0; d < (lattice->D + 1); d++) {
        forward[d] = &lattice_sites[lattice->get_site(s, d, +1)];
        backward[d] = &lattice_sites[lattice->get_site(s, d, -1)];
    }
}

su3_link su3_x_site::make_unitary(const su3_link& g) {
    
    // create 3-vectors from two random rows
    // the third row is unused
    int r1 = rand_int(0, 2);
    int r2 = (r1 + 1) % 3;
    int r3 = (r1 + 2) % 3;
    su3_vector u = g.row(r1);
    su3_vector v = g.row(r2);

    // normalize u
    u /= u.norm();

    // make v orthonormal to u
    v = v - u * (u.dot(v));
    v /= v.norm();
    
    // create a third orthonormal vector w
    su3_vector w = u.cross(v);
    
    su3_link U;

    U.row(r1) = u;
    U.row(r2) = v;
    U.row(r3) = w;

    return U;
}

su3_link su3_x_site::cayley_ham(const su3_link& Q) {
    // use cayley-hamilton algorithm to compute e^iQ
    // ref. Morningstar, Peardon, Phys. Rev. D 69, 054501 (2004)
    
    su3_link Q2 = Q * Q;
    double c0 = Q.determinant().real();
    double c0abs = abs(c0);
    double c1 = 0.5 * Q2.trace().real();
    if (c1 < 0) return (I * Q).exp(); // check for nan
    double c0_max = 2.0 * pow(c1 / 3.0, 1.5);
    double theta = acos(c0abs / c0_max);
    if (isnan(theta)) return (I * Q).exp(); // check for nan
    double u = sqrt(c1 / 3.0) * cos(theta / 3.0);
    if (isnan(u)) return (I * Q).exp(); // check for nan
    double w = sqrt(c1) * sin(theta / 3.0);
    if (isnan(w)) return (I * Q).exp(); // check for nan
    double u2 = u * u;
    double w2 = w * w;
    double cos_w = cos(w);
    double xi0;
    if (w <= 0.05) {
        xi0 = 1.0 - w2 / 6.0 * (1.0 - w2 / 20.0 * (1.0 - w2 / 42.0));
    } else {
        xi0 = sin(w) / w;
    }
    complex<double> h0 = polar(u2 - w2, 2.0 * u) + polar(8.0 * u2 * cos_w, -u) + polar(2.0 * u * (3.0 * u2 + w2) * xi0, -u + M_PI_2);
    complex<double> h1 = polar(2.0 * u, 2.0 * u) - polar(2.0 * u * cos_w, -u) + polar((3.0 * u2 - w2) * xi0, -u + M_PI_2);
    complex<double> h2 = polar(1.0, 2.0 * u) - polar(cos_w, -u) - polar(3.0 * u * xi0, -u + M_PI_2);
    complex<double> f0 = h0 / (9.0 * u2 - w2);
    complex<double> f1 = h1 / (9.0 * u2 - w2);
    complex<double> f2 = h2 / (9.0 * u2 - w2);
    if (c0 < 0) {
        f0 = conj(f0);
        f1 = -conj(f1);
        f2 = conj(f2);
    }
    if (isnan(real(f0)) || isnan(real(f1)) || isnan(real(f2))) return (I * Q).exp(); // check for nan
    return f0 * su3_identity + f1 * Q + f2 * Q2;
}

bool su3_x_site::lock() {
    if (is_locked || forward[0]->is_locked || backward[0]->is_locked) {
        this_thread::sleep_for(chrono::milliseconds(10));
        return true;
    }
    is_locked = true;
    forward[0]->is_locked = true;
    backward[0]->is_locked = true;
    return false;
}

void su3_x_site::unlock() {
    backward[0]->is_locked = false;
    forward[0]->is_locked = false;
    is_locked = false;
}

double su3_x_site::rand_double(double min, double max) {
    std::uniform_real_distribution<double> dist(min, max);
    return dist(*gen);
}

int su3_x_site::rand_int(int min, int max) {
    std::uniform_int_distribution<int> dist(min, max);
    return dist(*gen);
}

double su3_x_site::rand_normal(double mean, double stdev) {
    std::normal_distribution<double> dist(mean, stdev);
    return dist(*gen);
}

void su3_x_site::reset_links(bool cold) {
    for (int d = 0; d < lattice->D; d++) set_link(d, cold ? su3_identity : create_link(true));

    // always set the link in the extra dimension to the identity
    link[lattice->D] = su3_identity;
    link_inverse[lattice->D] = su3_identity;
}

void su3_x_site::copy_links(su3_x_site* site) {
    for (int d = 0; d < lattice->D; d++) {
        link[d] = su3_link(site->link[d]);
        link_inverse[d] = su3_link(site->link_inverse[d]);
    }

    // always set the link in the extra dimension to the identity
    link[lattice->D] = su3_identity;
    link_inverse[lattice->D] = su3_identity;
}

void su3_x_site::set_link(int d, const su3_link& value) {
    // never change the link in the extra dimension
    if (d == lattice->D) return;
    
    su3_link U = make_unitary(value);
    if (!U.allFinite()) {
        cout << "Non-finite matrix element found" << endl;
        return;
    }
    
    link[d] = U;
    link_inverse[d] = U.adjoint();
}

su2_link su3_x_site::create_su2(bool random) {
    
    double a, g;
    
    if (random) {
        a = rand_double(-1.0, 1.0);
        g = sqrt(1 - a * a);
    } else {
        double x;
        double z = lattice->z;
        
        do {
            x = rand_double(exp(-2.0 * z), 1.0);
            a = 1.0 + log(x) / z;
            g = sqrt(1 - a * a);
            // check for nan
        } while (isnan(g) || (rand_double() < g));
    }
    
    double theta = acos(rand_double(-1.0, 1.0));
    if (isnan(theta)) theta = 0.0;
    double phi = rand_double(0.0, 2.0 * M_PI);
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

su3_link su3_x_site::create_link(bool random) {
    
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
    if (rand_double() > 0.5) return RST.adjoint();
    return RST;
}

double su3_x_site::action() {
    // gauge action
    double S = 0.0;
    su3_link s;
    for (int d1 = 0; d1 < lattice->D; d1++) {
        for (int d2 = d1 + 1; d2 < lattice->D; d2++) {
            s = link[d1];
            s *= forward[d1]->link[d2];
            s *= forward[d2]->link_inverse[d1];
            s *= link_inverse[d2];
            S += (su3_identity - s).trace().real();
        }
    }
    return lattice->beta * S / 3.0;
}

double su3_x_site::hamiltonian() {
    // wilson action (including extra dimension)
    double S = 0.0;
    su3_link s;
    int d_max = forward_edge ? lattice->D : (lattice->D + 1);
    for (int d1 = 0; d1 < d_max; d1++) {
        for (int d2 = d1 + 1; d2 < d_max; d2++) {
            s = link[d1];
            s *= forward[d1]->link[d2];
            s *= forward[d2]->link_inverse[d1];
            s *= link_inverse[d2];
            if (d2 < lattice->D) {
                S += (su3_identity - s).trace().real();
            } else {
                S += (su3_identity - s).trace().real() * lattice->eps5;
            }
        }
    }
    S = lattice->beta * S / 3.0;

    // extra dimension never contributes to kinetic energy
    for (int d = 0; d < lattice->D; d++) S += (p_link[d] * p_link[d]).trace().real();
    return S;
}

void su3_x_site::init_momenta() {
    for (int d = 0; d < lattice->D; d++) {
        p_link[d] = su3_zero;
        for (int l = 0; l < 8; l++) p_link[d] += rand_normal() * lattice->lambda[l] * 0.5;
    }
    // momentum in extra dimension is always zero
    p_link[lattice->D] = su3_zero;
}

void su3_x_site::hmc_step_p(double frac) {
    for (int d = 0; d < lattice->D; d++) p_link[d] = p_link[d] - frac * lattice->dt * p_link_dot(d);
}

void su3_x_site::hmc_step_link() {
    for (int d = 0; d < lattice->D; d++) set_link(d, cayley_ham(lattice->dt * p_link[d]) * link[d]);
}

su3_link su3_x_site::p_link_dot(int d) {

    su3_link A = staple(d);
    su3_link U = link[d];
    su3_link UA = U * A;
    su3_link AU_dag = A.adjoint() * U.adjoint();
    su3_link X = (UA - AU_dag) - (UA - AU_dag).trace() * su3_identity / 3.0;
    return -I * lattice->beta * X / 12.0;
}

double su3_x_site::plaq() {
    // compute the average plaquette at this site
    double U = 0;
    su3_link u;
    for (int d1 = 0; d1 < lattice->D; d1++) {
        for (int d2 = d1 + 1; d2 < lattice->D; d2++) {
            u = link[d1];
            u *= forward[d1]->link[d2];
            u *= forward[d2]->link_inverse[d1];
            u *= link_inverse[d2];
            U += u.trace().real() / 3.0;
        }
    }
    return U / lattice->D / (lattice->D - 1) * 2.0;
}

su3_link su3_x_site::plaquette(int d1, int d2) {
    // compute the plaquette in the d1,d2 direction
    su3_link u1, u2, u3, u4;
    u1 = link[d1];
    u2 = forward[d1]->link[d2];
    u3 = forward[d2]->link_inverse[d1];
    u4 = link_inverse[d2];
    return u1 * u2 * u3 * u4;
}

su3_link su3_x_site::staple(int d1) {
    // compute the sum of staples attached to the plaquette pointing in the d1 direction
    su3_link U = su3_zero;
    su3_link u;
    for (int d2 = 0; d2 < lattice->D; d2++) {
        if (d1 == d2) continue;
        
        // forward staple
        u = forward[d1]->link[d2];
        u *= forward[d2]->link_inverse[d1];
        u *= link_inverse[d2];
        U += u;
        
        // backward staple
        u = backward[d2]->forward[d1]->link_inverse[d2];
        u *= backward[d2]->link_inverse[d1];
        u *= backward[d2]->link[d2];
        U += u;
    }
    
    if (d1 == lattice->D) return U;
    
    if (!forward_edge) {
        // forward staple in extra dimension
        u = forward[d1]->link[lattice->D];
        u *= forward[lattice->D]->link_inverse[d1];
        u *= link_inverse[lattice->D];
        U += u * lattice->eps5;
    }
    
    if (!backward_edge) {
        // backward staple in extra dimension
        u = backward[lattice->D]->forward[d1]->link_inverse[lattice->D];
        u *= backward[lattice->D]->link_inverse[d1];
        u *= backward[lattice->D]->link[lattice->D];
        U += u * lattice->eps5;
    }
    
    return U;
}

su3_link su3_x_site::cloverleaf(int d1, int d2) {
    su3_link U = su3_identity;
    U *= link[d1];
    U *= forward[d1]->link[d2];
    U *= forward[d2]->link_inverse[d1];
    U *= link_inverse[d2];
    
    U *= backward[d2]->link_inverse[d2];
    U *= backward[d2]->link[d1];
    U *= backward[d2]->forward[d1]->link[d2];
    U *= link_inverse[d1];
    
    U *= backward[d1]->link_inverse[d1];
    U *= backward[d1]->backward[d2]->link_inverse[d2];
    U *= backward[d1]->backward[d2]->link[d1];
    U *= backward[d2]->link[d2];
    
    U *= link[d2];
    U *= backward[d1]->forward[d2]->link_inverse[d1];
    U *= backward[d1]->link_inverse[d2];
    U *= backward[d1]->link[d1];
    
    return U;
}

double su3_x_site::wilson_loop(int a, int b) {
    // compute the average a x b wilson loop at this site
    double U = 0;
    su3_link u;
    su3_x_site* s = this;
    int x;

    // don't use the time direction on an asymmetric lattice
    int d_min = (lattice->N == lattice->T) ? 0 : 1;
    int n = 0; // number of measurements
    
    for (int d1 = d_min; d1 < lattice->D; d1++) {
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
            U += u.trace().real() / 3.0;
            n++;
        }
    }
    return U / double(n);
}

complex<double> su3_x_site::polyakov_loop(int d) {
    // compute the polyakov loops at this site in the d direction
    su3_x_site* s = this; // current site
    su3_link u1 = su3_identity;
    int t_max = (d == 0) ? lattice->T : lattice->N;
    
    for (int t = 0; t < t_max; t++) {
        u1 *= s->link[d];
        s = s->forward[d];
    }
    
    return u1.trace();
}

double su3_x_site::correlator(int T) {
    // correlator of duration T
    su3_link u = su3_identity;
    su3_x_site* s = this;
    int x;
    for (x = 0; x < T; x++) {
        u *= s->link[0];
        s = s->forward[0];  // time direction is always 0
    }
    return u.trace().real() / 3;
}

int levi_civita_4[] = {
    0,1,2,3,+1,
    0,1,3,2,-1,
    0,2,1,3,-1,
    0,2,3,1,+1,
    0,3,1,2,+1,
    0,3,2,1,-1,
    1,0,2,3,-1,
    1,0,3,2,+1,
    1,2,0,3,+1,
    1,2,3,0,-1,
    1,3,0,2,-1,
    1,3,2,0,+1,
    2,0,1,3,+1,
    2,0,3,1,-1,
    2,1,0,3,-1,
    2,1,3,0,+1,
    2,3,0,1,+1,
    2,3,1,0,-1,
    3,0,1,2,-1,
    3,0,2,1,+1,
    3,1,0,2,+1,
    3,1,2,0,-1,
    3,2,0,1,-1,
    3,2,1,0,+1
};

double su3_x_site::topological_charge() {

    // symmetric plaquettes
    su3_link u1, u2, u3, u4;
    double q = 0.0;

    for (int i = 0; i < 24; i++) {
        complex<double> C = 0.0;
        int l = i * 5;
        int d1 = levi_civita_4[l++];
        int d2 = levi_civita_4[l++];
        int d3 = levi_civita_4[l++];
        int d4 = levi_civita_4[l++];
        double parity = double(levi_civita_4[l++]);

        u1 = link[d1];
        u2 = forward[d1]->link[d2];
        u3 = forward[d2]->link_inverse[d1];
        u4 = link_inverse[d2];
        su3_link Up1p2 = u1 * u2 * u3 * u4;
        
        u1 = link[d1];
        u2 = forward[d1]->backward[d2]->link_inverse[d2];
        u3 = backward[d2]->link_inverse[d1];
        u4 = backward[d2]->link[d2];
        su3_link Up1n2 = u1 * u2 * u3 * u4;

        u1 = backward[d1]->link_inverse[d1];
        u2 = backward[d1]->link[d2];
        u3 = backward[d1]->forward[d2]->link[d1];
        u4 = link_inverse[d2];
        su3_link Un1p2 = u1 * u2 * u3 * u4;

        u1 = backward[d1]->link_inverse[d1];
        u2 = backward[d1]->backward[d2]->link_inverse[d2];
        u3 = backward[d1]->backward[d2]->link[d1];
        u4 = backward[d2]->link[d2];
        su3_link Un1n2 = u1 * u2 * u3 * u4;

        u1 = link[d3];
        u2 = forward[d3]->link[d4];
        u3 = forward[d4]->link_inverse[d3];
        u4 = link_inverse[d4];
        su3_link Up3p4 = u1 * u2 * u3 * u4;
        
        u1 = link[d3];
        u2 = forward[d3]->backward[d4]->link_inverse[d4];
        u3 = backward[d4]->link_inverse[d3];
        u4 = backward[d4]->link[d4];
        su3_link Up3n4 = u1 * u2 * u3 * u4;

        u1 = backward[d3]->link_inverse[d3];
        u2 = backward[d3]->link[d4];
        u3 = backward[d3]->forward[d4]->link[d3];
        u4 = link_inverse[d4];
        su3_link Un3p4 = u1 * u2 * u3 * u4;

        u1 = backward[d3]->link_inverse[d3];
        u2 = backward[d3]->backward[d4]->link_inverse[d4];
        u3 = backward[d3]->backward[d4]->link[d3];
        u4 = backward[d4]->link[d4];
        su3_link Un3n4 = u1 * u2 * u3 * u4;
        
        C += (Up1p2 * Up3p4).trace();
        C -= (Up1n2 * Up3p4).trace();
        C -= (Un1p2 * Up3p4).trace();
        C += (Un1n2 * Up3p4).trace();

        C -= (Up1p2 * Up3n4).trace();
        C += (Up1n2 * Up3n4).trace();
        C += (Un1p2 * Up3n4).trace();
        C -= (Un1n2 * Up3n4).trace();

        C -= (Up1p2 * Un3p4).trace();
        C += (Up1n2 * Un3p4).trace();
        C += (Un1p2 * Un3p4).trace();
        C -= (Un1n2 * Un3p4).trace();

        C += (Up1p2 * Un3n4).trace();
        C -= (Up1n2 * Un3n4).trace();
        C -= (Un1p2 * Un3n4).trace();
        C += (Un1n2 * Un3n4).trace();
        
        q += real(C) * parity;
    }

    return q / 16.0;
    
//    // single plaquette
//    complex<double> C = 0.0;
//    su2_link C01 = plaquette(0, 1).imag();
//    su2_link C10 = plaquette(1, 0).imag();
//    su2_link C23 = plaquette(2, 3).imag();
//    su2_link C32 = plaquette(3, 2).imag();
//    C += (C01 * C23).trace() - (C10 * C23).trace() - (C01 * C32).trace() + (C10 * C32).trace();
//    C += (C23 * C01).trace() - (C23 * C10).trace() - (C32 * C01).trace() + (C32 * C10).trace();
//    su2_link C02 = plaquette(0, 2).imag();
//    su2_link C20 = plaquette(2, 0).imag();
//    su2_link C13 = plaquette(1, 3).imag();
//    su2_link C31 = plaquette(3, 1).imag();
//    C -= (C02 * C13).trace() - (C20 * C13).trace() - (C02 * C31).trace() + (C20 * C31).trace();
//    C -= (C13 * C02).trace() - (C13 * C20).trace() - (C31 * C02).trace() + (C31 * C20).trace();
//    su2_link C03 = plaquette(0, 3).imag();
//    su2_link C30 = plaquette(3, 0).imag();
//    su2_link C12 = plaquette(1, 2).imag();
//    su2_link C21 = plaquette(2, 1).imag();
//    C += (C03 * C12).trace() - (C30 * C12).trace() - (C03 * C21).trace() + (C30 * C21).trace();
//    C += (C12 * C03).trace() - (C12 * C30).trace() - (C21 * C03).trace() + (C21 * C30).trace();
//    return C.real();

//    // cloverleaf
//    su3_link C01 = cloverleaf(0, 1).imag();
//    su3_link C02 = cloverleaf(0, 2).imag();
//    su3_link C03 = cloverleaf(0, 3).imag();
//    su3_link C10 = cloverleaf(1, 0).imag();
//    su3_link C12 = cloverleaf(1, 2).imag();
//    su3_link C13 = cloverleaf(1, 3).imag();
//    su3_link C20 = cloverleaf(2, 0).imag();
//    su3_link C21 = cloverleaf(2, 1).imag();
//    su3_link C23 = cloverleaf(2, 3).imag();
//    su3_link C30 = cloverleaf(3, 0).imag();
//    su3_link C31 = cloverleaf(3, 1).imag();
//    su3_link C32 = cloverleaf(3, 2).imag();
//
//    complex<double> C = 0.0;
//    C += (C01 * C23).trace() - (C01 * C32).trace();
//    C += (C02 * C31).trace() - (C02 * C13).trace();
//    C += (C03 * C12).trace() - (C03 * C21).trace();
//    C += (C12 * C03).trace() - (C12 * C30).trace();
//    C += (C13 * C20).trace() - (C13 * C02).trace();
//    C += (C10 * C32).trace() - (C10 * C23).trace();
//    C += (C23 * C01).trace() - (C23 * C10).trace();
//    C += (C20 * C13).trace() - (C20 * C31).trace();
//    C += (C21 * C30).trace() - (C21 * C03).trace();
//    C += (C30 * C21).trace() - (C30 * C12).trace();
//    C += (C31 * C02).trace() - (C31 * C20).trace();
//    C += (C32 * C10).trace() - (C32 * C01).trace();
//
//    return C.real();
}

double su3_x_site::mag_U() {
    return this->link[0].trace().real() / 3;
}

double su3_x_site::abs_U() {
    return abs(this->link[0].trace().real()) / 3;
}

double su3_x_site::four_point(int T, int R) {
    // four-point function of duration T and spacing R
    double U = 0;
    su3_link u;
    su3_x_site* s = this;
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

su3_link su3_x_site::overrelax(su3_link g) {
    
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

void su3_x_site::relax(bool coulomb) {
    while (!lock()) {}
    su3_link g = make_unitary(su3_identity - sum_G(coulomb) * I / 4);
//    g = overrelax(g);
    
    // apply the gauge transformation (including the time direction)
    for (int d = 0; d < lattice->D; d++) {
        set_link(d, g * link[d]);
        backward[d]->set_link(d, backward[d]->link[d] * g.adjoint());
    }
    unlock();
}

double su3_x_site::sum_landau() {
    double U = 0.0;
    for (int d = 0; d < lattice->D; d++) U += link[d].trace().real();
    return U / 3 / lattice->D;
}

double su3_x_site::sum_coulomb() {
    double U = 0.0;
    for (int d = 1; d < lattice->D; d++) U += link[d].trace().real();
    return U / 3 / (lattice->D - 1);
}

su3_link su3_x_site::sum_G(bool coulomb) {
    su3_link G = su3_zero;
    for (int d = coulomb ? 1 : 0; d < lattice->D; d++) {
        G += (link[d] - link_inverse[d]) / 2 / I;
        G -= (backward[d]->link[d] - backward[d]->link_inverse[d]) / 2 / I;
    }
    return G;
}

double su3_x_site::heat_bath() {
    double accept = 0;
    while (lock());
    for (int d = 0; d < lattice->D; d++) accept += heat_bath_link(d);
    unlock();
    return accept / lattice->D;
}

double su3_x_site::heat_bath_link(int d1) {
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
    
    if (!forward_edge) {
        u = forward[d1]->link[lattice->D];
        u *= forward[lattice->D]->link_inverse[d1];
        u *= link_inverse[lattice->D];
        A += u;
    }
    
    if (!backward_edge) {
        u = backward[lattice->D]->forward[d1]->link_inverse[lattice->D];
        u *= backward[lattice->D]->link_inverse[d1];
        u *= backward[lattice->D]->link[lattice->D];
        A += u;
    }
    
    double dS = -lattice->beta * ((new_link - old_link) * A).trace().real() / 3.0;
    
    // metropolis algorithm
    if (dS <= 0 || rand_double() < exp(-dS)) {
        set_link(d1, new_link);
        return 1.0;
    }
    return 0.0;
}

void su3_x_site::cool() {
    lock();
    for (int d = 0; d < lattice->D; d++) cool_link(d);
    unlock();
}

void su3_x_site::cool_link(int d1) {
    
    su3_link X = su3_zero;
    su3_link u1, u2, u3;
    for (int d2 = 0; d2 < lattice->D; d2++) {
        if (d1 == d2) continue;
        
        // forward staple
        u1 = link[d2];
        u2 = forward[d2]->link[d1];
        u3 = forward[d1]->link_inverse[d2];
        X += u1 * u2 * u3;
        
        // backward staple
        u1 = backward[d2]->link_inverse[d2];
        u2 = backward[d2]->link[d1];
        u3 = forward[d1]->backward[d2]->link[d2];
        X += u1 * u2 * u3;
    }
    set_link(d1, X);
}

void su3_x_site::wilson_flow(su3_x_site* target, double epsilon) {
    lock();
    for (int d = 0; d < lattice->D; d++) wilson_flow_link(target, epsilon, d);
    unlock();
}

void su3_x_site::wilson_flow_link(su3_x_site* target, double epsilon, int d1) {
    
    su3_link X = su3_zero;
    su3_link u1, u2, u3;
    for (int d2 = 0; d2 < lattice->D; d2++) {
        if (d1 == d2) continue;
        
        // forward staple
        u1 = link[d2];
        u2 = forward[d2]->link[d1];
        u3 = forward[d1]->link_inverse[d2];
        X += u1 * u2 * u3;
        
        // backward staple
        u1 = backward[d2]->link_inverse[d2];
        u2 = backward[d2]->link[d1];
        u3 = forward[d1]->backward[d2]->link[d2];
        X += u1 * u2 * u3;
    }
    X = X.adjoint().eval();

    su3_link W0 = link[d1];
    su3_link Z0 = W0 * X;
    Z0 = (Z0.adjoint() - Z0).eval();
    Z0 = 0.5 * Z0 - Z0.trace() * su3_identity / 6.0;

    su3_link W1 = cayley_ham(-I * 0.25 * Z0 * epsilon) * W0;
    su3_link Z1 = W1 * X;
    Z1 = (Z1.adjoint() - Z1).eval();
    Z1 = 0.5 * Z1 - Z1.trace() * su3_identity / 6.0;
    Z1 = (8.0/9.0) * Z1 - (17.0/36.0) * Z0;

    su3_link W2 = cayley_ham(-I * Z1 * epsilon) * W1;
    su3_link Z2 = W2 * X;
    Z2 = (Z2.adjoint() - Z2).eval();
    Z2 = 0.5 * Z2 - Z2.trace() * su3_identity / 6.0;
    Z2 = 0.75 * Z2 - Z1;

    su3_link W3 = cayley_ham(-I * Z2 * epsilon) * W2;
    target->set_link(d1, W3);
}

void su3_x_site::stout_smear(su3_x_site* target, double rho) {
    for (int d = 0; d < lattice->D; d++) stout_smear_link(target, rho, d);
}

void su3_x_site::stout_smear_link(su3_x_site* target, double rho, int d1) {
    su3_link C = su3_zero;
    su3_link u1, u2, u3;
    for (int d2 = 0; d2 < lattice->D; d2++) {
        if (d1 == d2) continue;
        
        // forward staple
        u1 = link[d2];
        u2 = forward[d2]->link[d1];
        u3 = forward[d1]->link_inverse[d2];
        C += u1 * u2 * u3;
        
        // backward staple
        u1 = backward[d2]->link_inverse[d2];
        u2 = backward[d2]->link[d1];
        u3 = forward[d1]->backward[d2]->link[d2];
        C += u1 * u2 * u3;
    }
    
    su3_link X = rho * C * link_inverse[d1];
    su3_link X_dag = X.adjoint();
    su3_link Z = X_dag - X;
    su3_link W = 0.5 * I * (Z - su3_identity * Z.trace() / 3.0);
    target->set_link(d1, cayley_ham(W) * link[d1]);
}

// MARK: su3_x_lattice

su3_x_lattice::su3_x_lattice(int N, int T, int N5, int D, double beta, double eps5, bool cold) {
    this->N = N;
    this->T = T;
    this->N5 = N5;
    this->D = D;
    this->beta = beta;
    this->eps5 = eps5;
    this->parallel = false;
    this->n_steps = 30;
    this->dt = 0.01;
    this->verbose = 0;

    init();

    for (int s = 0; s < n_sites; s++) site[s].reset_links(cold);

}

su3_x_lattice::su3_x_lattice(su3_x_lattice* lattice) {
    this->N = lattice->N;
    this->T = lattice->T;
    this->N5 = lattice->N5;
    this->D = lattice->D;
    this->beta = lattice->beta;
    this->eps5 = lattice->eps5;
    this->parallel = lattice->parallel;
    this->n_steps = lattice->n_steps;
    this->dt = lattice->dt;
    this->verbose = lattice->verbose;

    init();
    
    for (int s = 0; s < n_sites; s++) site[s].copy_links(&lattice->site[s]);
}

su3_x_lattice::~su3_x_lattice() {
}

void su3_x_lattice::init() {
    
    // initialize hmc and heat bath parameters
    z = 50;
    hmc_count = 0.0;
    hmc_accept = 0.0;

    // initialize pauli matrices
    sigma.resize(4);
    sigma[0] << 1, 0, 0, 1;
    sigma[1] << 0, 1, 1, 0;
    sigma[2] << 0, -I, I, 0;
    sigma[3] << 1, 0, 0, -1;

    // initialize gell-mann matrices
    lambda.resize(8);
    lambda[0] << 0,1,0, 1,0,0, 0,0,0;
    lambda[1] << 0,-I,0, I,0,0, 0,0,0;
    lambda[2] << 1,0,0, 0,-1,0, 0,0,0;
    lambda[3] << 0,0,1, 0,0,0, 1,0,0;
    lambda[4] << 0,0,-I, 0,0,0, I,0,0;
    lambda[5] << 0,0,0, 0,0,1, 0,1,0;
    lambda[6] << 0,0,0, 0,0,-I, 0,I,0;
    lambda[7] << SQRT1_3,0,0, 0,SQRT1_3,0, 0,0,-2.0*SQRT1_3;

    // initialize rng
    random_device rd;  // seed for the random number engine
    for (int t = 0; t < (T + 1); t++) gen.push_back(mt19937(rd())); // mersenne_twister_engine seeded with rd()

    // initialize sites
    n_sites = pow(N, D - 1) * T * N5;
    n_sites_5 = n_sites / N5;
    n_slice = n_sites / T;
    n5_center = (N5 - 1) / 2;
    site.resize(n_sites);
    hmc_site.resize(n_sites);

    for (int s = 0; s < n_sites; s++) {
        site[s].init(this, site.data(), s);
        hmc_site[s].init(this, hmc_site.data(), s);
    }
}

double su3_x_lattice::rand_double(double min, double max) {
    std::uniform_real_distribution<double> dist(min, max);
    return dist(gen[T]);
}

int su3_x_lattice::rand_int(int min, int max) {
    std::uniform_int_distribution<int> dist(min, max);
    return dist(gen[T]);
}

double su3_x_lattice::rand_normal(double mean, double stdev) {
    std::normal_distribution<double> dist(mean, stdev);
    return dist(gen[T]);
}

int su3_x_lattice::get_site(int s, int d, int n) {
    
    // get the site's coordinates
    int x[D + 1];
    x[0] = s / n_slice;
    for (int d = 1; d < D; d++) x[d] = (s / (int)pow(N, D - d - 1) / N5) % N;
    x[D] = s % N5;

    // go n steps in the d direction
    if (d == 0) {
        // time direction
        if (n > 0) {
            x[0] = (x[0] + n) % T;
        } else {
            x[0] = (x[0] + n + T) % T;
        }
    } else if (d == D) {
        // extra direction
        if (n > 0) {
            x[d] = (x[d] + n) % N5;
        } else {
            x[d] = (x[d] + n + N5) % N5;
        }
    } else {
        // spatial direction
        if (n > 0) {
            x[d] = (x[d] + n) % N;
        } else {
            x[d] = (x[d] + n + N) % N;
        }
    }

    // get the index of the new site
    s = 0;
    for (int d = 0; d < D; d++) s += x[d] * (int)pow(N, D - d - 1) * N5;
    s += x[D];
    return s;
}

double su3_x_lattice::plaq(int n5) {
    // compute the average plaquette of the lattice (normalized to 3)
    double P = 0;
    for (int s = n5; s < n_sites; s += N5) P += site[s].plaq();
    return P / double(n_sites_5) * 3.0;
}

double su3_x_lattice::action(int n5) {
    double S = 0;
    for (int s = n5; s < n_sites; s += N5) S += site[s].action();
    return S / double(n_sites_5);
}

double async_wilson_loop(su3_x_site* site, int n, int a, int b, int N5) {
    double sum = 0.0;
    for (int s = 0; s < n; s += N5) sum += site[s].wilson_loop(a, b);
    return sum;
}

double su3_x_lattice::wilson_loop(int a, int b, int n5) {
    // compute the average a x b wilson loop of the lattice
    double sum = 0;
    
    if (parallel) {
        future<double> async_slice[T];
        for (int t = 0; t < T; t++) async_slice[t] = async(async_wilson_loop, &site[t * n_slice + n5], n_slice, a, b, N5);
        for (int t = 0; t < T; t++) sum += async_slice[t].get();
    } else {
        for (int s = n5; s < n_sites; s += N5) sum += site[s].wilson_loop(a, b);
    }
    return sum / double(n_sites_5);
}

double su3_x_lattice::polyakov_loop(int R, int n5) {
    // compute the average product of polyakov loops spaced r apart
    
    // on a square lattice, any direction can be the time direction
    int d_max = (T == N) ? (D - 1) : 1;
    int n = 0;
    
    complex<double> p[N];
    complex<double> P, p1, p2;
    su3_x_site* s;
    for (int d1 = 0; d1 < d_max; d1++) { // time direction
        for (int d2 = d1 + 1; d2 < D; d2++) { // space directions

            // get polyakov loops
            s = &site[n5]; // start at (origin,n5)
            for (int r = 0; r < N; r++) {
                p[r] = s->polyakov_loop(d1);
                s = s->forward[d2];
            }
            for (int r = 0; r < N; r++) {
                p1 = p[r];
                p2 = p[(r + R) % N];
                P += p1 * conj(p2);
                n++;
            }
        }
    }

    return real(P) / double(n);
}

double su3_x_lattice::correlator(int T, int n5) {
    // compute the average correlator of length T
    double sum = 0;
    for (int s = n5; s < n_sites; s += N5) sum += site[s].correlator(T);
    return sum / double(n_sites_5);
}

double su3_x_lattice::four_point(int T, int R, int n5) {
    // compute the average T x R 4-point function of the lattice
    double sum = 0;
    for (int s = n5; s < n_sites; s += N5) sum += site[s].four_point(T, R);
    return sum / double(n_sites_5);
}

double su3_x_lattice::hamiltonian() {
    double H = 0;
    for (int s = 0; s < n_sites; s++) H += hmc_site[s].hamiltonian();
    return H;
}

void async_init_momenta(su3_x_site* site, int n) {
    for (int s = 0; s < n; s++) site[s].init_momenta();
}

void async_hmc_step_p(su3_x_site* site, int n, double frac) {
    for (int s = 0; s < n; s++) site[s].hmc_step_p(frac);
}

void async_hmc_step_link(su3_x_site* site, int n) {
    for (int s = 0; s < n; s++) site[s].hmc_step_link();
}

void su3_x_lattice::hmc(int n_sweeps, bool update_dt, bool no_metropolis) {
    
    if (n_sweeps) {
        for (int i = 0; i < n_sweeps; i++) hmc(0, update_dt, no_metropolis);
        return;
    }
    
    // copy lattice sites
    for (int s = 0; s < n_sites; s++) hmc_site[s].copy_links(&site[s]);

    // init conjugate momenta
    if (parallel) {
        future<void> async_slice[T];
        for (int t = 0; t < T; t++) async_slice[t] = async(async_init_momenta, &hmc_site[t * n_slice], n_slice);
        for (int t = 0; t < T; t++) async_slice[t].get();
    } else {
        for (int s = 0; s < n_sites; s++) hmc_site[s].init_momenta();
    }

    // calculate old hamiltonian
    double oldH = hamiltonian();
    
    // do the initial half-step
    if (parallel) {
        future<void> async_slice[T];
        for (int t = 0; t < T; t++) async_slice[t] = async(async_hmc_step_p, &hmc_site[t * n_slice], n_slice, 0.5);
        for (int t = 0; t < T; t++) async_slice[t].get();
    } else {
        for (int s = 0; s < n_sites; s++) hmc_site[s].hmc_step_p(0.5);
    }

    // do leapfrog steps
    for (int step = 0; step < (n_steps - 1); step++) {
        if (parallel) {
            future<void> async_slice[T];
            for (int t = 0; t < T; t++) async_slice[t] = async(async_hmc_step_link, &hmc_site[t * n_slice], n_slice);
            for (int t = 0; t < T; t++) async_slice[t].get();
            for (int t = 0; t < T; t++) async_slice[t] = async(async_hmc_step_p, &hmc_site[t * n_slice], n_slice, 1.0);
            for (int t = 0; t < T; t++) async_slice[t].get();
        } else {
            for (int s = 0; s < n_sites; s++) hmc_site[s].hmc_step_link();
            for (int s = 0; s < n_sites; s++) hmc_site[s].hmc_step_p(1.0);
        }
    }

    // do the final half-step
    if (parallel) {
        future<void> async_slice[T];
        for (int t = 0; t < T; t++) async_slice[t] = async(async_hmc_step_link, &hmc_site[t * n_slice], n_slice);
        for (int t = 0; t < T; t++) async_slice[t].get();
        for (int t = 0; t < T; t++) async_slice[t] = async(async_hmc_step_p, &hmc_site[t * n_slice], n_slice, 0.5);
        for (int t = 0; t < T; t++) async_slice[t].get();
    } else {
        for (int s = 0; s < n_sites; s++) hmc_site[s].hmc_step_link();
        for (int s = 0; s < n_sites; s++) hmc_site[s].hmc_step_p(0.5);
    }

    // calculate new hamiltonian
    double newH = hamiltonian();
    double dH = newH - oldH;

    // accept with metropolis algorithm
    bool accept = false;
    if (dH <= 0 || rand_double() < exp(-dH)) {
        accept = true;
        hmc_accept++;
    }
    
    if (accept || no_metropolis) {
        for (int s = 0; s < n_sites; s++) site[s].copy_links(&hmc_site[s]);
    }
    
    hmc_count++;
    if (hmc_count >= 20) {
        // update hmc statistics every 20 iterations
        double accept_rate = double(hmc_accept) / double(hmc_count);
        hmc_count = 0;
        hmc_accept = 0;
        
        // target a 70%-90% acceptance rate
        if (update_dt) {
            if (accept_rate > 0.91) {
                dt *= 1.1;
            } else if (accept_rate < 0.35) {
                dt *= 0.5;
            } else if (accept_rate < 0.69) {
                dt *= 0.9;
            }
        }
        if (verbose) {
            cout << setprecision(6) << fixed;
            cout << "accept_rate = " << accept_rate;
            cout << ", dt = " << dt;
            cout << ", plaq = " << plaq(n5_center) << endl;
        }
    }
}

double async_heat_bath(su3_x_site* site, int n) {
    double accept = 0.0;

    // start at a random spot in the time slice to mitigate parallel effects
    int start = site->rand_int(0, n - 1);
    int end = start + n;
    
    for (int s = start; s < end; s++) accept += site[s % n].heat_bath();

    return accept;
}

void su3_x_lattice::heat_bath(int n_sweeps) {
    
    if (n_sweeps) {
        for (int i = 0; i < n_sweeps; i++) heat_bath(0);
        return;
    }
    
    double accept = 0;
    
    if (parallel) {
        future<double> async_slice[T];
        for (int t = 0; t < T; t += 2) async_slice[t] = async(async_heat_bath, &site[t * n_slice], n_slice);
        for (int t = 0; t < T; t += 2) accept += async_slice[t].get();
        for (int t = 1; t < T; t += 2) async_slice[t] = async(async_heat_bath, &site[t * n_slice], n_slice);
        for (int t = 1; t < T; t += 2) accept += async_slice[t].get();
    } else {
        for (int s = 0; s < n_sites; s++) accept += site[s].heat_bath();
    }
    accept /= n_sites;

    if (verbose) {
        cout << setprecision(6);
        cout << "plaq = " << plaq(n5_center) << endl;
    }

    if (accept < 0.78) {
        z *= 1.1;
    } else if (accept > 0.82) {
        z *= 0.9;
    }
}

void async_cool(su3_x_site* site, int n, int N5) {
    for (int s = 0; s < n; s += N5) site[s].cool();
}

void su3_x_lattice::cool(int n5, int n_sweeps) {
    
    if (n_sweeps) {
        for (int i = 0; i < n_sweeps; i++) cool(n5, 0);
        return;
    }
    
    if (parallel) {
        future<void> async_slice[T];
        for (int t = 0; t < T; t += 2) async_slice[t] = async(async_cool, &site[t * n_slice + n5], n_slice, N5);
        for (int t = 0; t < T; t += 2) async_slice[t].get();
        for (int t = 1; t < T; t += 2) async_slice[t] = async(async_cool, &site[t * n_slice + n5], n_slice, N5);
        for (int t = 1; t < T; t += 2) async_slice[t].get();
    } else {
        for (int s = n5; s < n_sites; s += N5) site[s].cool();
    }
}

void async_wilson_flow(su3_x_site* site, su3_x_site* target, int n, int N5, double epsilon) {
    for (int s = 0; s < n; s += N5) site[s].wilson_flow(&target[s], epsilon);
}

void su3_x_lattice::wilson_flow(int n5, double epsilon, int n_sweeps) {
    if (n_sweeps) {
        for (int i = 0; i < n_sweeps; i++) wilson_flow(n5, epsilon, 0);
        return;
    }
    
    if (parallel) {
        future<void> async_slice[T];
        for (int t = 0; t < T; t++) async_slice[t] = async(async_wilson_flow, &site[t * n_slice + n5], &hmc_site[t * n_slice + n5], n_slice, N5, epsilon);
        for (int t = 0; t < T; t++) async_slice[t].get();
    } else {
        for (int s = n5; s < n_sites; s += N5) site[s].wilson_flow(&hmc_site[s], epsilon);
    }

    // copy lattice sites
    for (int s = 0; s < n_sites; s++) site[s].copy_links(&hmc_site[s]);
}

void async_stout_smear(su3_x_site* site, su3_x_site* target, int n, int N5, double rho) {
    for (int s = 0; s < n; s += N5) site[s].stout_smear(&target[s], rho);
}

void su3_x_lattice::stout_smear(int n5, double rho, int n_sweeps) {
    if (n_sweeps) {
        for (int i = 0; i < n_sweeps; i++) stout_smear(n5, rho, 0);
        return;
    }
    
    if (parallel) {
        future<void> async_slice[T];
        for (int t = 0; t < T; t++) async_slice[t] = async(async_stout_smear, &site[t * n_slice + n5], &hmc_site[t * n_slice + n5], n_slice, N5, rho);
        for (int t = 0; t < T; t++) async_slice[t].get();
    } else {
        for (int s = n5; s < n_sites; s += N5) site[s].stout_smear(&hmc_site[s], rho);
    }

    // copy lattice sites
    for (int s = 0; s < n_sites; s++) site[s].copy_links(&hmc_site[s]);
}

double async_topological_charge(su3_x_site* site, int n, int n5) {
    double Q = 0.0;
    for (int s = 0; s < n; s += n5) Q += site[s].topological_charge();
    return Q;
}

double su3_x_lattice::topological_charge(int n5) {
    double Q = 0.0;

    if (parallel) {
        future<double> sweep_slice[T];
        for (int t = 0; t < T; t++) sweep_slice[t] = async(async_topological_charge, &site[t * n_slice + n5], n_slice, N5);
        for (int t = 0; t < T; t++) Q += sweep_slice[t].get();

    } else {
        for (int s = n5; s < n_sites; s += N5) Q += site[s].topological_charge();
    }
    return (Q / 32.0 / M_PI / M_PI);
}

int su3_x_lattice::thermalize(int n_min, int n_max) {
    
    // create a cold start model
    su3_x_lattice coldStart = su3_x_lattice(N, T, N5, D, beta, eps5, true);
    coldStart.parallel = parallel;
    coldStart.verbose = verbose;
    double hot_plaquette = plaq(n5_center);
    double cold_plaquette = coldStart.plaq(n5_center);
//    cout << "P_hot = ";
//    for (int n = 0; n < N5; n++) cout << plaq(n) << ", ";
//    cout << "P_cold = " << cold_plaquette << endl;
    int m;
    for (m = 0; hot_plaquette < cold_plaquette; m++) {
        // sweep both models until the average plaquette reaches the same value
//        heat_bath();
//        coldStart.heat_bath();
        hmc(0, true, true);
        coldStart.hmc(0, true, true);
        hot_plaquette = plaq(n5_center);
        cold_plaquette = coldStart.plaq(n5_center);
//        cout << "dP = " << abs(cold_plaquette - hot_plaquette) << ", P_hot = ";
//        for (int n = 0; n < N5; n++) cout << plaq(n) << ", ";
//        cout << "P_cold = " << cold_plaquette << endl;
//        if (n_max && (m > n_max)) break;
    }
    
    // sweep a few more times
    if (verbose) {
        cout << "hot_plaquette = " << hot_plaquette;
        cout << ", cold_plaquette = " << cold_plaquette << endl;
    }
    int n_therm = max(n_min, m * 2);
//    heat_bath(n_therm - m);
    hmc(n_therm - m, true, false);

    // return value is total number of sweeps
    if (verbose) cout << "n_therm = " << n_therm << endl;
    return n_therm;
}

double su3_x_lattice::ave_link_t(int n5) {
    // average trace of links in the t direction
    // this should be ~0 in Coulomb gauge
    
    double U = 0;
    
    for (int s = n5; s < n_sites; s += N5) U += site[s].link[0].trace().real() / 3;
    
    return U / double(n_sites);
}

double su3_x_lattice::theta(bool coulomb, int n5) {
    su3_link G1, G2;
    double th = 0.0;
    
    for (int s = n5; s < n_sites; s += N5) {
        G1 = site[s].sum_G(coulomb);
        G1 -= (G1.trace() * su3_identity) / 3;
        G2 = G1.adjoint();
        th += (G1 * G2).trace().real() / 3;
    }
    
    return th;
}

double async_relax(su3_x_site* site, int n, double error_target, bool coulomb, int N5) {
    double avg_link = 0.0;
    double th = 1.0;
    while (th > error_target) {
        
        // start at a random spot in the time slice to mitigate parallel effects
        int start = site->rand_int(0, n - 1);
        int end = start + n;
        
        for (int s = start; s < end; s++) site[s % n].relax(coulomb);
        
        // add up all the links
        su3_link G1, G2;
        th = 0.0;
        for (int s = 0; s < n; s += N5) {
            G1 = site[s].sum_G(coulomb);
            G1 -= (G1.trace() * su3_identity) / 3;
            G2 = G1.adjoint();
            th += (G1 * G2).trace().real() / 3;
        }
        avg_link = 0.0;
        for (int s = 0; s < n; s++) avg_link += site[s].sum_coulomb();
        avg_link /= (n * N5);
    }
    return avg_link;
}

long double su3_x_lattice::relax(long double error_target, bool coulomb, int n5) {
    
    double avg_link = 0.0;
    double th = 1.0;

    if (parallel) {
        future<double> avg_link_slice[T];
        for (int t = 0; t < T; t++)
            avg_link_slice[t] = async(async_relax, &site[t * n_slice + n5], n_slice, error_target / T, coulomb, N5);
        
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
                    G2 = G1.adjoint();
                    th += (G1 * G2).trace().real() / 3;
                }
                avg_link_slice = 0.0;
                for (int s = s1; s < s2; s++) avg_link_slice += site[s].sum_coulomb();
                avg_link_slice /= n_slice;
//                cout << "avg_link = " << avg_link_slice << ", ";
//                cout << "theta = " << th << endl;
            }
            avg_link += avg_link_slice;
        }
        return avg_link / T;
    } else {
        // for landau gauge, relax the entire lattice all at once
        while (th > error_target) {
            
            for (int s = 0; s < n_sites; s++) site[s].relax(coulomb);
            
            th = theta(coulomb, n5);

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
//            cout << "plaq = " << plaq() << endl;
        }
        return avg_link;
    }
}

void su3_x_lattice::write_four_point(const char* filename, bool coulomb, int n5) {
    
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
            file << "\t" << four_point(T, R, n5);
        }
    }
    file << "\n";
    file.close();
}

void su3_x_lattice::write_correlator(const char* filename, bool coulomb, int n5) {
    
    ofstream file;
    file.open(filename, ofstream::app);
    
    int t_prop = T / 2;
    file << setprecision(12);
    file << N << "^" << D - 1 << "*" << T << "\t";
    file << beta << "\t";
    file << (coulomb ? "C*" : "L*");
    for (int T = 1; T <= t_prop; T++) {
        file << "\t" << correlator(T, n5);
    }
    file << "\n";
    file.close();
}
