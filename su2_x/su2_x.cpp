//
//  su2_x.cpp
//
//  Created by Evan Owen on 4/3/21.
//  Copyright © 2021 Evan Owen. All rights reserved.
//

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <future>
#include <unsupported/Eigen/MatrixFunctions>

#include "su2_x.hpp"

using namespace std;

// MARK: su2_x_site

void su2_x_site::init(su2_x_lattice* lattice, su2_x_site* lattice_sites, int s) {
    this->lattice = lattice;

    int t = s / lattice->n_slice;
    gen = &lattice->gen[t];
    int n5 = s % lattice->N5;
//    eps = 1.0 - abs(lattice->n5_center - n5) * lattice->eps5;
    eps = 1.0;
    forward_edge = n5 == (lattice->N5 - 1);
    backward_edge = n5 == 0;
    reset_links(true);
    is_locked = false;

    for (int d = 0; d < (lattice->D + 1); d++) {
        forward[d] = &lattice_sites[lattice->get_site(s, d, +1)];
        backward[d] = &lattice_sites[lattice->get_site(s, d, -1)];
    }
}

su2_link su2_x_site::make_unitary(const su2_link& g) {

    double det = g.determinant().real();
    return g / det;

}

// compute the exponential of an anti-hermitian traceless 2x2 matrix

su2_link su2_x_site::exp_su2(const su2_link& g) {

    double a1 = imag(g(0,1));
    double a2 = real(g(0,1));
    double a3 = imag(g(0,0));

    double alpha = sqrt(a1 * a1 + a2 * a2 + a3 * a3);
    a1 /= alpha;
    a2 /= alpha;
    a3 /= alpha;

    su2_link m;
    m(0,0) = complex<double>(cos(alpha), a3 * sin(alpha));
    m(0,1) = complex<double>(a2 * sin(alpha), a1 * sin(alpha));
    m(1,0) = complex<double>(-a2 * sin(alpha), a1 * sin(alpha));
    m(1,1) = complex<double>(cos(alpha), -a3 * sin(alpha));
    return m;
}

bool su2_x_site::lock() {
    if (is_locked || forward[0]->is_locked || backward[0]->is_locked) {
        this_thread::sleep_for(chrono::milliseconds(10));
        return true;
    }
    is_locked = true;
    forward[0]->is_locked = true;
    backward[0]->is_locked = true;
    return false;
}

void su2_x_site::unlock() {
    backward[0]->is_locked = false;
    forward[0]->is_locked = false;
    is_locked = false;
}

double su2_x_site::rand_double(double min, double max) {
    std::uniform_real_distribution<double> dist(min, max);
    return dist(*gen);
}

int su2_x_site::rand_int(int min, int max) {
    std::uniform_int_distribution<int> dist(min, max);
    return dist(*gen);
}

double su2_x_site::rand_normal(double mean, double stdev) {
    std::normal_distribution<double> dist(mean, stdev);
    return dist(*gen);
}

void su2_x_site::reset_links(bool cold) {
    for (int d = 0; d < lattice->D; d++) {
        set_link(d, cold ? su2_identity : create_random_link());
    }

    // always set the link in the extra dimension to the identity
    link[lattice->D] = su2_identity;
    link_inverse[lattice->D] = su2_identity;
}

void su2_x_site::copy_links(su2_x_site* site) {
    for (int d = 0; d < lattice->D; d++) {
        link[d] = su2_link(site->link[d]);
        link_inverse[d] = su2_link(site->link_inverse[d]);
    }

    // always set the link in the extra dimension to the identity
    link[lattice->D] = su2_identity;
    link_inverse[lattice->D] = su2_identity;
}

void su2_x_site::read_links(ifstream& ckptFile, bool bigEndian) {

    for (int d = 0; d < lattice->D; d++) {
        // read real and imaginary parts of SU(2) matrix (8 doubles)
        double rawValues[8];
        if (bigEndian) {
            char bytes1[8];
            char bytes2[8];
            for (int i = 0; i < 8; i++) {
                ckptFile.read(bytes1, 8);
                // reverse byte order
                for (int j = 0; j < 8; j++) {
                    bytes2[j] = bytes1[7 - j];
                }
                memcpy((char*)(&rawValues[i]), bytes2, 8);
            }
        } else {
            ckptFile.read((char*)rawValues, 18 * 8);
        }
        complex<double> u11(rawValues[0], rawValues[1]);
        complex<double> u12(rawValues[2], rawValues[3]);
        complex<double> u21(rawValues[4], rawValues[5]);
        complex<double> u22(rawValues[6], rawValues[7]);

        su2_link U;
        U << u11, u12,
             u21, u22;

        // NERSC format stores lorentz indices backwards from how i do it
//        set_link(lattice->D - d - 1, U);
        set_link(d, U);
    }

    // always set the link in the extra dimension to the identity
    link[lattice->D] = su2_identity;
    link_inverse[lattice->D] = su2_identity;
}

void su2_x_site::set_link(int d, const su2_link& value) {
    // never change the link in the extra dimension
    if (d == lattice->D) return;

    su2_link U = make_unitary(value);
    if (!U.allFinite()) {
        cout << "Non-finite matrix element found" << endl;
        return;
    }

    link[d] = U;
    link_inverse[d] = U.adjoint();
}

su2_link su2_x_site::create_random_link() {

    double a = rand_double(-1.0, 1.0);
    double g = sqrt(1 - a * a);
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

su2_link su2_x_site::create_link(double a) {

    double g = sqrt(1 - a * a);
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

double su2_x_site::action() {
    // gauge action
    double S = 0.0;
    su2_link s;
    double n = 0.0;
    for (int d1 = 0; d1 < lattice->D - 1; d1++) {
        for (int d2 = d1 + 1; d2 < lattice->D; d2++) {
            s = link[d1];
            s *= forward[d1]->link[d2];
            s *= forward[d2]->link_inverse[d1];
            s *= link_inverse[d2];
            S += (su2_identity - s).trace().real() / 2.0;
            n += 1.0;
        }
    }
    return eps * lattice->beta * S / n;
}

double su2_x_site::hamiltonian() {
    // wilson action (including extra dimension)
    double S = 0.0;
    su2_link s;
    int d_max = forward_edge ? lattice->D : (lattice->D + 1);
    for (int d1 = 0; d1 < d_max - 1; d1++) {
        for (int d2 = d1 + 1; d2 < d_max; d2++) {
            s = link[d1];
            s *= forward[d1]->link[d2];
            s *= forward[d2]->link_inverse[d1];
            s *= link_inverse[d2];
            if (d2 == lattice->D) {
                S += (su2_identity - s).trace().real() * lattice->eps5;
            } else {
                S += (su2_identity - s).trace().real() * eps;
            }
        }
    }
    S = lattice->beta * S / 2.0;

    // extra dimension never contributes to kinetic energy
    for (int d = 0; d < lattice->D; d++) S += (p_link[d] * p_link[d]).trace().real();
    return S;
}

void su2_x_site::init_momenta() {
    for (int d = 0; d < lattice->D; d++) {
        p_link[d] = su2_zero;
        for (int i = 1; i < 4; i++) {
            p_link[d] += rand_normal() * lattice->sigma[i] * 0.5;
        }
    }
    // momentum in extra dimension is always zero
    p_link[lattice->D] = su2_zero;
}

void su2_x_site::hmc_step_p(double frac) {
    for (int d = 0; d < lattice->D; d++) {
        p_link[d] -= frac * lattice->dt * p_link_dot(d);
    }
}

void su2_x_site::hmc_step_link() {
    for (int d = 0; d < lattice->D; d++) {
        set_link(d, exp_su2(I * lattice->dt * p_link[d]) * link[d]);
    }
}

su2_link su2_x_site::p_link_dot(int d) {

    su2_link A = staple(d) * eps + staple_x(d) * lattice->eps5;
    su2_link U = link[d];
    su2_link UA = U * A;
    su2_link AU_dag = A.adjoint() * U.adjoint();
    su2_link X = (UA - AU_dag) - (UA - AU_dag).trace() * su2_identity * 0.5;
    return -I * lattice->beta * X * 0.125;
}

double su2_x_site::link_trace() {
    // compute the average link trace at this site
    double Tr = 0;
    for (int d1 = 0; d1 < lattice->D; d1++) {
        Tr += link[d1].trace().real();
    }
    return Tr / lattice->D / 2.0;
}

double su2_x_site::plaq() {
    // compute the average plaquette at this site
    double U = 0;
    su2_link u;
    for (int d1 = 0; d1 < lattice->D; d1++) {
        for (int d2 = d1 + 1; d2 < lattice->D; d2++) {
            u = link[d1];
            u *= forward[d1]->link[d2];
            u *= forward[d2]->link_inverse[d1];
            u *= link_inverse[d2];
            U += u.trace().real() / 2.0;
        }
    }
    return U / lattice->D / (lattice->D - 1) * 2.0;
}

su2_link su2_x_site::plaquette(int d1, int d2) {
    // compute the plaquette in the d1,d2 direction
    su2_link u1, u2, u3, u4;
    u1 = link[d1];
    u2 = forward[d1]->link[d2];
    u3 = forward[d2]->link_inverse[d1];
    u4 = link_inverse[d2];
    return u1 * u2 * u3 * u4;
}

su2_link su2_x_site::staple(int d1) {
    // compute the sum of staples attached to the plaquette pointing in the d1 direction
    su2_link U = su2_zero;
    su2_link u;
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

    return U;
}

su2_link su2_x_site::reverse_staple(int d1) {
    // compute the sum of staples attached to the plaquette pointing in the d1 direction
    su2_link U = su2_zero;
    su2_link u;
    for (int d2 = 0; d2 < lattice->D; d2++) {
        if (d1 == d2) continue;

        // forward staple
        u = link[d2];
        u *= forward[d2]->link[d1];
        u *= forward[d1]->link_inverse[d2];
        U += u;

        // backward staple
        u = backward[d2]->link_inverse[d2];
        u *= backward[d2]->link[d1];
        u *= backward[d2]->forward[d1]->link[d2];
        U += u;
    }

    return U;
}

su2_link su2_x_site::staple_x(int d1) {
    // compute the sum of extra dimension staples attached to the plaquette pointing in the d1 direction
    su2_link U = su2_zero;

    if (d1 == lattice->D) return U;

    if (!forward_edge) {
        // forward staple in extra dimension
        U += forward[lattice->D]->link_inverse[d1];
    }

    if (!backward_edge) {
        // backward staple in extra dimension
        U += backward[lattice->D]->link_inverse[d1];
    }

    return U;
}

su2_link su2_x_site::reverse_staple_x(int d1) {
    // compute the sum of extra dimension staples attached to the plaquette pointing in the d1 direction
    su2_link U = su2_zero;

    if (d1 == lattice->D) return U;

    if (!forward_edge) {
        // forward staple in extra dimension
        U += forward[lattice->D]->link[d1];
    }

    if (!backward_edge) {
        // backward staple in extra dimension
        U += backward[lattice->D]->link[d1];
    }

    return U;
}

su2_link su2_x_site::cloverleaf(int d1, int d2) {
    su2_link U = su2_identity;
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

double su2_x_site::wilson_loop(int a, int b) {
    // compute the average a x b wilson loop at this site
    double U = 0;
    su2_link u;
    su2_x_site* s = this;
    int x;

    // don't use the time direction on an asymmetric lattice
    int d_min = (lattice->N == lattice->T) ? 0 : 1;
    int n = 0; // number of measurements

    for (int d1 = d_min; d1 < lattice->D; d1++) {
        for (int d2 = d1 + 1; d2 < lattice->D; d2++) {
            u = su2_identity;
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
            U += u.trace().real() / 2.0;
            n++;
        }
    }
    return U / double(n);
}

complex<double> su2_x_site::polyakov_loop(int d) {
    // compute the polyakov loops at this site in the d direction
    su2_x_site* s = this; // current site
    su2_link u1 = su2_identity;
    int t_max = (d == 0) ? lattice->T : lattice->N;

    for (int t = 0; t < t_max; t++) {
        u1 *= s->link[d];
        s = s->forward[d];
    }

    return u1.trace();
}

double su2_x_site::correlator(int T) {
    // correlator of duration T
    su2_link u = su2_identity;
    su2_x_site* s = this;
    int x;
    for (x = 0; x < T; x++) {
        u *= s->link[0];
        s = s->forward[0];  // time direction is always 0
    }
    return u.trace().real() / 2.0;
}

double su2_x_site::field_strength() {

    //            / +-->--+    +-->--+                          \
    //           /  |     |    |     |         u3           u4   \
    //           |  |     |    |     |                           |
    // U = 1/4 * | (x)-<--+ +  +--<-(x)  +  +-->-(x)  + (x)->--+ |
    //           |                          |     |      |     | |
    //           \     u1         u2        |     |      |     | /
    //            \                         +--<--+      +--<--+/
    //
    // F = (U - U^dag) / 2
    // X = F - Tr(F) / 2
    // E = -X * X / 4

    double E = 0.0;
    su2_link u1, u2, u3, u4;

    for (int d1 = 0; d1 < lattice->D; d1++) {
        for (int d2 = 0; d2 < lattice->D; d2++) {
            if (d1 == d2) continue;
            u1 = link[d2];
            u1 *= forward[d2]->link[d1];
            u1 *= forward[d1]->link_inverse[d2];
            u1 *= link_inverse[d1];

            u2 = backward[d1]->link_inverse[d1];
            u2 *= backward[d1]->link[d2];
            u2 *= backward[d1]->forward[d2]->link[d1];
            u2 *= link_inverse[d2];

            u3 = backward[d2]->link_inverse[d2];
            u3 *= backward[d2]->backward[d1]->link_inverse[d1];
            u3 *= backward[d2]->backward[d1]->link[d2];
            u3 *= backward[d1]->link[d1];

            u4 = link[d1];
            u4 *= backward[d2]->forward[d1]->link_inverse[d2];
            u4 *= backward[d2]->link_inverse[d1];
            u4 *= backward[d2]->link[d2];

            // this matches how grid does it, gives the same result
//            u1 = link[d1];
//            u1 *= forward[d1]->link[d2];
//            u1 *= forward[d2]->link_inverse[d1];
//            u1 *= link_inverse[d2];
//
//            u2 = link[d1];
//            u2 *= backward[d2]->forward[d1]->link_inverse[d2];
//            u2 *= backward[d2]->link_inverse[d1];
//            u2 *= backward[d2]->link[d2];
//
//            u3 = backward[d2]->link_inverse[d2];
//            u3 *= backward[d1]->backward[d2]->link_inverse[d1];
//            u3 *= backward[d1]->backward[d2]->link[d2];
//            u3 *= backward[d1]->link[d1];
//
//            u4 = link[d2];
//            u4 *= backward[d1]->forward[d2]->link_inverse[d1];
//            u4 *= backward[d1]->link_inverse[d2];
//            u4 *= backward[d1]->link[d1];

            su2_link U = (u1 + u2 + u3 + u4) / 4.0; // average of 4 plaquettes
            su2_link F = (U - U.adjoint()) / 2.0; // make anti-hermitian
            su2_link X = F - F.trace() * su2_identity / 2.0; // make traceless
            E -= real((X * X).trace()) / 4.0; // add field strength
        }
    }

    return E;
}

//double su2_x_site::field_strength_x() {
//
//    // same as field_strength(), but includes links in the extra dimension
//
//    double E = 0.0;
//    su2_link u1, u2, u3, u4;
//
//    for (int d1 = 0; d1 < lattice->D; d1++) {
//        for (int d2 = 0; d2 < lattice->D; d2++) {
//            u1 = link[d1];
//            u1 *= forward[d1]->link[d2];
//            u1 *= forward[d2]->link_inverse[d1];
//            u1 *= link_inverse[d2];
//
//            u2 = backward[d2]->link_inverse[d2];
//            u2 *= backward[d2]->link[d1];
//            u2 *= backward[d2]->forward[d1]->link[d2];
//            u2 *= link_inverse[d1];
//
//            u3 = backward[d1]->link_inverse[d1];
//            u3 *= backward[d1]->backward[d2]->link_inverse[d2];
//            u3 *= backward[d1]->backward[d2]->link[d1];
//            u3 *= backward[d2]->link[d2];
//
//            u4 = link[d2];
//            u4 *= backward[d1]->forward[d2]->link_inverse[d1];
//            u4 *= backward[d1]->link_inverse[d2];
//            u4 *= backward[d1]->link[d1];
//
//            su2_link U = (u1 + u2 + u3 + u4) / 4.0;
//            su2_link F = (U - U.adjoint()) / 2.0;
//            E += -real((F * F).trace()) / 4.0;
//        }
//
//        u1 = su2_zero;
//        u2 = su2_zero;
//        u3 = su2_zero;
//        u4 = su2_zero;
//
//        if (!forward_edge) {
//            u1 = link[d1] * forward[lattice->D]->link_inverse[d1];
//            u2 = link_inverse[d1] * forward[lattice->D]->link[d1];
//        }
//
//        if (!backward_edge) {
//            u3 = link[d1] * forward[lattice->D]->link_inverse[d1];
//            u4 = link_inverse[d1] * forward[lattice->D]->link[d1];
//        }
//
//        su2_link U = (u1 + u2 + u3 + u4) * lattice->eps5 / 4.0;
//        su2_link F = (U - U.adjoint()) / 2.0;
//        E += -real((U * U).trace()) / 4.0;
//    }
//
//    return E;
//}

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

double su2_x_site::topological_charge() {

    // symmetric plaquettes
    su2_link u1, u2, u3, u4;
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
        su2_link Up1p2 = u1 * u2 * u3 * u4;

        u1 = link[d1];
        u2 = forward[d1]->backward[d2]->link_inverse[d2];
        u3 = backward[d2]->link_inverse[d1];
        u4 = backward[d2]->link[d2];
        su2_link Up1n2 = u1 * u2 * u3 * u4;

        u1 = backward[d1]->link_inverse[d1];
        u2 = backward[d1]->link[d2];
        u3 = backward[d1]->forward[d2]->link[d1];
        u4 = link_inverse[d2];
        su2_link Un1p2 = u1 * u2 * u3 * u4;

        u1 = backward[d1]->link_inverse[d1];
        u2 = backward[d1]->backward[d2]->link_inverse[d2];
        u3 = backward[d1]->backward[d2]->link[d1];
        u4 = backward[d2]->link[d2];
        su2_link Un1n2 = u1 * u2 * u3 * u4;

        u1 = link[d3];
        u2 = forward[d3]->link[d4];
        u3 = forward[d4]->link_inverse[d3];
        u4 = link_inverse[d4];
        su2_link Up3p4 = u1 * u2 * u3 * u4;

        u1 = link[d3];
        u2 = forward[d3]->backward[d4]->link_inverse[d4];
        u3 = backward[d4]->link_inverse[d3];
        u4 = backward[d4]->link[d4];
        su2_link Up3n4 = u1 * u2 * u3 * u4;

        u1 = backward[d3]->link_inverse[d3];
        u2 = backward[d3]->link[d4];
        u3 = backward[d3]->forward[d4]->link[d3];
        u4 = link_inverse[d4];
        su2_link Un3p4 = u1 * u2 * u3 * u4;

        u1 = backward[d3]->link_inverse[d3];
        u2 = backward[d3]->backward[d4]->link_inverse[d4];
        u3 = backward[d3]->backward[d4]->link[d3];
        u4 = backward[d4]->link[d4];
        su2_link Un3n4 = u1 * u2 * u3 * u4;

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
//    su2_link C01 = cloverleaf(0, 1).imag();
//    su2_link C02 = cloverleaf(0, 2).imag();
//    su2_link C03 = cloverleaf(0, 3).imag();
//    su2_link C10 = cloverleaf(1, 0).imag();
//    su2_link C12 = cloverleaf(1, 2).imag();
//    su2_link C13 = cloverleaf(1, 3).imag();
//    su2_link C20 = cloverleaf(2, 0).imag();
//    su2_link C21 = cloverleaf(2, 1).imag();
//    su2_link C23 = cloverleaf(2, 3).imag();
//    su2_link C30 = cloverleaf(3, 0).imag();
//    su2_link C31 = cloverleaf(3, 1).imag();
//    su2_link C32 = cloverleaf(3, 2).imag();
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

double su2_x_site::mag_U() {
    return this->link[0].trace().real() / 2.0;
}

double su2_x_site::abs_U() {
    return abs(this->link[0].trace().real()) / 2.0;
}

double su2_x_site::four_point(int T, int R) {
    // four-point function of duration T and spacing R
    double U = 0;
    su2_link u;
    su2_x_site* s = this;
    int x;
    for (int d = 1; d < lattice->D; d++) { // d is the spacing direction
        u = su2_identity;
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
        U += u.trace().real() / 2.0;
    }
    return U / (lattice->D - 1);
}

su2_link su2_x_site::overrelax(su2_link g) {

    double w = 1.7;
    int n_max = 5;

    su2_link g1 = g - su2_identity; // original relaxation matrix minus identity
    su2_link gn = su2_identity; // original g raised to the nth power
    su2_link gw = su2_zero; // overrelaxation matrix

    for (int n = 0; n <= n_max; n++) {
        double wn = tgamma(w + 1) / tgamma(w + 1 - n) / tgamma(n + 1);
        gw += wn * gn;
        gn *= g1;
    }
    return make_unitary(gw);
}

void su2_x_site::relax(bool coulomb) {
    while (!lock()) {}
    su2_link g = make_unitary(su2_identity - sum_G(coulomb) * I / 4);
   g = overrelax(g);

    // apply the gauge transformation (including the time direction)
    for (int d = 0; d < lattice->D; d++) {
        set_link(d, g * link[d]);
        backward[d]->set_link(d, backward[d]->link[d] * g.adjoint());
    }
    unlock();
}

double su2_x_site::sum_landau() {
    double U = 0.0;
    for (int d = 0; d < lattice->D; d++) U += link[d].trace().real();
    return U / 2.0 / lattice->D;
}

double su2_x_site::sum_coulomb() {
    double U = 0.0;
    for (int d = 1; d < lattice->D; d++) U += link[d].trace().real();
    return U / 2.0 / (lattice->D - 1);
}

su2_link su2_x_site::sum_G(bool coulomb) {
    su2_link G = su2_zero;
    for (int d = coulomb ? 1 : 0; d < lattice->D; d++) {
        G += (link[d] - link_inverse[d]) / 2.0 / I;
        G -= (backward[d]->link[d] - backward[d]->link_inverse[d]) / 2.0 / I;
    }
    return G;
}

void su2_x_site::heat_bath() {
    while (lock());
    for (int d = 0; d < lattice->D; d++) heat_bath_link(d);
    unlock();
}

void su2_x_site::heat_bath_link(int d1) {
    su2_link old_link = link[d1];

    su2_link u;
    su2_link A = su2_zero;
    for (int d2 = 0; d2 < lattice->D; d2++) {
        if (d1 == d2) continue;

        // forward staple
        u = forward[d1]->link[d2];
        u *= forward[d2]->link_inverse[d1];
        u *= link_inverse[d2];
        A += u * eps;

        // backward staple
        u = backward[d2]->forward[d1]->link_inverse[d2];
        u *= backward[d2]->link_inverse[d1];
        u *= backward[d2]->link[d2];
        A += u * eps;
    }

    if (!forward_edge) {
        u = forward[d1]->link[lattice->D];
        u *= forward[lattice->D]->link_inverse[d1];
        u *= link_inverse[lattice->D];
        A += u * lattice->eps5;
    }

    if (!backward_edge) {
        u = backward[lattice->D]->forward[d1]->link_inverse[lattice->D];
        u *= backward[lattice->D]->link_inverse[d1];
        u *= backward[lattice->D]->link[lattice->D];
        A += u * lattice->eps5;
    }

    // heat bath algorithm
    double k = sqrt(A.determinant().real());
    double x, a, g;

    while (true) {
        x = rand_double(exp(-2 * k * lattice->beta), 1.0);
        a = 1 + log(x) / k / lattice->beta;
        g = sqrt(1 - a * a);
        if (rand_double() < g) break;
    }
    su2_link new_link = create_link(a) * (A / k).inverse();
    set_link(d1, new_link);
}

void su2_x_site::cool() {
    lock();
    for (int d = 0; d < lattice->D; d++) cool_link(d);
    unlock();
}

void su2_x_site::cool_link(int d1) {

    su2_link X = su2_zero;
    su2_link u1, u2, u3;
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

//void su2_x_site::wilson_flow(su2_x_site* target, double dtau) {
//    lock();
//    for (int d = 0; d < lattice->D; d++) wilson_flow_link(target, dtau, d);
//    unlock();
//}
//
//void su2_x_site::wilson_flow_link(su2_x_site* target, double dtau, int d1) {
//
//    su2_link X;
//    su2_link u1, u2, u3;
//    for (int d2 = 0; d2 < lattice->D; d2++) {
//        if (d1 == d2) continue;
//
//        // forward staple
//        u1 = link[d2];
//        u2 = forward[d2]->link[d1];
//        u3 = forward[d1]->link_inverse[d2];
//        X += u1 * u2 * u3;
//
//        // backward staple
//        u1 = backward[d2]->link_inverse[d2];
//        u2 = backward[d2]->link[d1];
//        u3 = forward[d1]->backward[d2]->link[d2];
//        X += u1 * u2 * u3;
//    }
//    X = X.adjoint().eval();
//
//    su2_link W0 = link[d1];
//    su2_link Z0 = W0 * X;
//    Z0 = (Z0.adjoint() - Z0).eval();
//    Z0 = 0.5 * Z0 - Z0.trace() * su2_identity / 6.0;
//
//    su2_link W1 = exp_su2(-I * 0.25 * Z0 * dtau) * W0;
//    su2_link Z1 = W1 * X;
//    Z1 = (Z1.adjoint() - Z1).eval();
//    Z1 = 0.5 * Z1 - Z1.trace() * su2_identity / 6.0;
//    Z1 = (8.0/9.0) * Z1 - (17.0/36.0) * Z0;
//
//    su2_link W2 = exp_su2(-I * Z1 * dtau) * W1;
//    su2_link Z2 = W2 * X;
//    Z2 = (Z2.adjoint() - Z2).eval();
//    Z2 = 0.5 * Z2 - Z2.trace() * su2_identity / 6.0;
//    Z2 = 0.75 * Z2 - Z1;
//
//    su2_link W3 = exp_su2(-I * Z2 * dtau) * W2;
//    target->set_link(d1, W3);
//}

void su2_x_site::wilson_flow(su2_x_site* target, int step) {
    target->lock();
    for (int d = 0; d < lattice->D; d++) wilson_flow_link(target, step, d);
    target->unlock();
}

void su2_x_site::wilson_flow_link(su2_x_site* target, int step, int d) {

//    // get z from reverse staple, including extra dimension
//    su2_link x = reverse_staple(d) + reverse_staple_x(d) * lattice->eps5;

    // get z from reverse staple, excluding extra dimension
    su2_link x = reverse_staple(d);
    x.adjointInPlace();
    su2_link w = link[d];
    su2_link z = w * x;

    // get traceless, antihermitian part
    z = (z - z.adjoint().eval()).eval() / 2.0;
    z = z - z.trace() * su2_identity / 2.0;

    if (step == 1) {
        z = z / 4.0;
    } else if (step == 2) {
        z = (8.0/9.0) * z - (17.0/9.0) * wf_z[d];
    } else if (step == 3) {
        z = (3.0/4.0) * z - wf_z[d];
    }
    target->wf_z[d] = z;
    target->set_link(d, exp_su2(I * lattice->wf_dt * z) * w);
}

void su2_x_site::stout_smear(su2_x_site* target, double rho) {
    for (int d = 0; d < lattice->D; d++) stout_smear_link(target, rho, d);
}

void su2_x_site::stout_smear_link(su2_x_site* target, double rho, int d1) {
    su2_link C = su2_zero;
    su2_link u1, u2, u3;
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

    su2_link X = rho * C * link_inverse[d1];
    su2_link X_dag = X.adjoint();
    su2_link Z = X_dag - X;
    su2_link W = 0.5 * I * (Z - su2_identity * Z.trace() * 0.5);
    target->set_link(d1, exp_su2(W) * link[d1]);
}

// MARK: su2_x_lattice

su2_x_lattice::su2_x_lattice(int N, int T, int N5, int D, double beta, double eps5, bool cold) {
    this->N = N;
    this->T = T;
    this->N5 = N5;
    this->D = D;
    this->beta = beta;
    this->eps5 = eps5;
    this->parallel = false;
    this->n_steps = 30;
    this->dt = 0.03333333;
    this->verbose = 0;

    init();

    for (int s = 0; s < n_sites; s++) site[s].reset_links(cold);
}

su2_x_lattice::su2_x_lattice(su2_x_lattice* lattice) {
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

su2_x_lattice::su2_x_lattice(int N, int T, int N5, int D, double beta, double eps5, ifstream& ckptFile, bool isNersc) {

    this->N = N;
    this->T = T;
    this->N5 = N5;
    this->D = D;
    this->beta = beta;
    this->eps5 = eps5;
    this->parallel = false;
    this->n_steps = 30;
    this->dt = 0.03333333;
    this->verbose = 0;

    // read lattice data from a Grid checkpoint file
    string line;
    bool isHdr = true;
    double hdrLinkTrace = 0.0;
    double hdrPlaq = 0.0;
    int checksum = 0;
    bool bigEndian = false;

    ckptFile.seekg(0);

    // initialize the lattice
    init();

    for (int s = 0; s < n_sites; s++) site[s].read_links(ckptFile, bigEndian);

    // check average plaquette
    printf("Link Trace (Calculated) = %lf\n", link_trace());
    for (int n5 = 0; n5 < N5; n5++) {
        printf("Plaquette(%d) = %lf\n", n5, plaq(n5) / 2.0);
    }
}

su2_x_lattice::~su2_x_lattice() {
}

void su2_x_lattice::init() {

    // initialize hmc and heat bath parameters
    hmc_count = 0.0;
    hmc_accept = 0.0;

    // wilson flow time step
    wf_dt = 0.02;

    // initialize pauli matrices
    sigma.resize(4);
    sigma[0] << 1, 0, 0, 1;
    sigma[1] << 0, 1, 1, 0;
    sigma[2] << 0, -I, I, 0;
    sigma[3] << 1, 0, 0, -1;

    // initialize rng
    random_device rd;  // seed for the random number engine
    for (int t = 0; t < (T + 1); t++) {
        // mersenne_twister_engine seeded with rd()}
        gen.push_back(mt19937(rd()));
    }

    // initialize sites
    n_sites = pow(N, D - 1) * T * N5;
    n_sites_5 = n_sites / N5;
    n_slice = n_sites / T;
    n5_center = (N5 - 1) / 2;
    site.resize(n_sites);
    site_1.resize(n_sites);

    for (int s = 0; s < n_sites; s++) {
        site[s].init(this, site.data(), s);
        site_1[s].init(this, site_1.data(), s);
    }
}

double su2_x_lattice::rand_double(double min, double max) {
    std::uniform_real_distribution<double> dist(min, max);
    return dist(gen[T]);
}

int su2_x_lattice::rand_int(int min, int max) {
    std::uniform_int_distribution<int> dist(min, max);
    return dist(gen[T]);
}

double su2_x_lattice::rand_normal(double mean, double stdev) {
    std::normal_distribution<double> dist(mean, stdev);
    return dist(gen[T]);
}

int su2_x_lattice::get_site(int s, int d, int n) {

    // get the site's coordinates
    int x[D + 1];
    x[0] = s / n_slice;
    for (int d1 = 1; d1 < D; d1++) x[d1] = (s / (int)pow(N, D - d1 - 1) / N5) % N;
    x[D] = s % N5;

    // go n steps in the d direction
    if (d == 0) {
        // time direction
        x[0] += n;
        while (x[0] < 0) x[0] += T;
        x[0] %= T;
    } else if (d == D) {
        // extra direction
        x[D] += n;
        while (x[D] < 0) x[D] += N5;
        x[D] %= N5;
    } else {
        // spatial direction
        x[d] += n;
        while (x[d] < 0) x[d] += N;
        x[d] %= N;
    }

    // get the index of the new site
    s = 0;
    for (int d1 = 0; d1 < D; d1++) s += x[d1] * (int)pow(N, D - d1 - 1) * N5;
    s += x[D];
    return s;
}

double su2_x_lattice::link_trace() {
    // compute the average link trace of the lattice (normalized to 1)
    double P = 0;
    for (int s = 0; s < n_sites; s++) P += site[s].link_trace();
    return P / double(n_sites);
}

double su2_x_lattice::plaq(int n5) {
    // compute the average plaquette of the lattice
    double P = 0;
    for (int s = n5; s < n_sites; s += N5) P += site[s].plaq();
    return P / double(n_sites_5);
}

double su2_x_lattice::action(int n5) {
    double S = 0;
    for (int s = n5; s < n_sites; s += N5) S += site[s].action();
    return S / double(n_sites_5);
}

double async_wilson_loop(su2_x_site* site, int n, int a, int b, int N5) {
    double sum = 0.0;
    for (int s = 0; s < n; s += N5) sum += site[s].wilson_loop(a, b);
    return sum;
}

double su2_x_lattice::wilson_loop(int a, int b, int n5) {
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

double su2_x_lattice::polyakov_loop(int R, int n5) {
    // compute the average product of polyakov loops spaced r apart

    // on a square lattice, any direction can be the time direction
    int d_max = (T == N) ? (D - 1) : 1;
    int n = 0;

    complex<double> p[N];
    complex<double> P, p1, p2;
    su2_x_site* s;
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

double su2_x_lattice::correlator(int T, int n5) {
    // compute the average correlator of length T
    double sum = 0;
    for (int s = n5; s < n_sites; s += N5) sum += site[s].correlator(T);
    return sum / double(n_sites_5);
}

double su2_x_lattice::four_point(int T, int R, int n5) {
    // compute the average T x R 4-point function of the lattice
    double sum = 0;
    for (int s = n5; s < n_sites; s += N5) sum += site[s].four_point(T, R);
    return sum / double(n_sites_5);
}

double su2_x_lattice::hamiltonian() {
    double H = 0;
    for (int s = 0; s < n_sites; s++) H += site_1[s].hamiltonian();
    return H;
}

void async_init_momenta(su2_x_site* site, int n) {
    for (int s = 0; s < n; s++) site[s].init_momenta();
}

void async_hmc_step_p(su2_x_site* site, int n, double frac) {
    for (int s = 0; s < n; s++) site[s].hmc_step_p(frac);
}

void async_hmc_step_link(su2_x_site* site, int n) {
    for (int s = 0; s < n; s++) site[s].hmc_step_link();
}

void su2_x_lattice::hmc(int n_sweeps, bool update_dt, bool no_metropolis) {

    if (n_sweeps) {
        for (int i = 0; i < n_sweeps; i++) hmc(0, update_dt, no_metropolis);
        return;
    }

    // copy lattice sites
    for (int s = 0; s < n_sites; s++) site_1[s].copy_links(&site[s]);

    // init conjugate momenta
    if (parallel) {
        future<void> async_slice[T];
        for (int t = 0; t < T; t++) async_slice[t] = async(async_init_momenta, &site_1[t * n_slice], n_slice);
        for (int t = 0; t < T; t++) async_slice[t].get();
    } else {
        for (int s = 0; s < n_sites; s++) site_1[s].init_momenta();
    }

    // calculate old hamiltonian
    double oldH = hamiltonian();

    // do the initial half-step
    if (parallel) {
        future<void> async_slice[T];
        for (int t = 0; t < T; t++) async_slice[t] = async(async_hmc_step_p, &site_1[t * n_slice], n_slice, 0.5);
        for (int t = 0; t < T; t++) async_slice[t].get();
    } else {
        for (int s = 0; s < n_sites; s++) site_1[s].hmc_step_p(0.5);
    }

    // do leapfrog steps
    for (int step = 0; step < (n_steps - 1); step++) {
        if (parallel) {
            future<void> async_slice[T];
            for (int t = 0; t < T; t++) async_slice[t] = async(async_hmc_step_link, &site_1[t * n_slice], n_slice);
            for (int t = 0; t < T; t++) async_slice[t].get();
            for (int t = 0; t < T; t++) async_slice[t] = async(async_hmc_step_p, &site_1[t * n_slice], n_slice, 1.0);
            for (int t = 0; t < T; t++) async_slice[t].get();
        } else {
            for (int s = 0; s < n_sites; s++) site_1[s].hmc_step_link();
            for (int s = 0; s < n_sites; s++) site_1[s].hmc_step_p(1.0);
        }
    }

    // do the final half-step
    if (parallel) {
        future<void> async_slice[T];
        for (int t = 0; t < T; t++) async_slice[t] = async(async_hmc_step_link, &site_1[t * n_slice], n_slice);
        for (int t = 0; t < T; t++) async_slice[t].get();
        for (int t = 0; t < T; t++) async_slice[t] = async(async_hmc_step_p, &site_1[t * n_slice], n_slice, 0.5);
        for (int t = 0; t < T; t++) async_slice[t].get();
    } else {
        for (int s = 0; s < n_sites; s++) site_1[s].hmc_step_link();
        for (int s = 0; s < n_sites; s++) site_1[s].hmc_step_p(0.5);
    }

    // calculate new hamiltonian
    double newH = hamiltonian();
    double dH = newH - oldH;
    // printf("oldH: %.6f, newH: %.6f, dH: %.6f\n", oldH, newH, dH);

    // accept with metropolis algorithm
    bool accept = false;
    if (dH <= 0 || rand_double() < exp(-dH)) {
        accept = true;
        hmc_accept++;
    }

    if (accept || no_metropolis) {
        for (int s = 0; s < n_sites; s++) site[s].copy_links(&site_1[s]);
    }

    hmc_count++;
//    if (hmc_count >= 20) {
//        // update hmc statistics every 20 iterations
//        double accept_rate = double(hmc_accept) / double(hmc_count);
//        hmc_count = 0;
//        hmc_accept = 0;
//
//        // target a 70%-90% acceptance rate
//        if (update_dt) {
//            if (accept_rate > 0.91) {
//                dt *= 1.1;
//            } else if (accept_rate < 0.35) {
//                dt *= 0.5;
//            } else if (accept_rate < 0.69) {
//                dt *= 0.9;
//            }
//        }
//        if (verbose) {
//            cout << setprecision(6) << fixed;
//            cout << "accept_rate: " << accept_rate;
//            cout << ", dt: " << dt;
//            cout << ", plaq: " << plaq(n5_center) << endl;
//        }
//    }
}

void async_heat_bath(su2_x_site* site, int n) {

    for (int s = 0; s < n; s++) site[s].heat_bath();
}

void su2_x_lattice::heat_bath(int n_sweeps) {

    if (n_sweeps) {
        for (int i = 0; i < n_sweeps; i++) heat_bath(0);
        return;
    }

    if (parallel) {
        future<void> async_slice[T];
        for (int t = 0; t < T; t += 2) async_slice[t] = async(async_heat_bath, &site[t * n_slice], n_slice);
        for (int t = 0; t < T; t += 2) async_slice[t].get();
        for (int t = 1; t < T; t += 2) async_slice[t] = async(async_heat_bath, &site[t * n_slice], n_slice);
        for (int t = 1; t < T; t += 2) async_slice[t].get();
    } else {
        for (int s = 0; s < n_sites; s++) site[s].heat_bath();
    }

    if (verbose) {
        cout << setprecision(6);
        cout << "plaq: " << plaq(n5_center) << endl;
    }
}

void async_cool(su2_x_site* site, int n, int N5) {
    for (int s = 0; s < n; s += N5) site[s].cool();
}

void su2_x_lattice::cool(int n5, int n_sweeps) {

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

//void async_wilson_flow(su2_x_site* site, su2_x_site* target, int n, int step) {
//    for (int s = 0; s < n; s++) site[s].wilson_flow(&target[s], step);
//}
//
//void su2_x_lattice::wilson_flow(int n_sweeps) {
//    if (n_sweeps) {
//        for (int i = 0; i < n_sweeps; i++) wilson_flow();
//        return;
//    }
//
//    if (parallel) {
//        future<void> async_slice[T];
//        for (int t = 0; t < T; t++) async_slice[t] = async(async_wilson_flow, &site[t * n_slice], &site_1[t * n_slice], n_slice, 1);
//        for (int t = 0; t < T; t++) async_slice[t].get();
//        for (int t = 0; t < T; t++) async_slice[t] = async(async_wilson_flow, &site_1[t * n_slice], &site[t * n_slice], n_slice, 2);
//        for (int t = 0; t < T; t++) async_slice[t].get();
//        for (int t = 0; t < T; t++) async_slice[t] = async(async_wilson_flow, &site[t * n_slice], &site_1[t * n_slice], n_slice, 3);
//        for (int t = 0; t < T; t++) async_slice[t].get();
//    } else {
//        for (int s = 0; s < n_sites; s++) site[s].wilson_flow(&site_1[s], 1);
//        for (int s = 0; s < n_sites; s++) site_1[s].wilson_flow(&site[s], 2);
//        for (int s = 0; s < n_sites; s++) site[s].wilson_flow(&site_1[s], 3);
//    }
//
//    // copy lattice sites
//    for (int s = 0; s < n_sites; s++) site[s].copy_links(&site_1[s]);
//}

void async_wilson_flow(su2_x_site* site, su2_x_site* target, int n, int N5, int step) {
    for (int s = 0; s < n; s += N5) site[s].wilson_flow(&target[s], step);
}

void su2_x_lattice::wilson_flow(int n5, int n_sweeps) {
    if (n_sweeps) {
        for (int i = 0; i < n_sweeps; i++) wilson_flow(n5, 0);
        return;
    }

    if (parallel) {
        future<void> async_slice[T];
        for (int t = 0; t < T; t++) async_slice[t] = async(async_wilson_flow, &site[t * n_slice + n5], &site_1[t * n_slice + n5], n_slice, N5, 1);
        for (int t = 0; t < T; t++) async_slice[t].get();
        for (int t = 0; t < T; t++) async_slice[t] = async(async_wilson_flow, &site_1[t * n_slice + n5], &site[t * n_slice + n5], n_slice, N5, 2);
        for (int t = 0; t < T; t++) async_slice[t].get();
        for (int t = 0; t < T; t++) async_slice[t] = async(async_wilson_flow, &site[t * n_slice + n5], &site_1[t * n_slice + n5], n_slice, N5, 3);
        for (int t = 0; t < T; t++) async_slice[t].get();
    } else {
        for (int s = n5; s < n_sites; s += N5) site[s].wilson_flow(&site_1[s], 1);
        for (int s = n5; s < n_sites; s += N5) site_1[s].wilson_flow(&site[s], 2);
        for (int s = n5; s < n_sites; s += N5) site[s].wilson_flow(&site_1[s], 3);
    }

    // copy lattice sites
    for (int s = 0; s < n_sites; s++) site[s].copy_links(&site_1[s]);
}

void async_stout_smear(su2_x_site* site, su2_x_site* target, int n, int N5, double rho) {
    for (int s = 0; s < n; s += N5) site[s].stout_smear(&target[s], rho);
}

void su2_x_lattice::stout_smear(int n5, double rho, int n_sweeps) {
    if (n_sweeps) {
        for (int i = 0; i < n_sweeps; i++) stout_smear(n5, rho, 0);
        return;
    }

    if (parallel) {
        future<void> async_slice[T];
        for (int t = 0; t < T; t++) async_slice[t] = async(async_stout_smear, &site[t * n_slice + n5], &site_1[t * n_slice + n5], n_slice, N5, rho);
        for (int t = 0; t < T; t++) async_slice[t].get();
    } else {
        for (int s = n5; s < n_sites; s += N5) site[s].stout_smear(&site_1[s], rho);
    }

    // copy lattice sites
    for (int s = 0; s < n_sites; s++) site[s].copy_links(&site_1[s]);
}

double async_field_strength(su2_x_site* site, int n, int N5) {
    double E = 0.0;
    for (int s = 0; s < n; s += N5) E += site[s].field_strength();
    return E;
}

double su2_x_lattice::field_strength(int n5) {
    double E = 0.0;

    if (parallel) {
        future<double> sweep_slice[T];
        for (int t = 0; t < T; t++) sweep_slice[t] = async(async_field_strength, &site[t * n_slice + n5], n_slice, N5);
        for (int t = 0; t < T; t++) E += sweep_slice[t].get();

    } else {
        for (int s = n5; s < n_sites; s += N5) E += site[s].field_strength();
    }
    return E / double(n_sites_5);
}

//double async_field_strength_x(su2_x_site* site, int n) {
//    double E = 0.0;
//    for (int s = 0; s < n; s++) E += site[s].field_strength_x();
//    return E;
//}
//
//double su2_x_lattice::field_strength_x(int n5) {
//    double E = 0.0;
//
//    if (parallel) {
//        future<double> sweep_slice[T];
//        for (int t = 0; t < T; t++) sweep_slice[t] = async(async_field_strength_x, &site[t * n_slice], n_slice);
//        for (int t = 0; t < T; t++) E += sweep_slice[t].get();
//
//    } else {
//        for (int s = n5; s < n_sites; s++) E += site[s].field_strength_x();
//    }
//    return E / double(n_sites);
//}

//double async_field_strength_x(su2_x_site* site, int n, int N5) {
//    double E = 0.0;
//    for (int s = 0; s < n; s += N5) E += site[s].field_strength_x();
//    return E;
//}
//
//double su2_x_lattice::field_strength_x(int n5) {
//    double E = 0.0;
//
//    if (parallel) {
//        future<double> sweep_slice[T];
//        for (int t = 0; t < T; t++) sweep_slice[t] = async(async_field_strength_x, &site[t * n_slice + n5], n_slice, N5);
//        for (int t = 0; t < T; t++) E += sweep_slice[t].get();
//
//    } else {
//        for (int s = n5; s < n_sites; s += N5) E += site[s].field_strength_x();
//    }
//    return E / double(n_sites_5);
//}

double async_topological_charge(su2_x_site* site, int n, int N5) {
    double Q = 0.0;
    for (int s = 0; s < n; s += N5) Q += site[s].topological_charge();
    return Q;
}

double su2_x_lattice::topological_charge(int n5) {
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

int su2_x_lattice::thermalize(int n_min, int n_max) {

    // create a cold start model
    su2_x_lattice coldStart = su2_x_lattice(N, T, N5, D, beta, eps5, true);
    coldStart.parallel = parallel;
    coldStart.verbose = verbose;
    double hot_plaquette = plaq(n5_center);
    double cold_plaquette = coldStart.plaq(n5_center);
    if (verbose) {
        cout << "P_hot: ";
        for (int n = 0; n < N5; n++) cout << plaq(n) << " ";
        cout << "P_cold: " << cold_plaquette << endl;
        cout << "S_hot: ";
        for (int n = 0; n < N5; n++) cout << action(n) << " ";
        cout << "S_cold: " << coldStart.action(n5_center) << endl;
        cout << "E_hot: ";
        for (int n = 0; n < N5; n++) cout << field_strength(n) << " ";
        cout << "E_cold: " << coldStart.field_strength(n5_center) << endl;
    }
    int m;
    for (m = 0; hot_plaquette < cold_plaquette; m++) {
        // sweep both models until the average plaquette reaches the same value
        heat_bath();
        coldStart.heat_bath();
        // hmc(0, true, true);
        // coldStart.hmc(0, true, true);
        hot_plaquette = plaq(n5_center);
        cold_plaquette = coldStart.plaq(n5_center);
        if (verbose) {
            cout << "dP: " << abs(cold_plaquette - hot_plaquette) << " P_hot: ";
            for (int n = 0; n < N5; n++) cout << plaq(n) << " ";
            cout << "P_cold: " << cold_plaquette << endl;
        }
//        if (n_max && (m > n_max)) break;
    }

    // sweep a few more times
    if (verbose) {
        cout << "hot_plaquette: " << hot_plaquette;
        cout << " cold_plaquette: " << cold_plaquette << endl;
    }
    int n_therm = max(n_min, m * 2);
    heat_bath(n_therm - m);
    // hmc(n_therm - m, true, false);

    // return value is total number of sweeps
    if (verbose) cout << "n_therm: " << n_therm << endl;
    return n_therm;
}

double su2_x_lattice::ave_link_t(int n5) {
    // average trace of links in the t direction
    // this should be ~0 in Coulomb gauge

    double U = 0;

    for (int s = n5; s < n_sites; s += N5) U += site[s].link[0].trace().real() / 2.0;

    return U / double(n_sites);
}

double su2_x_lattice::theta(bool coulomb, int n5) {
    su2_link G1, G2;
    double th = 0.0;

    for (int s = n5; s < n_sites; s += N5) {
        G1 = site[s].sum_G(coulomb);
        G1 -= (G1.trace() * su2_identity) / 2.0;
        G2 = G1.adjoint();
        th += (G1 * G2).trace().real() / 2.0;
    }

    return th;
}

double async_relax(su2_x_site* site, int n, double error_target, bool coulomb, int N5) {
    double avg_link = 0.0;
    double th = 1.0;
    while (th > error_target) {

        // start at a random spot in the time slice to mitigate parallel effects
        int start = site->rand_int(0, n - 1);
        int end = start + n;

        for (int s = start; s < end; s++) site[s % n].relax(coulomb);

        // add up all the links
        su2_link G1, G2;
        th = 0.0;
        for (int s = 0; s < n; s += N5) {
            G1 = site[s].sum_G(coulomb);
            G1 -= (G1.trace() * su2_identity) / 2.0;
            G2 = G1.adjoint();
            th += (G1 * G2).trace().real() / 2.0;
        }
        avg_link = 0.0;
        for (int s = 0; s < n; s++) avg_link += site[s].sum_coulomb();
        avg_link /= (n * N5);
    }
    return avg_link;
}

long double su2_x_lattice::relax(long double error_target, bool coulomb, int n5) {

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
                su2_link G1, G2;
                th = 0.0;
                for (int s = s1; s < s2; s++) {
                    G1 = site[s].sum_G(coulomb);
                    G1 -= (G1.trace() * su2_identity) / 2.0;
                    G2 = G1.adjoint();
                    th += (G1 * G2).trace().real() / 2.0;
                }
                avg_link_slice = 0.0;
                for (int s = s1; s < s2; s++) avg_link_slice += site[s].sum_coulomb();
                avg_link_slice /= n_slice;
//                cout << "avg_link: " << avg_link_slice << ", ";
//                cout << "theta: " << th << endl;
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
//            cout << "avg_link: " << avg_link << ", ";
//            cout << "theta: " << th << ", ";
//            cout << "plaq: " << plaq() << endl;
        }
        return avg_link;
    }
}