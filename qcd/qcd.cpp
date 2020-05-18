//
//  qcd.cpp
//  qcd_test
//
//  Created by Evan Owen on 2/28/20.
//  Copyright © 2020 Evan Owen. All rights reserved.
//

#include <iostream>
#include <random>
#include <cmath>
#include <future>
#include "Eigen/Dense"
#include "qcd.hpp"

using namespace std;

qcd_spinor::qcd_spinor() {
    for (int c = 0; c < 24; c++) component[c] = 0.0;
    component[0] = 1.0;
//    component[6] = 1.0;
//    component[12] = 1.0;
//    component[18] = 1.0;
}

qcd_spinor::qcd_spinor(qcd_spinor& psi) {
    for (int c = 0; c < 24; c++) component[c] = psi.component[c];
}

su3_vector qcd_spinor::psi(int s) {
    complex<double> c1(component[s * 6 + 0], component[s * 6 + 1]);
    complex<double> c2(component[s * 6 + 2], component[s * 6 + 3]);
    complex<double> c3(component[s * 6 + 4], component[s * 6 + 5]);
    return su3_vector(c1,c2,c3);
}

su3_vector qcd_spinor::psi_conj(int s) {
    complex<double> c1(component[s * 6 + 0], -component[s * 6 + 1]);
    complex<double> c2(component[s * 6 + 2], -component[s * 6 + 3]);
    complex<double> c3(component[s * 6 + 4], -component[s * 6 + 5]);
    return su3_vector(c1,c2,c3);
}

u1 qcd_spinor::contract(qcd_spinor psi) {
    u1 a = psi_conj(0).transpose() * psi.psi(2);
    u1 b = psi_conj(1).transpose() * psi.psi(3);
    u1 c = psi_conj(2).transpose() * psi.psi(0);
    u1 d = psi_conj(3).transpose() * psi.psi(1);
    return (a + b + c + d);
}

u1 qcd_spinor::contract(su3 U, qcd_spinor psi) {
    u1 a = psi_conj(0).transpose() * U * psi.psi(2);
    u1 b = psi_conj(1).transpose() * U * psi.psi(3);
    u1 c = psi_conj(2).transpose() * U * psi.psi(0);
    u1 d = psi_conj(3).transpose() * U * psi.psi(1);
    return (a + b + c + d);
}

u1 qcd_spinor::contract(const int gamma, su3 U, qcd_spinor psi) {
    u1 u = 0;
    if (gamma == 0) {
        u1 a = psi_conj(0).transpose() * U * psi.psi(0);
        u1 b = psi_conj(1).transpose() * U * psi.psi(1);
        u1 c = psi_conj(2).transpose() * U * psi.psi(2);
        u1 d = psi_conj(3).transpose() * U * psi.psi(3);
        return (a + b + c + d) * I;
    } else if (gamma == 1) {
        u1 a = psi_conj(0).transpose() * U * psi.psi(1);
        u1 b = psi_conj(1).transpose() * U * psi.psi(0);
        u1 c = psi_conj(2).transpose() * U * psi.psi(3);
        u1 d = psi_conj(3).transpose() * U * psi.psi(2);
        u = (-a - b + c + d) * I;
    } else if (gamma == 2) {
        u1 a = psi_conj(0).transpose() * U * psi.psi(1);
        u1 b = psi_conj(1).transpose() * U * psi.psi(0);
        u1 c = psi_conj(2).transpose() * U * psi.psi(3);
        u1 d = psi_conj(3).transpose() * U * psi.psi(2);
        u = (a - b - c + d);
    } else if (gamma == 3) {
        u1 a = psi_conj(0).transpose() * U * psi.psi(0);
        u1 b = psi_conj(1).transpose() * U * psi.psi(1);
        u1 c = psi_conj(2).transpose() * U * psi.psi(2);
        u1 d = psi_conj(3).transpose() * U * psi.psi(3);
        u = (-a + b + c - d) * I;
    } else if (gamma == 5) {
        u1 a = psi_conj(0).transpose() * U * psi.psi(2);
        u1 b = psi_conj(1).transpose() * U * psi.psi(3);
        u1 c = psi_conj(2).transpose() * U * psi.psi(0);
        u1 d = psi_conj(3).transpose() * U * psi.psi(1);
        return (a + b - c - d);
    }
    return u;
}

bool qcd_site::lock() {
//    return false;
    if (is_locked || forward[0]->is_locked || backward[0]->is_locked) {
        this_thread::sleep_for(chrono::milliseconds(10));
        return true;
    }
    is_locked = true;
    forward[0]->is_locked = true;
    backward[0]->is_locked = true;
    return false;
}

void qcd_site::unlock() {
    backward[0]->is_locked = false;
    forward[0]->is_locked = false;
    is_locked = false;
}

void qcd_site::reset_site(bool cold) {
    for (int d = 0; d < lattice->D; d++) set_link(d, cold ? su3_identity : create_link(true));
    for (int f = 0; f < lattice->F; f++) psi[f] = create_spinor(!cold);
}

void qcd_site::copy_site(qcd_site* site) {
    for (int d = 0; d < lattice->D; d++) set_link(d, site->link[d]);
    for (int f = 0; f < lattice->F; f++) psi[f] = qcd_spinor(site->psi[f]);
}

su2 qcd_site::create_su2(bool random, double z) {
    
    double a, g;
    
    if (random) {
        a = lattice->rand(-1.0, 1.0);
        g = sqrt(1 - a * a);
    } else {
        double x;
        
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
    
    su2 m;
    m(0,0) = u1(a0, a3);
    m(0,1) = u1(a2, a1);
    m(1,0) = u1(-a2, a1);
    m(1,1) = u1(a0, -a3);
    return m;
}

su3 qcd_site::create_link(bool random, double z) {
    
    su3 R = su3_identity;
    R.block<2,2>(0,0) = create_su2(random, z);
    
    su3 S = su3_identity;
    su2 s = create_su2(random, z);
    S(0,0) = s(0,0);
    S(0,2) = s(0,1);
    S(2,0) = s(1,0);
    S(2,2) = s(1,1);
    
    su3 T = su3_identity;
    T.block<2,2>(1,1) = create_su2(random, z);
    
    su3 RST = R * S * T;
    if (lattice->rand() > 0.5) return RST.transpose().conjugate();
    return RST;
}

u1 qcd_site::create_u1(bool random) {
    if (!random) return u1(1, 0);
    double phi = lattice->rand(-M_PI, M_PI);
    return u1(cos(phi), sin(phi));
}

su3_vector qcd_site::random_vector() {
    u1 a1 = create_u1(true);
    u1 a2 = create_u1(true);
    u1 a3 = create_u1(true);
    double a = sqrt(abs(a1) + abs(a2) + abs(a3));
    a1 /= a; a2 /= a; a3 /= a;
    return su3_vector(a1, a2, a3);
}

qcd_spinor qcd_site::create_spinor(bool random, double z) {
    qcd_spinor psi;
//    if (random) {
//        for (int c = 0; c < 24; c++) psi.component[c] = lattice->rand(-1.0, 1.0);
//    }
    return psi;
}

void qcd_site::set_link(int d, su3 value) {
    u1 det = value.determinant();
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

double qcd_site::action() {
    
    double S = gauge_action();
    for (int f = 0; f < lattice->F; f++) S += real(fermion_action(f));
    return S;
}

double qcd_site::gauge_action() {
    // gauge action
    double S_G = 0.0;
    su3 s;
    for (int d1 = 0; d1 < lattice->D; d1++) {
        for (int d2 = d1 + 1; d2 < lattice->D; d2++) {
            s = link[d1];
            s *= forward[d1]->link[d2];
            s *= forward[d2]->link_inverse[d1];
            s *= link_inverse[d2];
            S_G += s.trace().real() / 3.0;
        }
    }
//    S_G = S_G / 6.0;
    return lattice->beta * (1.0 - S_G);
}

u1 qcd_site::psi_bar_psi(int f) {
    return psi[f].contract(psi[f]);
}

u1 qcd_site::fermion_action(int f) {
    u1 S_F = 0;
    for (int d = 0; d < lattice->D; d++) {
        S_F += psi[f].contract(link[d], forward[d]->psi[f]) * forward_edge[d];
        S_F -= psi[f].contract(d, link[d], forward[d]->psi[f]) * forward_edge[d];
        S_F += forward[d]->psi[f].contract(link_inverse[d], psi[f]) * forward_edge[d];
        S_F += forward[d]->psi[f].contract(d, link_inverse[d], psi[f]) * forward_edge[d];
    }
    
    return psi[f].contract(psi[f]) - lattice->kappa[f] * S_F;
}

double qcd_site::site_fermion_action(int f, qcd_spinor t_psi) {
    
//    qcd_spinor s_psi = psi[f];
//    psi[f] = t_psi;
//    u1 S_F = fermion_action(f); for (int d = 0; d < lattice->D; d++) S_F += backward[d]->fermion_action(f);
//    psi[f] = s_psi;
//
//    return real(S_F);

    u1 S_F = 0;
    for (int d = 0; d < lattice->D; d++) {
        S_F += t_psi.contract(link[d], forward[d]->psi[f]) * forward_edge[d];
        S_F -= t_psi.contract(d, link[d], forward[d]->psi[f]) * forward_edge[d];
        S_F += t_psi.contract(backward[d]->link_inverse[d], backward[d]->psi[f]) * backward_edge[d];
        S_F += t_psi.contract(d, backward[d]->link_inverse[d], backward[d]->psi[f]) * backward_edge[d];

        S_F += backward[d]->psi[f].contract(backward[d]->link[d], t_psi) * backward_edge[d];
        S_F -= backward[d]->psi[f].contract(d, backward[d]->link[d], t_psi) * backward_edge[d];
        S_F += forward[d]->psi[f].contract(link_inverse[d], t_psi) * forward_edge[d];
        S_F += forward[d]->psi[f].contract(d, link_inverse[d], t_psi) * forward_edge[d];
    }

    return real(t_psi.contract(t_psi) - lattice->kappa[f] * S_F);
}

double qcd_site::sweep_gauge() {
    double accept = 0;
    while (lock());
    for (int d = 0; d < lattice->D; d++) accept += sweep_link(d);
//    for (int f = 0; f < lattice->F; f++) sweep_fermion(f);
    unlock();
    return accept / lattice->D;
}

double qcd_site::sweep_link(int d1) {
    su3 old_link = link[d1];
    su3 old_link_inverse = link_inverse[d1];
    su3 new_link = old_link * create_link(false, lattice->z_G);
    su3 new_link_inverse = new_link.transpose().conjugate();

    su3 u;
    su3 A = su3_zero;
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
    
    double dS = -lattice->beta * ((new_link - old_link) * A).trace().real() / 3.0;
    
    for (int f = 0; f < lattice->F; f++) {
        u1 dS_F = 0;
        dS_F += psi[f].contract(new_link, forward[d1]->psi[f]) * forward_edge[d1];
        dS_F -= psi[f].contract(d1, new_link, forward[d1]->psi[f]) * forward_edge[d1];
        dS_F += forward[d1]->psi[f].contract(new_link_inverse, psi[f]) * forward_edge[d1];
        dS_F += forward[d1]->psi[f].contract(d1, new_link_inverse, psi[f]) * forward_edge[d1];

        dS_F -= psi[f].contract(old_link, forward[d1]->psi[f]) * forward_edge[d1];
        dS_F += psi[f].contract(d1, old_link, forward[d1]->psi[f]) * forward_edge[d1];
        dS_F -= forward[d1]->psi[f].contract(old_link_inverse, psi[f]) * forward_edge[d1];
        dS_F -= forward[d1]->psi[f].contract(d1, old_link_inverse, psi[f]) * forward_edge[d1];

        dS -= lattice->kappa[f] * real(dS_F);
    }

    // metropolis algorithm
    if (dS <= 0 || lattice->rand() < exp(-dS)) {
        set_link(d1, new_link);
        return 1.0;
    }
    return 0.0;
}

double qcd_site::sweep_fermion(int f) {
    while (lock());
    bool verbose = (lattice->rand() * lattice->n_sites < 5.0);

    // copy the spinor value
    qcd_spinor t_psi(psi[f]);

    // step distance
//    double h = 1.0;
    double h = exp(-lattice->z_F[f]);

    double steps = 0.0;
    double old_S = site_fermion_action(f, t_psi);

    double S = old_S;
    Eigen::Matrix<double, 24, 1> d0;
    Eigen::Matrix<double, 24, 1> gradient;
    Eigen::Matrix<double, 24, 24> inverse_hessian;
    for (int c = 0; c < 12; c++) inverse_hessian(c,c+12) = inverse_hessian(c+12,c) = 0.5;

//    while (true) {
        double S_b[24];
        double S_f[24];
        
        // set the initial spinor values
        for (int c = 0; c < 24; c++) d0(c) = t_psi.component[c];

        // compute the gradient
        for (int c = 0; c < 24; c++) {
            t_psi.component[c] = d0(c) - h;
            S_b[c] = site_fermion_action(f, t_psi);
            t_psi.component[c] = d0(c) + h;
            S_f[c] = site_fermion_action(f, t_psi);
            t_psi.component[c] = d0(c);
            
            gradient(c) = (S_f[c] - S_b[c]) * 0.5;
        }

//        if (verbose) cout << "step = " << int(steps) << ", gradient_norm = " << gradient.norm() / h << "\n";
//        if (verbose) cout << "gradient = \n" << gradient.transpose() / h << "\n";

        // exit the loop when the gradient gets small enough
//        if (gradient.norm() / h < 1.0 || steps > 5.0) break;

        // set the optimized spinor values
        d0 -= (inverse_hessian * gradient) / h * 0.1;
//        h *= 0.5;
        for (int c = 0; c < 24; c++) t_psi.component[c] = d0(c);

        // calculate the new action
        S = site_fermion_action(f, t_psi);
        steps += 1.0;
//    }

    // metropolis algorithm
    double dS = real(S - old_S);
    if (dS <= 0 || lattice->rand() < exp(-dS)) {
        psi[f] = t_psi;
//        unlock();
//        return 1.0;
    }
    
    unlock();
    return gradient.norm();
}

qcd_lattice::qcd_lattice(int N, int T, int D, int F, double kappa, double beta, double* m, bool cold) {
    this->N = N;
    this->T = T;
    this->D = D;
    this->F = F;
    this->beta = beta;
    for (int f = 0; f < this->F; f++) this->m[f] = m[f];
    for (int f = 0; f < this->F; f++) this->kappa[f] = kappa / (m[f] - 8.0 * kappa * (m[f] - 1.0));
    this->z_G = 50.0;
    for (int f = 0; f < this->F; f++) this->z_F[f] = 2.3;

    random_device rd;  // seed for the random number engine
    gen = mt19937(rd()); // mersenne_twister_engine seeded with rd()
    
    n_sites = pow(N, D - 1) * T;
    n_slice = n_sites / T; // number of sites in each time slice
    site = new qcd_site[n_sites];
    
    for (int s = 0; s < n_sites; s++) {
        site[s].lattice = this;
        site[s].is_locked = NULL;
        site[s].reset_site(cold);
        
        for (int d = 0; d < D; d++) {
            update_i(s);
            int edge = (d == 0 ? T : N) - 1;
            site[s].forward_edge[d] = (i[d] == edge ? -1.0 : 1.0);
            site[s].backward_edge[d] = (i[d] == 0 ? -1.0 : 1.0);
            move(d, +1);
            site[s].forward[d] = &site[update_s()];
            move(d, -2);
            site[s].backward[d] = &site[update_s()];
        }
    }
}

qcd_lattice::qcd_lattice(qcd_lattice* lattice) {
    this->N = lattice->N;
    this->T = lattice->T;
    this->D = lattice->D;
    this->F = lattice->F;
    this->beta = lattice->beta;
    for (int f = 0; f < this->F; f++) this->m[f] = lattice->m[f];
    for (int f = 0; f < this->F; f++) this->kappa[f] = lattice->kappa[f];
    this->parallel = lattice->parallel;
    this->z_G = lattice->z_G;
    for (int f = 0; f < this->F; f++) this->z_F[f] = lattice->z_F[f];

    random_device rd;  // seed for the random number engine
    gen = mt19937(rd()); // mersenne_twister_engine seeded with rd()
    
    n_sites = pow(N, D - 1) * T;
    n_slice = n_sites / T; // number of sites in each time slice
    site = new qcd_site[n_sites];
    
    for (int s = 0; s < n_sites; s++) {
        site[s].lattice = this;
        site[s].is_locked = NULL;
        site[s].copy_site(&lattice->site[s]);
        
        for (int d = 0; d < D; d++) {
            update_i(s);
            int edge = (d == 0 ? T : N) - 1;
            site[s].forward_edge[d] = (i[d] == edge ? -1.0 : 1.0);
            site[s].backward_edge[d] = (i[d] == 0 ? -1.0 : 1.0);
            move(d, +1);
            site[s].forward[d] = &site[update_s()];
            move(d, -2);
            site[s].backward[d] = &site[update_s()];
        }
    }
}

qcd_lattice::~qcd_lattice() {
    delete[] site;
}

double qcd_lattice::rand(double min, double max) {
    return ((double)gen()) / 4294967295.0 * (max - min) + min;
}

void qcd_lattice::update_i(int s) {
    i[0] = s / (int)pow(N, D - 1);
    for (int d = 1; d < D; d++) i[d] = s / (int)pow(N, D - d - 1) % N;
}

int qcd_lattice::update_s() {
    int s = 0;
    for (int d = 0; d < D; d++) s += i[d] * (int)pow(N, D - d - 1);
    return s;
}

void qcd_lattice::move(int d, int n) {
    
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

double qcd_lattice::plaquette() {
    return (1 - gauge_action() / beta) / 2;
}

double qcd_lattice::action() {
    double S = 0;
    for (int s = 0; s < n_sites; s++) S += site[s].action();
    return S / double(n_sites);
}

double qcd_lattice::gauge_action() {
    double S = 0;
    for (int s = 0; s < n_sites; s++) S += site[s].gauge_action();
    return S / double(n_sites);
}

double qcd_lattice::psi_bar_psi(int f) {
    u1 P = 0;
    for (int s = 0; s < n_sites; s++) P += site[s].psi_bar_psi(f);
//    cout << "psi_bar_psi = " << abs(imag(F) / real(F)) << "\n";
    return real(P) / double(n_sites);
}

double qcd_lattice::fermion_action() {
    u1 S = 0;
    for (int s = 0; s < n_sites; s++) {
        for (int f = 0; f < F; f++) S += site[s].fermion_action(f);
    }
//    cout << "F = " << abs(imag(S) / real(S)) << "\n";
    return real(S) / double(n_sites * F);
}

double qcd_lattice::xi_avg() {
    double xi = 0.0;
    for (int s = 0; s < n_sites; s++) {
        for (int f = 0; f < F; f++) {
            double psi = 0.0;
            for (int c = 0; c < 24; c++) {
                psi += site[s].psi[f].component[c] * site[s].psi[f].component[c];
            }
            xi += sqrt(psi);
        }
    }
    return xi / n_sites / F;
}

double async_sweep_G(qcd_site* site, int n) {
    double accept = 0.0;

    // start at a random spot in the time slice to mitigate parallel effects
    int start = rand() % n;
    int end = start + n;
    
    for (int s = start; s < end; s++) accept += site[s % n].sweep_gauge();

    return accept;
}

double async_sweep_F(qcd_site* site, int n, int f) {
    double accept = 0.0;

    // start at a random spot in the time slice to mitigate parallel effects
    int start = rand() % n;
    int end = start + n;
    
    for (int s = start; s < end; s++) accept += site[s % n].sweep_fermion(f);

    return accept;
}

void qcd_lattice::sweep(int n_sweeps, bool gauge_only) {
    sweep_gauge();
    sweep_f();
}
void qcd_lattice::sweep_gauge(int n_sweeps) {
    
    if (n_sweeps != 1) {
        for (int i = 0; i < n_sweeps; i++) sweep_gauge(1);
        return;
    }
    
    double accept_G = 0;

    if (parallel) {
        future<double> sweep_slice_G[T];
        for (int t = 0; t < T; t++) sweep_slice_G[t] = async(async_sweep_G, &site[t * n_slice], n_slice);
        for (int t = 0; t < T; t++) accept_G += sweep_slice_G[t].get();
        
    } else {
        for (int s = 0; s < n_sites; s++) accept_G += site[s].sweep_gauge();
    }
    accept_G /= n_sites;

//    cout << "accept_G = " << accept_G << ", z_G = " << z_G << "\n";
    if (accept_G < 0.78) {
        z_G *= 1.1;
    } else if (accept_G > 0.82) {
        z_G *= 0.9;
    }
}

double qcd_lattice::sweep_f(int n_sweeps) {
    
    if (n_sweeps != 1) {
        double a = 0.0;
        for (int i = 0; i < n_sweeps; i++) a += sweep_f(1);
        return a;
    }
    
    double accept[F];
    for (int f = 0; f < F; f++) {
        accept[f] = 0.0;
    }

    if (parallel) {
        future<double> sweep_slice[T][F];
        for (int t = 0; t < T; t++) for (int f = 0; f < F; f++) sweep_slice[t][f] = async(async_sweep_F, &site[t * n_slice], n_slice, f);
        for (int t = 0; t < T; t++) for (int f = 0; f < F; f++) accept[f] += sweep_slice[t][f].get();
    } else {
        for (int s = 0; s < n_sites; s++) {
            for (int f = 0; f < F; f++) accept[f] += site[s].sweep_fermion(f);
        }
    }
    
    for (int f = 0; f < F; f++) accept[f] /= n_sites;
//    cout << "accept_F = " << accept[0] << ", z_F = " << z_F[0] << "\n";
//    for (int f = 0; f < F; f++) {
//
//        if (accept[f] < 0.78) {
//            z_F[f] *= 0.9;
//        } else if (accept[f] > 0.82) {
//            z_F[f] *= 1.1;
//        }
//    }
    return accept[0];
}

int qcd_lattice::thermalize(int n_min, int n_max) {
    
    // create a cold start model
    qcd_lattice coldStart = qcd_lattice(N, T, D, F, kappa[0], beta, m, true);
    double hot_plaquette = plaquette();
    double cold_plaquette = coldStart.plaquette();
    cout << "Average Hot Plaquette: " << hot_plaquette << "\n";
    cout << "Average Cold Plaquette: " << cold_plaquette << "\n";
    double hot_fermion = psi_bar_psi(0);
    double cold_fermion = coldStart.psi_bar_psi(0);
    cout << "Average Hot Fermion: " << hot_fermion << "\n";
    cout << "Average Cold Fermion: " << cold_fermion << "\n";
    
    // sweep the gauge fields on both models until the average plaquette reaches the same value
    int m;
    for (m = 0; hot_plaquette < cold_plaquette; m++) {
        sweep_gauge();
        coldStart.sweep_gauge();
        hot_plaquette = plaquette();
        cold_plaquette = coldStart.plaquette();
        hot_fermion = fermion_action();
        cold_fermion = coldStart.fermion_action();
        cout << "dP = " << abs(hot_plaquette - cold_plaquette);
        cout << ", P_hot = " << hot_plaquette;
        cout << ", P_cold = " << cold_plaquette;
        cout << ", F_hot = " << hot_fermion;
        cout << ", F_cold = " << cold_fermion << "\n";
        if (n_max && (m > n_max)) break;
    }
    
    // sweep the fermion field
    int n_therm = max(n_min, m);
    for (int i = 0; i < n_therm; i++) {
        double poo = sweep_f();
        if (poo < 0.1) sweep_gauge();
//        sweep_f();
//        sweep_f();
//        sweep_f();
//        sweep_f();
//        sweep_f();
//        sweep_gauge();
        cout << "P = " << plaquette();
        cout << "grad = " << poo;
        cout << ", xi = " << xi_avg();
        cout << ", psi_bar_psi = " << psi_bar_psi(0);
        cout << ", F = " << fermion_action() << "\n";
    }

    // sweep the gauge field and the fermion field
    for (int i = 0; i < n_therm; i++) {
        sweep_gauge();
        sweep_f();
        cout << "P = " << plaquette();
        cout << ", xi = " << xi_avg();
        cout << ", psi_bar_psi = " << psi_bar_psi(0);
        cout << ", F = " << fermion_action() << "\n";
    }

//    cout << "aveS = " << action() << "\n";
    
    // return value is total number of sweeps
    return n_therm;
}
