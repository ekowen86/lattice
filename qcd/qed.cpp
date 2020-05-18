//
//  qed.cpp
//  Lattice
//
//  Created by Evan Owen on 2/16/20.
//  Copyright © 2020 Evan Owen. All rights reserved.
//

#include <iostream>
#include <random>
#include <cmath>
#include <future>
#include "Eigen/Dense"
#include "qed.hpp"

using namespace std;

qed_spinor::qed_spinor() {
}

qed_spinor::qed_spinor(qed_psi* psi) {
    set_psi(psi);
}

void qed_spinor::set_psi(qed_psi* psi) {
    this->psi[0] = psi[0];
    this->psi[1] = psi[1];
    this->psi[2] = psi[2];
    this->psi[3] = psi[3];
    this->psi_bar[0] = psi[2].conj();
    this->psi_bar[1] = psi[3].conj();
    this->psi_bar[2] = psi[0].conj();
    this->psi_bar[3] = psi[1].conj();
}

//qed_spinor qed_spinor::operator+ (qed_spinor spinor) {
//    qed_psi psi[4];
//    for (int i = 0; i < 4; i++) psi[i] = this->psi[i] + spinor.psi[i];
//    return qed_spinor(psi);
//}

qed_spinor qed_spinor::operator* (qed_link u) {
    qed_psi psi[4];
    for (int i = 0; i < 4; i++) psi[i] = this->psi[i] * u;
    return qed_spinor(psi);
}

qed_spinor qed_spinor::operator* (int gamma) {
    qed_psi psi[4];
    switch (gamma) {
        case 0:
            psi[0] = this->psi[2];
            psi[1] = this->psi[3];
            psi[2] = this->psi[0];
            psi[3] = this->psi[1];
            break;
        case 1:
            psi[0] = this->psi[3];
            psi[1] = this->psi[2];
            psi[2] = -this->psi[1];
            psi[3] = -this->psi[0];
            break;
        case 2:
            psi[0] = this->psi[3] * -I;
            psi[1] = this->psi[2] * I;
            psi[2] = this->psi[1] * I;
            psi[3] = this->psi[0] * -I;
            break;
        case 3:
            psi[0] = this->psi[2];
            psi[1] = -this->psi[3];
            psi[2] = -this->psi[0];
            psi[3] = this->psi[1];
            break;
        case 5:
            psi[0] = -this->psi[0];
            psi[1] = -this->psi[1];
            psi[2] = this->psi[2];
            psi[3] = this->psi[3];
            break;
    }
    return qed_spinor(psi);
}

qed_link qed_spinor::operator* (qed_spinor spinor) {
    qed_link u = 0.0;
    for (int i = 0; i < 4; i++) u += psi_bar[i] * spinor.psi[i];
    return u;
}

bool qed_site::lock() {
    if (is_locked || forward[0]->is_locked || backward[0]->is_locked) {
        this_thread::sleep_for(chrono::milliseconds(10));
        return false;
    }
    is_locked = true;
    forward[0]->is_locked = true;
    backward[0]->is_locked = true;
    return true;
}

void qed_site::unlock() {
    is_locked = false;
    forward[0]->is_locked = false;
    backward[0]->is_locked = false;
}

void qed_site::reset_site(bool cold) {
    for (int d = 0; d < lattice->D; d++) set_link(d, cold ? qed_identity : create_link(true));
    for (int f = 0; f < lattice->F; f++) spinor[f] = create_spinor(!cold);
}

void qed_site::copy_site(qed_site* site) {
    for (int d = 0; d < lattice->D; d++) set_link(d, site->link[d]);
    for (int f = 0; f < lattice->F; f++) spinor[f].set_psi(site->spinor[f].psi);
}

qed_link qed_site::create_link(bool random, double z) {
    if (!random) return qed_link(1, 0);
    double phi = lattice->rand(-M_PI * z, M_PI * z);
    return qed_link(cos(phi), sin(phi));
}

//qed_vector qed_site::create_vector(bool random, double z) {
//    if (random) {
//        return create_link(true, z);
//    } else {
//        return qed_zero;
//    }
//}

qed_spinor qed_site::create_spinor(bool random, double z) {
    qed_psi psi[4];
    if (random) {
        psi[0] = qed_psi(create_link(true, z), create_link(true, z));
        psi[1] = qed_psi(create_link(true, z), create_link(true, z));
        psi[2] = qed_psi(create_link(true, z), create_link(true, z));
        psi[3] = qed_psi(create_link(true, z), create_link(true, z));
    } else {
        psi[0] = qed_psi(1, 1);
        psi[1] = qed_psi(1, 1);
        psi[2] = qed_psi(1, 1);
        psi[3] = qed_psi(1, 1);
    }
    return qed_spinor(psi);
}

void qed_site::set_link(int d, qed_link value) {
    double v = norm(value);
    if (abs(v - 1.0) > 1e-10) {
        // re-unitarize the value
        value /= v;
    }
    link[d] = value;
    link_inverse[d] = conj(value);
}

qed_link qed_site::action() {
    
    qed_link S = gauge_action();
    for (int f = 0; f < lattice->F; f++) S += fermion_action(f);
    return S;
}

qed_link qed_site::gauge_action() {
    // gauge action
    double S_G = 0;
    qed_link s;
    for (int d1 = 0; d1 < lattice->D; d1++) {
        for (int d2 = d1 + 1; d2 < lattice->D; d2++) {
            s = link[d1];
            s *= forward[d1]->link[d2];
            s *= forward[d2]->link_inverse[d1];
            s *= link_inverse[d2];
            S_G += s.real();
        }
    }
    S_G = S_G / lattice->D / (lattice->D - 1) * 2;
    return lattice->beta * (1 - S_G);
}

qed_link qed_site::fermion_action(int f) {
    qed_link S_F = 0;
    for (int d = 0; d < lattice->D; d++) {
        S_F += spinor[f] * d * link[d] * forward[d]->spinor[f] * forward_edge[d];
        S_F -= spinor[f] * link[d] * forward[d]->spinor[f] * forward_edge[d];
        S_F -= spinor[f] * d * backward[d]->link_inverse[d] * backward[d]->spinor[f] * backward_edge[d];
        S_F -= spinor[f] * backward[d]->link_inverse[d] * backward[d]->spinor[f] * backward_edge[d];
    }
    double a = lattice->a;
    S_F /= (2 * a);
    S_F += (lattice->m[f] + 4 / a) * (spinor[f] * spinor[f]);
    return S_F * a * a * a * a;
}

double qed_site::sweep() {
    double accept = 0;
    lock();
    for (int d = 0; d < lattice->D; d++) accept += sweep_link(d);
    // sweep fermion
//    accept += sweep_fermion();
    unlock();
    return accept / lattice->D;
}

double qed_site::sweep_link(int d1) {
    qed_link old_link = link[d1];
    qed_link new_link = old_link * create_link(true, lattice->z_G);
    
    qed_link u;
    qed_link A = qed_zero;
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
    
    double dE = real((new_link - old_link) * A);
    
    // metropolis algorithm
    if (dE >= 0 || lattice->rand() < exp(lattice->beta * dE)) {
        set_link(d1, new_link);
        return 1;
    }
    return 0;
}

double qed_site::sweep_fermion(int f) {
    lock();
    qed_spinor old_spinor = spinor[f];
    qed_spinor new_spinor(spinor[f].psi);
    qed_psi* new_psi = new_spinor.psi;
    new_psi[0].v1 *= create_link(true, lattice->z_F[f]);
    new_psi[0].v2 *= create_link(true, lattice->z_F[f]);
    new_psi[1].v1 *= create_link(true, lattice->z_F[f]);
    new_psi[1].v2 *= create_link(true, lattice->z_F[f]);
    new_psi[2].v1 *= create_link(true, lattice->z_F[f]);
    new_psi[2].v2 *= create_link(true, lattice->z_F[f]);
    new_psi[3].v1 *= create_link(true, lattice->z_F[f]);
    new_psi[3].v2 *= create_link(true, lattice->z_F[f]);
    new_spinor.set_psi(new_psi);
    
    qed_link old_S = fermion_action(f);
    for (int d = 0; d < lattice->D; d++) {
        old_S += forward[d]->fermion_action(f);
        old_S += backward[d]->fermion_action(f);
    }
    spinor[f] = new_spinor;
    qed_link new_S = fermion_action(f);
    for (int d = 0; d < lattice->D; d++) {
        new_S += forward[d]->fermion_action(f);
        new_S += backward[d]->fermion_action(f);
    }

    double dS = real(new_S - old_S);
    
    // metropolis algorithm
    if (dS >= 0 || lattice->rand() < exp(dS)) {
        unlock();
        return 1;
    } else {
        spinor[f] = old_spinor;
        unlock();
        return 0;
    }
}

qed_lattice::qed_lattice(int N, int T, int D, int F, double a, double beta, double* m, bool cold) {
    this->N = N;
    this->T = T;
    this->D = D;
    this->F = F;
    this->a = a;
    this->beta = beta;
    for (int f = 0; f < this->F; f++) this->m[f] = m[f];
    this->z_G = 0.5;
    for (int f = 0; f < this->F; f++) this->z_F[f] = 0.5;

    random_device rd;  // seed for the random number engine
    gen = mt19937(rd()); // mersenne_twister_engine seeded with rd()
    
    n_sites = pow(N, D - 1) * T;
    n_slice = n_sites / T; // number of sites in each time slice
    site = new qed_site[n_sites];
    
    for (int s = 0; s < n_sites; s++) {
        site[s].lattice = this;
        site[s].is_locked = NULL;
        site[s].reset_site(cold);
        
        for (int d = 0; d < D; d++) {
            update_i(s);
            int edge = (d == 0 ? T : N) - 1;
            site[s].forward_edge[d] = (i[d] == edge ? -1 : 1);
            site[s].backward_edge[d] = (i[d] == 0 ? -1 : 1);
            move(d, +1);
            site[s].forward[d] = &site[update_s()];
            move(d, -2);
            site[s].backward[d] = &site[update_s()];
        }
    }
}

qed_lattice::qed_lattice(qed_lattice* lattice) {
    this->N = lattice->N;
    this->T = lattice->T;
    this->D = lattice->D;
    this->F = lattice->F;
    this->a = lattice->a;
    this->beta = lattice->beta;
    for (int f = 0; f < this->F; f++) this->m[f] = lattice->m[f];
    this->parallel = lattice->parallel;
    this->z_G = 0.5;
    for (int f = 0; f < this->F; f++) this->z_F[f] = 0.5;

    random_device rd;  // seed for the random number engine
    gen = mt19937(rd()); // mersenne_twister_engine seeded with rd()
    
    n_sites = pow(N, D - 1) * T;
    n_slice = n_sites / T; // number of sites in each time slice
    site = new qed_site[n_sites];
    
    for (int s = 0; s < n_sites; s++) {
        site[s].lattice = this;
        site[s].is_locked = NULL;
        site[s].copy_site(&lattice->site[s]);
        
        for (int d = 0; d < D; d++) {
            update_i(s);
            int edge = (d == 0 ? T : N) - 1;
            site[s].forward_edge[d] = (i[d] == edge ? -1 : 1);
            site[s].backward_edge[d] = (i[d] == 0 ? -1 : 1);
            move(d, +1);
            site[s].forward[d] = &site[update_s()];
            move(d, -2);
            site[s].backward[d] = &site[update_s()];
        }
    }
}

qed_lattice::~qed_lattice() {
    delete[] site;
}

double qed_lattice::rand(double min, double max) {
    return ((double)gen()) / 4294967295.0 * (max - min) + min;
}

void qed_lattice::update_i(int s) {
    i[0] = s / (int)pow(N, D - 1);
    for (int d = 1; d < D; d++) i[d] = s / (int)pow(N, D - d - 1) % N;
}

int qed_lattice::update_s() {
    int s = 0;
    for (int d = 0; d < D; d++) s += i[d] * (int)pow(N, D - d - 1);
    return s;
}

void qed_lattice::move(int d, int n) {
    
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

double qed_lattice::action() {
    qed_link S = 0;
    for (int s = 0; s < n_sites; s++) S += site[s].action();
//    cout << abs(S.imag() / S.real()) << "\n";
    return real(S) / n_sites;
}

double async_sweep_G(qed_site* site, int n) {
    double accept = 0.0;

    // start at a random spot in the time slice to mitigate parallel effects
    int start = rand() % n;
    int end = start + n;
    
    for (int s = start; s < end; s++) accept += site[s % n].sweep();

    return accept;
}

double async_sweep_F(qed_site* site, int n, int f) {
    double accept = 0.0;

    // start at a random spot in the time slice to mitigate parallel effects
    int start = rand() % n;
    int end = start + n;
    
    for (int s = start; s < end; s++) accept += site[s % n].sweep_fermion(f);

    return accept;
}

void qed_lattice::sweep(int n_sweeps) {
    
    if (n_sweeps) {
        for (int i = 0; i < n_sweeps; i++) sweep(0);
        return;
    }
    
    double accept_G = 0;
    double accept_F[F];
    for (int f = 0; f < F; f++) accept_F[f] = 0.0;

    if (parallel) {
        future<double> sweep_slice_G[T];
        future<double> sweep_slice_F[T][F];
        for (int t = 0; t < T; t++) {
            sweep_slice_G[t] = async(async_sweep_G, &site[t * n_slice], n_slice);
            for (int f = 0; f < F; f++) sweep_slice_F[t][f] = async(async_sweep_F, &site[t * n_slice], n_slice, f);
        }
        for (int t = 0; t < T; t++) {
            accept_G += sweep_slice_G[t].get();
            for (int f = 0; f < F; f++) accept_F[f] += sweep_slice_F[t][f].get();
        }
    } else {
        for (int s = 0; s < n_sites; s++) {
            accept_G += site[s].sweep();
            for (int f = 0; f < F; f++) accept_F[f] += site[s].sweep_fermion(f);
        }
    }
    accept_G /= n_sites;
    for (int f = 0; f < F; f++) accept_F[f] /= n_sites;

//    cout << "accept_G = " << accept_G << ", z_G = " << z_G << "\n";
    if (accept_G < 0.78) {
        z_G *= 0.95;
    } else if (accept_G > 0.82) {
        z_G *= 1.05;
        if (z_G > 1.0) z_G = 1.0;
    }
    
//    cout << "accept_F = " << accept_F << ", z_F = " << z_F << "\n";
    for (int f = 0; f < F; f++) {
        
        if (accept_F[f] < 0.78) {
            z_F[f] *= 0.95;
        } else if (accept_F[f] > 0.82) {
            z_F[f] *= 1.05;
            if (z_F[f] > 1.0) z_F[f] = 1.0;
        }
    }
}

int qed_lattice::thermalize(int n_min, int n_max) {
    
    // create a cold start model
    qed_lattice coldStart = qed_lattice(N, T, D, F, a, beta, m, true);
    double hot_action = action();
    double cold_action = coldStart.action();
    int m;
    for (m = 0; hot_action > cold_action; m++) {
        // sweep both models until the average plaquette reaches the same value
        sweep();
        coldStart.sweep();
        hot_action = action();
        cold_action = coldStart.action();
        cout << "dE = " << abs(cold_action - hot_action) << "\n";
        if (n_max && (m > n_max)) break;
    }
    
    // sweep a few more times
    int n_therm = max(n_min, m * 2);
    for (int i = m; i < n_therm; i++) sweep();
    
    cout << "Eave = " << action() << "\n";
    
    // return value is total number of sweeps
    return n_therm;
}
