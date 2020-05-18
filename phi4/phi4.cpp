//
//  phi4.cpp
//
//  Created by Evan Owen on 11/7/19.
//  Copyright © 2019 Evan Owen. All rights reserved.
//

#include <iostream>
#include <iomanip>
#include <fstream>
#include <random>
#include <cmath>
#include <future>
#include "phi4.hpp"

using namespace std;

bool phi4_site::lock() {
    if (is_locked || forward[0]->is_locked || backward[0]->is_locked) {
        this_thread::sleep_for(chrono::milliseconds(10));
        return false;
    }
    is_locked = true;
    forward[0]->is_locked = true;
    backward[0]->is_locked = true;
    return true;
}

void phi4_site::unlock() {
    is_locked = false;
    forward[0]->is_locked = false;
    backward[0]->is_locked = false;
}

void phi4_site::reset_value(bool cold) {
    value = cold ? 0.0 : lattice->rand(0, 1);
}

//void phi4_site::copy_links(phi4_site* site) {
//    for (int d = 0; d < lattice->D; d++) set_link(d, phi4_link(site->link[d]));
//}

void phi4_site::set_value(double value) {
    this->value = value;
}

//phi4_link phi4_site::create_link(bool random) {
//
//    su3_link R = su3_identity;
//    R.block<2,2>(0,0) = create_su2(random);
//
//    su3_link S = su3_identity;
//    su2_link s = create_su2(random);
//    S(0,0) = s(0,0);
//    S(0,2) = s(0,1);
//    S(2,0) = s(1,0);
//    S(2,2) = s(1,1);
//
//    su3_link T = su3_identity;
//    T.block<2,2>(1,1) = create_su2(random);
//
//    su3_link RST = R * S * T;
//    if (rand() > 0.5) return RST.transpose().conjugate();
//    return RST;
//}

double phi4_site::action() {
    // compute the contribution to the lattice action at this site
    double S = 0;
    double s = 0;
    for (int d = 0; d < lattice->D; d++) {
        s = this->forward[d]->value - value;
        S += s * s / 2;
    }
//    S *= lattice->beta * S / 2;
    S += lattice->mu * value * value / 2;
    S += lattice->lambda * value * value * value * value / 2;
    return S;
}

double phi4_site::sweep() {
    double accept = 0;
    lock();
    for (int d = 0; d < lattice->D; d++) accept += sweep_link(d);
    unlock();
    return accept / lattice->D;
}

double phi4_site::sweep_link(int d1) {
    double old_val = value;
    double old_val2 = old_val * old_val;
    double old_val4 = old_val2 * old_val2;

    double new_val = lattice->rand(0, 1);
//    double new_val = lattice->rand(-1, 1) * lattice->z + old_val;
    if (new_val < 0) new_val = 0.0;
    if (new_val > 1) new_val = 1.0;
    double new_val2 = new_val * new_val;
    double new_val4 = new_val2 * new_val2;

    double delta = new_val - old_val;
    double delta2 = new_val2 - old_val2;
    double delta4 = new_val4 - old_val4;

    double dS = 0;
    for (int d = 0; d < lattice->D; d++) {
        dS -= this->forward[d]->value;
        dS -= this->backward[d]->value;
    }
    dS *= delta;
    dS += lattice->D * delta2;
//    dS *= lattice->beta;
    dS += lattice->mu * delta2 / 2;
    dS += lattice->lambda * delta4 / 2;
    
    // metropolis algorithm
    if (dS <= 0 || lattice->rand(0, 1) < exp(-dS)) {
        set_value(new_val);
        return 1;
    } else {
        return 0;
    }
}

phi4_lattice::phi4_lattice(int N, int T, int D, double beta, double mu, double lambda, bool cold) {
    this->N = N;
    this->T = T;
    this->D = D;
    this->beta = beta;
    this->mu = mu * mu;
    this->lambda = lambda;
    this->parallel = false;
    this->z = 0.5;

    random_device rd;  // seed for the random number engine
    gen = mt19937(rd()); // mersenne_twister_engine seeded with rd()
    
    n_sites = pow(N, D - 1) * T;
    n_slice = n_sites / T; // number of sites in each time slice
    site = new phi4_site[n_sites];
    
    for (int s = 0; s < n_sites; s++) {
        site[s].lattice = this;
        site[s].is_locked = NULL;
        site[s].reset_value(cold);
        
        for (int d = 0; d < D; d++) {
            update_i(s);
            move(d, +1);
            site[s].forward[d] = &site[update_s()];
            move(d, -2);
            site[s].backward[d] = &site[update_s()];
        }
    }
}

phi4_lattice::phi4_lattice(phi4_lattice* lattice) {
    this->N = lattice->N;
    this->T = lattice->T;
    this->D = lattice->D;
    this->beta = lattice->beta;
    this->mu = lattice->mu;
    this->lambda = lattice->lambda;
    this->parallel = lattice->parallel;
    this->z = 1.0;

    random_device rd;  // seed for the random number engine
    gen = mt19937(rd()); // mersenne_twister_engine seeded with rd()
    
    n_sites = pow(N, D - 1) * T;
    n_slice = n_sites / T; // number of sites in each time slice
    site = new phi4_site[n_sites];
    
    for (int s = 0; s < n_sites; s++) {
        site[s].lattice = this;
        site[s].is_locked = NULL;
        site[s].set_value(lattice->site[s].value);
        
        for (int d = 0; d < D; d++) {
            update_i(s);
            move(d, +1);
            site[s].forward[d] = &site[update_s()];
            move(d, -2);
            site[s].backward[d] = &site[update_s()];
        }
    }
}

phi4_lattice::~phi4_lattice() {
    delete[] site;
}

double phi4_lattice::rand(double min, double max) {
    return ((double)gen()) / 4294967295.0 * (max - min) + min;
}

void phi4_lattice::update_i(int s) {
    i[0] = s / (int)pow(N, D - 1);
    for (int d = 1; d < D; d++) i[d] = s / (int)pow(N, D - d - 1) % N;
}

int phi4_lattice::update_s() {
    int s = 0;
    for (int d = 0; d < D; d++) s += i[d] * (int)pow(N, D - d - 1);
    return s;
}

void phi4_lattice::move(int d, int n) {
    
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

double phi4_lattice::action() {
    // compute the average action of the lattice
    double U = 0;
    for (int s = 0; s < n_sites; s++) U += site[s].action();
    return U / n_sites;
}

double phi4_lattice::phi_ave() {
    // compute the average site value
    double phi = 0;
    for (int s = 0; s < n_sites; s++) phi += site[s].value;
    return phi / n_sites;
}

double async_sweep(phi4_site* site, int n) {
    double accept = 0.0;

    // start at a random spot in the time slice to mitigate parallel effects
    int start = rand() % n;
    int end = start + n;
    
    for (int s = start; s < end; s++) accept += site[s % n].sweep();

    return accept;
}

void phi4_lattice::sweep(int n_sweeps) {
    
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
    
    if (accept < 0.78) {
        z *= 1.05;
        if (z > 0.5) z = 0.5;
    } else if (accept > 0.82) {
        z *= 0.95;
    }
}

int phi4_lattice::thermalize(int n_min, int n_max) {
    
    // create a cold start model
    phi4_lattice coldStart = phi4_lattice(N, T, D, beta, mu, lambda, true);
    double hot_action = action();
    double cold_action = coldStart.action();
    int m;
    for (m = 0; hot_action > cold_action; m++) {
        // sweep both models until the average plaquette reaches the same value
        sweep();
        coldStart.sweep();
        hot_action = action();
        cold_action = coldStart.action();
//        cout << "dE = " << abs(cold_action - hot_action) << "\n";
        if (n_max && (m > n_max)) break;
    }
    
    // sweep a few more times
    int n_therm = max(n_min, m * 2);
    for (int i = m; i < n_therm; i++) sweep();
    
    cout << "Eave = " << action() / n_sites << "\n";
    
    // return value is total number of sweeps
    return n_therm;
}
