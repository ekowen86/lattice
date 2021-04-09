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

// MARK: phi4_site

void phi4_site::init(phi4_lattice* lattice, phi4_site* lattice_sites, int s) {
    this->lattice = lattice;
    
    int t = s / lattice->n_slice;
    gen = &lattice->gen[t];
    reset_value(true);
    is_locked = false;
    is_clustered = false;
    
    for (int d = 0; d < lattice->D; d++) {
        forward[d] = &lattice_sites[lattice->get_site(s, d, +1)];
        backward[d] = &lattice_sites[lattice->get_site(s, d, -1)];
    }
}

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

double phi4_site::rand_double(double min, double max) {
    std::uniform_real_distribution<double> dist(min, max);
    return dist(*gen);
}

int phi4_site::rand_int(int min, int max) {
    std::uniform_int_distribution<int> dist(min, max);
    return dist(*gen);
}

double phi4_site::rand_normal(double mean, double stdev) {
    std::normal_distribution<double> dist(mean, stdev);
    return dist(*gen);
}

void phi4_site::reset_value(bool cold) {
    if (cold) {
        set_value(0.0);
    } else {
        set_value(rand_double(-2.0, 2.0));
    }
}

void phi4_site::copy_value(phi4_site* site) {
    set_value(site->value);
}

void phi4_site::set_value(double value) {
    this->value = value;
    value2 = value * value;
    value4 = value2 * value2;
}

double phi4_site::action() {
    // compute the contribution to the lattice action at this site
    double S = 0;
    double s = 0;
    for (int d = 0; d < lattice->D; d++) {
        s = this->forward[d]->value - value;
        S += s * s / 2.0;
    }
    S += lattice->mu2 * value2 / 2.0;
    S += lattice->lambda * value4 / 2.0;
    return S;
}

double phi4_site::heat_bath() {
    double accept = 0;
    lock();
    for (int d = 0; d < lattice->D; d++) accept += heat_bath_link(d);
    unlock();
    return accept / lattice->D;
}

double phi4_site::heat_bath_link(int d1) {
    // choose a new value
    double new_val = value + lattice->rand_double(-1, 1) * lattice->z;
    double new_val2 = new_val * new_val;
    double new_val4 = new_val2 * new_val2;

    double delta = new_val - value;
    double delta2 = new_val2 - value2;
    double delta4 = new_val4 - value4;

    double dS = 0;
    for (int d = 0; d < lattice->D; d++) {
        dS -= this->forward[d]->value;
        dS -= this->backward[d]->value;
    }
    dS *= delta;
    dS += lattice->D * delta2;
    dS += lattice->mu2 * delta2 / 2;
    dS += lattice->lambda * delta4 / 2;
    
    // metropolis algorithm
    if (dS <= 0 || lattice->rand_double(0, 1) < exp(-dS)) {
        set_value(new_val);
        return 1;
    } else {
        return 0;
    }
}

bool phi4_site::wolff_test(double test_value) {
    // fail if already part of the cluster
    if (is_clustered) return false;
    
    // fail if signs don't match
    if (signbit(test_value) != signbit(value)) return false;
    
    if (rand_double() < (1 - exp(-2 * test_value * value))) {
        is_clustered = true;
    }
    return is_clustered;
}

void phi4_site::wolff_add_neighbors() {
    for (int d = 0; d < lattice->D; d++) {
        if (forward[d]->wolff_test(value)) forward[d]->wolff_add_neighbors();
        if (backward[d]->wolff_test(value)) backward[d]->wolff_add_neighbors();
    }
}

void phi4_site::wolff_flip() {
    if (!is_clustered) return;
    set_value(-value);
    is_clustered = false;
}

// MARK: phi4_lattice

phi4_lattice::phi4_lattice(int L, int T, int D, double mu2, double lambda, bool cold) {
    this->L = L;
    this->T = T;
    this->D = D;
    this->mu2 = mu2;
    this->lambda = lambda;
    this->parallel = false;
    this->z = 0.5;
    this->n_wolff = 5;

    init();

    for (int s = 0; s < n_sites; s++) site[s].reset_value(cold);
}

phi4_lattice::phi4_lattice(phi4_lattice* lattice) {
    this->L = lattice->L;
    this->T = lattice->T;
    this->D = lattice->D;
    this->mu2 = lattice->mu2;
    this->lambda = lattice->lambda;
    this->parallel = lattice->parallel;
    this->z = lattice->z;
    this->n_wolff = lattice->n_wolff;

    init();
    
    for (int s = 0; s < n_sites; s++) site[s].copy_value(&lattice->site[s]);
}

phi4_lattice::~phi4_lattice() {
}

void phi4_lattice::init() {
    
    // initialize rng
    random_device rd;  // seed for the random number engine
    for (int t = 0; t < (T + 1); t++) gen.push_back(mt19937(rd())); // mersenne_twister_engine seeded with rd()

    // initialize sites
    n_sites = pow(L, D - 1) * T;
    n_slice = n_sites / T;
    site.resize(n_sites);

    for (int s = 0; s < n_sites; s++) {
        site[s].init(this, site.data(), s);
    }
}

double phi4_lattice::rand_double(double min, double max) {
    std::uniform_real_distribution<double> dist(min, max);
    return dist(gen[T]);
}

int phi4_lattice::rand_int(int min, int max) {
    std::uniform_int_distribution<int> dist(min, max);
    return dist(gen[T]);
}

double phi4_lattice::rand_normal(double mean, double stdev) {
    std::normal_distribution<double> dist(mean, stdev);
    return dist(gen[T]);
}

int phi4_lattice::get_site(int s, int d, int n) {
    
    // get the site's coordinates
    int x[D];
    x[0] = s / n_slice;
    for (int d1 = 1; d1 < D; d1++) x[d1] = (s / (int)pow(L, D - d1 - 1)) % L;

    // go n steps in the d direction
    if (d == 0) {
        // time direction
        x[0] += n;
        while (x[0] < 0) x[0] += T;
        x[0] %= T;
//        if (n > 0) {
//            x[0] = (x[0] + n) % T;
//        } else {
//            x[0] = (x[0] + n + T) % T;
//        }
    } else {
        // spatial direction
        x[d] += n;
        while (x[d] < 0) x[d] += L;
        x[d] %= L;
//        if (n > 0) {
//            x[d] = (x[d] + n) % L;
//        } else {
//            x[d] = (x[d] + n + L) % L;
//        }
    }

    // get the index of the new site
    s = 0;
    for (int d1 = 0; d1 < D; d1++) s += x[d1] * (int)pow(L, D - d1 - 1);
    return s;
}

//void phi4_lattice::update_i(int s) {
//    i[0] = s / (int)pow(N, D - 1);
//    for (int d = 1; d < D; d++) i[d] = s / (int)pow(N, D - d - 1) % N;
//}
//
//int phi4_lattice::update_s() {
//    int s = 0;
//    for (int d = 0; d < D; d++) s += i[d] * (int)pow(N, D - d - 1);
//    return s;
//}
//
//void phi4_lattice::move(int d, int n) {
//    
//    if (d == 0) {
//        if (n > 0) {
//            i[0] = (i[0] + n) % T;
//        } else {
//            i[0] = (i[0] + n + T) % T;
//        }
//    } else {
//        if (n > 0) {
//            i[d] = (i[d] + n) % N;
//        } else {
//            i[d] = (i[d] + n + N) % N;
//        }
//    }
//}

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

double async_heat_bath(phi4_site* site, int n) {
    double accept = 0.0;

    // start at a random spot in the time slice to mitigate parallel effects
    int start = rand() % n;
    int end = start + n;
    
    for (int s = start; s < end; s++) accept += site[s % n].heat_bath();

    return accept;
}

void phi4_lattice::heat_bath(int n_sweeps) {
    
    if (n_sweeps) {
        for (int i = 0; i < n_sweeps; i++) heat_bath(0);
        return;
    }
    
    double accept = 0;
    
    if (parallel) {
        future<double> sweep_slice[T];
        for (int t = 0; t < T; t++) sweep_slice[t] = async(async_heat_bath, &site[t * n_slice], n_slice);
        for (int t = 0; t < T; t++) accept += sweep_slice[t].get();
    } else {
        for (int s = 0; s < n_sites; s++) accept += site[s].heat_bath();
    }
    accept /= n_sites;
    
    if (accept < 0.78) {
        z *= 0.95;
//        if (z > 0.5) z = 0.5;
    } else if (accept > 0.82) {
        z *= 1.05;
    }
//    cout << "z = " << z << endl;
}

void phi4_lattice::wolff_update() {
    
    // do heat bath sweeps first
    heat_bath(n_wolff);
    
    // add one site to the cluster
    phi4_site seed = site[rand_int(0, n_sites - 1)];
    seed.is_clustered = true;
    
    // recursively add neighbors
    seed.wolff_add_neighbors();
    
    //
    for (int s = 0; s < n_sites; s++) site[s].wolff_flip();
}

int phi4_lattice::thermalize(int n_min, int n_max) {
    
    // create a cold start model
    phi4_lattice coldStart = phi4_lattice(L, T, D, mu2, lambda, true);
    double hot_action = action();
    double cold_action = coldStart.action();
    
    printf("hot_S = %6e\n", hot_action);
    printf("cold_S = %6e\n", cold_action);

    int m;
    for (m = 0; hot_action > cold_action; m++) {
        // sweep both models until the average plaquette reaches the same value
        heat_bath();
        coldStart.heat_bath();
        hot_action = action();
        cold_action = coldStart.action();
//        printf("dS = %6e\n", abs(cold_action - hot_action));
        if (n_max && (m > n_max)) break;
    }
    
    // sweep a few more times
    int n_therm = max(n_min, m * 2);
    for (int i = m; i < n_therm; i++) heat_bath();
    
    printf("S = %6e\n", action() / n_sites);
    
    // return value is total number of sweeps
    return n_therm;
}
