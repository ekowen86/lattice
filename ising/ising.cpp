//
//  ising.cpp
//
//  Created by Evan Owen on 4/3/18.
//  Copyright © 2018 Evan Owen. All rights reserved.
//

#include <string.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <random>
#include "ising.hpp"

ising::ising(int N, int D, double beta, bool cold) {
    this->N = N;
    this->D = D;
    this->beta = beta;
    
    std::random_device rd;  // seed for the random number engine
    gen = std::mt19937(rd()); // mersenne_twister_engine seeded with rd()
    rng = std::uniform_real_distribution<>(0.0, 1.0);

    i = new int[D];
    n_sites = pow(N, D);
    site = new double[n_sites];
    
    if (cold) {
        iterate(&ising::init_cold);
    } else {
        iterate(&ising::init_hot);
    }
}

ising::~ising() {
    delete i;
    delete site;
}

void ising::init_hot(int s) {
    if (rng(gen) > 0.5) {
        site[s] = +1;
    } else {
        site[s] = -1;
    }
}

void ising::init_cold(int s) {
    site[s] = +1;
}

void ising::log_site(int s) {

    update_i(s);
    for (int d = 0; d < D; d++) {
        std::cout << i[d];
    }

    std::cout << ": " << site[s] << "\n";
}

void ising::log() {
    iterate(&ising::log_site);
}

void ising::update_i(int s) {
    for (int d = 0; d < D; d++) {
        i[d] = s / (int)pow(N, D - d - 1) % N;
    }
}

int ising::update_s() {
    int s = 0;
    for (int d = 0; d < D; d++) {
        s += i[d] * (int)pow(N, D - d - 1);
    }
    return s;
}

double ising::get_site() {
    int s = update_s();
    return site[s];
}

void ising::set_site(double value) {
    int s = update_s();
    site[s] = value;
}

void ising::iterate(void (ising::*f)(int)) {
    
    memset(i, 0, D);
    for (int s = 0; s < n_sites; s++) {
        (this->*f)(s);
    }
}

double ising::iterate_sum(double (ising::*f)(int)) {
    double sum = 0;
    
    memset(i, 0, D);
    for (int s = 0; s < n_sites; s++) {
        sum += (this->*f)(s);
    }
    return sum;
}

int ising::thermalize() {
    // create a cold start model
    ising coldStart = ising(N, D, beta, true);
    double E = energy();
    double coldE = coldStart.energy();
    int m = 0;
    for (m = 0; E < coldE; m++) {
        // sweep both models until the energies reach the same value
        E = energy();
        coldE = coldStart.energy();
        sweep();
        coldStart.sweep();
    }
    
    // sweep a few more times
    int i = 0;
    for (i = 0; i < (m / 2); i++) {
        sweep();
    }
    
    // return value is total number of sweeps
    return m + i;
}

void ising::sweep() {
    iterate(&ising::sweep_site);
}

void ising::sweep_site(int s) {
    
    double dH = 0;
    update_i(s);
    for (int d = 0; d < D; d++) {
        move(d, 1);
        dH += get_site();
        move(d, -2);
        dH += get_site();
        move(d, 1);
    }
    dH *= 2 * get_site();

    if (dH >= 0 || rng(gen) < exp(beta * dH)) {
        // flip the site
        if (get_site() == +1) {
            set_site(-1);
        } else {
            set_site(+1);
        }
    }
}

void ising::move(int d, int n) {
    if (n > 0) {
        i[d] = (i[d] + n) % N;
    } else {
        i[d] = (i[d] + n + N) % N;
    }
}

double ising::site_energy(int s) {
    update_i(s);
    double u = 0;
    for (int d = 0; d < D; d++) {
        move(d, 1);
        u += get_site();
        move(d, -1);
    }
    return u * get_site();
}

double ising::energy() {
    return iterate_sum(&ising::site_energy) / n_sites;
}

