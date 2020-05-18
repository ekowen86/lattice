//
//  su2_a.cpp
//
//  Created by Evan Owen on 4/9/18.
//  Copyright © 2018 Evan Owen. All rights reserved.
//

#include <iostream>
#include <iomanip>
#include <fstream>
#include <random>
#include <cmath>
#include "Eigen/Dense"
#include "su2_a.hpp"

using namespace std;

void su2_site::reset_links(bool cold) {
    su2_link link = su2_identity;
    for (int d = 0; d < lattice->D; d++) {
        if (!cold) link = random_link();
        set_link(d, link);
    }
}

void su2_site::copy_links(su2_site* site) {
    for (int d = 0; d < lattice->D; d++) set_link(d, su2_link(site->link[d]));
}

void su2_site::set_link(int d, su2_link value) {
    value /= sqrt(value.determinant().real()); // make sure the determinant is 1
    link[d] = value;
    link_inverse[d] = value.inverse();
}

su2_link su2_site::create_link(double a, double b, double c) {
    double g = sqrt(1 - a * a);
    double theta = acos(2 * b - 1);
    double phi = c * 2 * M_PI;
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

su2_link su2_site::random_link() {
    return create_link(1 - 2 * lattice->rand(), lattice->rand(), lattice->rand());
}

double su2_site::plaq() {
    // compute the average plaquette at this site
    double U = 0;
    su2_link u;
    for (int d1 = 0; d1 < lattice->D; d1++) {
        for (int d2 = d1 + 1; d2 < lattice->D; d2++) {
            u = link[d1];
            u *= forward[d1]->link[d2];
            u *= forward[d2]->link_inverse[d1];
            u *= link_inverse[d2];
            U += 1 - u.trace().real() / 2;
        }
    }
    return U / lattice->D / (lattice->D - 1);
}

double su2_site::wilson_loop(int a, int b) {
    // compute the average a x b wilson loop at this site
    double U = 0;
    su2_link u;
    su2_site* s = this;
    int x;
    for (int d1 = 0; d1 < lattice->D; d1++) {
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
            U += u.trace().real() / 2;
        }
    }
    return U / lattice->D / (lattice->D - 1);
}

double su2_site::correlator(int T) {
    // correlator of duration T
    su2_link u = su2_identity;
    su2_site* s = this;
    int x;
    for (x = 0; x < T; x++) {
        u *= s->link[0];
        s = s->forward[0];  // time direction is always 0
    }
    return u.trace().real() / 2;
}

double su2_site::four_point(int T, int R) {
    // four-point function of duration T and spacing R
    double U = 0;
    su2_link u;
    su2_site* s = this;
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
        U += u.trace().real() / 2;
    }
    return U / (lattice->D - 1);
}

su2_link su2_site::overrelax(su2_link g) {
    
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
    return gw / sqrt(gw.determinant().real());
}

void su2_site::relax_landau() {
    su2_link G = su2_zero;
    
    // sum all links connected to this site
    for (int d = 0; d < lattice->D; d++) {
        G += link[d];
        G += backward[d]->link_inverse[d];
    }
    
    su2_link ginv = G / sqrt(G.determinant().real());
    su2_link g = ginv.inverse();
    g = overrelax(g);
    ginv = g.inverse();
    
    // apply the gauge transformation
    su2_link u;
    for (int d = 0; d < lattice->D; d++) {
        u = link[d];
        set_link(d, g * u);
        u = backward[d]->link[d];
        backward[d]->set_link(d, u * ginv);
    }
}

void su2_site::relax_coulomb() {
    su2_link G = su2_zero;

    // sum all links connected to this site excluding the time direction 0
    for (int d = 1; d < lattice->D; d++) {
        G += link[d];
        G += backward[d]->link_inverse[d];
    }

    su2_link ginv = G / sqrt(G.determinant().real());
    su2_link g = ginv.inverse();
    g = overrelax(g);
    ginv = g.inverse();

    // apply the gauge transformation (including the time direction)
    su2_link u;
    for (int d = 0; d < lattice->D; d++) {
        u = link[d];
        set_link(d, g * u);
        u = backward[d]->link[d];
        backward[d]->set_link(d, u * ginv);
    }
}

double su2_site::sum_landau() {
    su2_link U = su2_zero;
    for (int d = 0; d < lattice->D; d++) U += link[d];
    return U.trace().real() / 2 / lattice->D;
}

double su2_site::sum_coulomb() {
    su2_link U = su2_zero;
    for (int d = 1; d < lattice->D; d++) U += link[d];
    return U.trace().real() / 2 / (lattice->D - 1);
}

void su2_site::sweep() {
    for (int d = 0; d < lattice->D; d++) sweep_link(d);
}

void su2_site::sweep_link(int d1) {
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
    
//    // metropolis algorithm
//    su2_link old_link = link[d1];
//    su2_link new_link = create_link(1 - 2 * lattice->rand(), lattice->rand(), lattice->rand());
//    double E1 = (old_link * U).trace().real() / 2;
//    double E2 = (new_link * U).trace().real() / 2;
//    double dE = E2 - E1;
//    if (dE >= 0 || lattice->rand() < exp(lattice->beta * dE)) {
//        set_link(d1, new_link);
//    }
    
    // heat bath algorithm
    double k = sqrt(U.determinant().real());
    double x, a0, g;
    
    while (true) {
        x = lattice->rand(exp(-2 * k * lattice->beta), 1.0);
        a0 = 1 + log(x) / k / lattice->beta;
        g = sqrt(1 - a0 * a0);
        if (lattice->rand() < g) break;
    }
    su2_link new_link = create_link(a0, lattice->rand(), lattice->rand());
    set_link(d1, new_link * (U / k).inverse());
}

su2_a::su2_a(int N, int T, int D, double beta, bool cold) {
    this->N = N;
    this->T = T;
    this->D = D;
    this->beta = beta;
    
    random_device rd;  // seed for the random number engine
    gen = mt19937(rd()); // mersenne_twister_engine seeded with rd()

    n_sites = pow(N, D - 1) * T;
    n_links = n_sites * D;
    n_plaquettes = n_links * (D - 1) / 2;
    site = new su2_site[n_sites];
    
    for (int s = 0; s < n_sites; s++) {
        site[s].lattice = this;
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

su2_a::su2_a(su2_a* lattice) {
    this->N = lattice->N;
    this->T = lattice->T;
    this->D = lattice->D;
    this->beta = lattice->beta;
    
    random_device rd;  // seed for the random number engine
    gen = mt19937(rd()); // mersenne_twister_engine seeded with rd()

    n_sites = pow(N, D - 1) * T;
    n_links = n_sites * D;
    n_plaquettes = n_links * (D - 1) / 2;
    site = new su2_site[n_sites];
    
    for (int s = 0; s < n_sites; s++) {
        site[s].lattice = this;
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

su2_a::~su2_a() {
    delete[] site;
}

double su2_a::rand(double min, double max) {
    return ((double)gen()) / 4294967295.0 * (max - min) + min;
}

void su2_a::update_i(int s) {
    i[0] = s / (int)pow(N, D - 1);
    for (int d = 1; d < D; d++) {
        i[d] = s / (int)pow(N, D - d - 1) % N;
    }
}

int su2_a::update_s() {
    int s = 0;
    for (int d = 0; d < D; d++) {
        s += i[d] * (int)pow(N, D - d - 1);
    }
    return s;
}

void su2_a::move(int d, int n) {
    
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

double su2_a::plaq() {
    // compute the average plaquette of the lattice
    double U = 0;
    for (int s = 0; s < n_sites; s++) U += site[s].plaq();
    return U / n_sites;
}

double su2_a::wilson_loop(int a, int b) {
    // compute the average a x b wilson loop of the lattice
    double sum = 0;
    for (int s = 0; s < n_sites; s++) sum += site[s].wilson_loop(a, b);
    return sum / n_sites;
}

double su2_a::correlator(int T) {
    // compute the average correlator of length T
    double sum = 0;
    for (int s = 0; s < n_sites; s++) sum += site[s].correlator(T);
    return sum / n_sites;
}

double su2_a::four_point(int T, int R) {
    // compute the average T x R 4-point function of the lattice
    double sum = 0;
    for (int s = 0; s < n_sites; s++) sum += site[s].four_point(T, R);
    return sum / n_sites;
}

void su2_a::sweep(int n_sweeps) {
    
    if (n_sweeps) {
        for (int i = 0; i < n_sweeps; i++) sweep(0);
        return;
    }
    
    for (int s = 0; s < n_sites; s++) site[s].sweep();
}

int su2_a::thermalize(int n_min) {
    
    // create a cold start model
    su2_a coldStart = su2_a(N, T, D, beta, true);
    double E = plaq();
    double Ec = coldStart.plaq();
    int m;
    for (m = 0; E > Ec; m++) {
        // sweep both models until the energies reach the same value
        sweep();
        coldStart.sweep();
        E = plaq();
        Ec = coldStart.plaq();
        cout << "dE = " << abs(Ec - E) << "\n";
    }
    
    // sweep a few more times
    int n_therm = max(n_min, m * 2);
    for (int i = m; i < n_therm; i++) {
        sweep();
    }
    
    // return value is total number of sweeps
    return n_therm;
}

double su2_a::ave_link_t() {
    // average trace of links in the t direction
    // this should be ~0 in Coulomb gauge
    
    double U = 0;
    
    for (int s = 0; s < n_sites; s++) U += site[s].link[0].trace().real() / 2;

    return U / n_sites;
}

long double su2_a::relax(bool coulomb) {
    
    // choose the relaxation function
    void (su2_site::*relax_site)() = &su2_site::relax_landau;
    if (coulomb) relax_site = &su2_site::relax_coulomb;

    double (su2_site::*sum_site)() = &su2_site::sum_landau;
    if (coulomb) sum_site = &su2_site::sum_coulomb;

    // iterate over each site
    for (int s = 0; s < n_sites; s++) (&site[s]->*relax_site)();
    
    // add up all the links
    long double U = 0;
    for (int s = 0; s < n_sites; s++) U += (long double)(&site[s]->*sum_site)();
    
    return U / n_sites;
}

void su2_a::write_four_point(const char* filename, bool coulomb) {

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

void su2_a::write_correlator(const char* filename, bool coulomb) {
    
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
