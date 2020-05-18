//
//  phi4.hpp
//
//  Created by Evan Owen on 11/7/19.
//  Copyright © 2019 Evan Owen. All rights reserved.
//

#ifndef phi4_hpp
#define phi4_hpp

#include <random>
#include <complex>

class phi4_lattice;

#define D_MAX 6

class phi4_site {
public:
    // variables
    phi4_lattice* lattice; // parent lattice
    double value;
//    phi4_link link[D_MAX]; // link values in each direction
//    phi4_link link_inverse[D_MAX]; // link inverse values
    phi4_site* forward[D_MAX]; // adjacent sites in the forward direction
    phi4_site* backward[D_MAX]; // adjacent sites in the backward direction
    bool is_locked;
    
    // methods
    bool lock();
    void unlock();
    void reset_value(bool cold);
    void set_value(double value);
    double action();
    double sweep();
    double sweep_link(int d1);
};

class phi4_lattice {
public:
    // variables
    int N; // array size (spacial)
    int T; // array size (time)
    int D; // number of dimensions
    double beta; // coupling factor
    double mu; // mass
    double lambda; // interation coefficient
    bool parallel;
    phi4_site* site; // site values
    int n_sites; // number of sites
    int n_slice; // number of sites in each time slice
    int i[D_MAX]; // current position
    std::mt19937 gen; // random number generator
    double z; // metropolis factor
    
    // methods
    phi4_lattice(int N, int T, int D, double beta, double mu, double lambda, bool cold = false);
    phi4_lattice(phi4_lattice* lattice);
    ~phi4_lattice();
    double rand(double min = 0.0, double max = 1.0);
    void update_i(int s);
    int update_s();
    void move(int d, int n);
    double action();
    double phi_ave();
    void sweep(int n_sweeps = 0);
    int thermalize(int n_min = 0, int n_max = 0);
};

#endif /* phi4_hpp */
