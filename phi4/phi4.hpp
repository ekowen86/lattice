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
    double value2;
    double value4;
    phi4_site* forward[D_MAX]; // adjacent sites in the forward direction
    phi4_site* backward[D_MAX]; // adjacent sites in the backward direction
    std::mt19937* gen; // random number generator for this time-slice
    bool is_locked;
    bool is_clustered;
    
    // methods
    void init(phi4_lattice* lattice, phi4_site* lattice_sites, int s);
    bool lock();
    void unlock();
    double rand_double(double min = 0.0, double max = 1.0);
    int rand_int(int min, int max);
    double rand_normal(double mean = 0.0, double stdev = 1.0);
    void reset_value(bool cold);
    void copy_value(phi4_site* site);
    void set_value(double value);
    double action();
    double heat_bath();
    double heat_bath_link(int d1);
    bool wolff_test(double test_value);
    void wolff_add_neighbors();
    void wolff_flip();
};

class phi4_lattice {
public:
    // variables
    int L; // array size (spacial)
    int T; // array size (time)
    int D; // number of dimensions
    double mu2; // mass squared
    double lambda; // coupling
    bool parallel;
    std::vector<phi4_site> site; // site values
    int n_sites; // number of sites
    int n_slice; // number of sites in each time slice
    int i[D_MAX]; // current position
    std::vector<std::mt19937> gen; // array of random number generators (one for each time-slice)
    double z; // metropolis factor
    int n_wolff; // number of heat bath sweeps per wolff flip
    
    // methods
    phi4_lattice(int N, int T, int D, double mu2, double lambda, bool cold = false);
    phi4_lattice(phi4_lattice* lattice);
    ~phi4_lattice();
    void init();
    double rand_double(double min = 0.0, double max = 0.0);
    int rand_int(int min, int max);
    double rand_normal(double mean = 0.0, double stdev = 1.0);
    int get_site(int s, int d, int n);
    void move(int d, int n);
    double action();
    double phi_ave();
    void heat_bath(int n_sweeps = 0);
    void wolff_update();
    int thermalize(int n_min = 0, int n_max = 0);
};

#endif /* phi4_hpp */
