//
//  ising.hpp
//
//  Created by Evan Owen on 4/3/18.
//  Copyright © 2018 Evan Owen. All rights reserved.
//

#ifndef ising_hpp
#define ising_hpp

#include <random>

class ising {
public:
    // variables
    int N; // array size
    int D; // number of dimensions
    double beta; // coupling factor
    double* site; // site values
    int n_sites; // number of sites
    int* i; // current position
    
    // random number generator stuff
    std::mt19937 gen;
    std::uniform_real_distribution<> rng;

    // methods
    ising(int N, int D, double beta, bool cold = false);
    ~ising();
    void init_hot(int s);
    void init_cold(int s);
    void log_site(int s);
    double get_site();
    void set_site(double value);
    void log();
    void move(int d, int n);
    double site_energy(int s);
    double energy();
    void update_i(int s);
    int update_s();
    int thermalize();
    void sweep();
    void sweep_site(int s);
    void iterate(void (ising::*f)(int));
    double iterate_sum(double (ising::*f)(int));
};

#endif /* ising_lattice_hpp */
