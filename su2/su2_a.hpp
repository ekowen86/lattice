//
//  su2_a.hpp
//
//  Created by Evan Owen on 4/9/18.
//  Copyright © 2018 Evan Owen. All rights reserved.
//

#ifndef su2_a_hpp
#define su2_a_hpp

#include <random>
#include "Eigen/Dense"

class su2_a;
typedef Eigen::Matrix<std::complex<double>, 2, 2> su2_link;
#define su2_identity su2_link::Identity()
#define su2_zero su2_link::Zero()
#define D_MAX 4

class su2_site {
public:
    // variables
    su2_a* lattice; // parent lattice
    su2_link link[D_MAX]; // link values in each direction
    su2_link link_inverse[D_MAX]; // link inverse values
    su2_site* forward[D_MAX]; // adjacent sites in the forward direction
    su2_site* backward[D_MAX]; // adjacent sites in the backward direction

    // methods
    void reset_links(bool cold);
    void copy_links(su2_site* site);
    void set_link(int d, su2_link value);
    su2_link create_link(double a, double b, double c);
    su2_link random_link();
    double plaq();
    double wilson_loop(int a, int b);
    double correlator(int T);
    double four_point(int T, int R);
    su2_link overrelax(su2_link g);
    void relax_landau();
    void relax_coulomb();
    double sum_landau();
    double sum_coulomb();
    double theta();
    void sweep();
    void sweep_link(int d1);
};

class su2_a {
public:
    // variables
    int N; // array size (spacial)
    int T; // array size (time)
    int D; // number of dimensions
    double beta; // coupling factor
    su2_site* site; // site values
    int n_sites; // number of sites
    int n_links; // number of links
    int n_plaquettes; // number of plaquettes
    int i[D_MAX]; // current position
    std::mt19937 gen; // random number generator

    // methods
    su2_a(int N, int T, int D, double beta, bool cold = false);
    su2_a(su2_a* lattice);
    ~su2_a();
    double rand(double min = 0.0, double max = 1.0);
    void update_i(int s);
    int update_s();
    void move(int d, int n);
    double plaq();
    double wilson_loop(int a, int b);
    double correlator(int T);
    double four_point(int T, int R);
    void sweep(int n_sweeps = 0);
    int thermalize(int n_min = 0);
    double ave_link_t();
    long double relax(bool coulomb);
    void write_four_point(const char* filename, bool coulomb);
    void write_correlator(const char* filename, bool coulomb);
};

#endif /* su2_a_hpp */
