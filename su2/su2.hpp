//
//  su2.hpp
//
//  Created by Evan Owen on 4/9/18.
//  Copyright © 2018 Evan Owen. All rights reserved.
//

#ifndef su2_hpp
#define su2_hpp

#include <random>
#include "Eigen/Dense"

class su2;
typedef Eigen::Matrix<std::complex<double>, 2, 2> su2_link;
#define su2_identity su2_link::Identity()
#define su2_zero su2_link::Zero()
#define D_max 4

class su2_site {
public:
    su2* lattice;
    su2_link link[D_max];
    su2_link link_inverse[D_max];
    
    su2_site* forward[D_max];
    su2_site* backward[D_max];
    
    void reset_links(bool cold);
    void copy_links(su2_site* site);
    void set_link(int d, su2_link value);
    su2_link create_link(double a, double b, double c);
    double energy();
    double correlator_landau(int T);
    double correlator_coulomb(int T);
    double loop_energy(int a, int b);
    double four_point_landau(int T, int R);
    double four_point_coulomb(int T, int R);
    su2_link overrelax(su2_link g);
    void relax_landau();
    void relax_coulomb();
    double sum_landau();
    double sum_coulomb();
    void sweep();
    void sweep_link(int d1);
};

class su2 {
public:
    // variables
    int N; // array size
    int D; // number of dimensions
    double beta; // coupling factor
    su2_site* site; // site values
    int n_sites; // number of sites
    int n_links; // number of links
    int n_plaquettes; // number of plaquettes
    int i[D_max]; // current position
    std::mt19937 gen;

    // methods
    su2(int N, int D, double beta, bool cold = false);
    su2(su2* lattice);
    ~su2();
    double rand(double min = 0.0, double max = 1.0);
    su2_link get_link(int d);
    su2_link get_link_inverse(int d);
    void set_link(int d, su2_link value);
    void move(int d, int n);
    double energy();
    double loop_energy(int a, int b);
    double wilson_energy(int a, int b);
    double correlator(int T, bool coulomb);
    double four_point(int T, int R, bool coulomb);
    void update_i(int s);
    int update_s();
    int thermalize(int n_min = 0);
    void sweep(int n_sweeps = 0);
    double ave_link_t();
    double relax(bool coulomb);
    void write_correlator(const char* filename, int n_prop, bool coulomb);
    void write_four_point(const char* filename, int n_prop, bool coulomb);
};

#endif /* su2_hpp */
