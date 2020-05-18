//
//  su3_p.hpp
//  Lattice
//
//  Created by Evan Owen on 11/16/18.
//  Copyright © 2018 Evan Owen. All rights reserved.
//

#ifndef su3_hpp
#define su3_hpp

#include <random>
#include <complex>
#include "Eigen/Dense"

class su3_lattice;
typedef Eigen::Matrix<std::complex<double>, 3, 3> su3_link;
typedef Eigen::Matrix<std::complex<double>, 1, 3> su3_vector;
#define su3_identity su3_link::Identity()
#define su3_zero su3_link::Zero()

typedef Eigen::Matrix<std::complex<double>, 2, 2> su2_link;
#define su2_identity su2_link::Identity()
#define su2_zero su2_link::Zero()

#define I complex<double>(0,1)
#define D_MAX 4

class su3_site {
public:
    // variables
    su3_lattice* lattice; // parent lattice
    su3_link link[D_MAX]; // link values in each direction
    su3_link link_inverse[D_MAX]; // link inverse values
    su3_site* forward[D_MAX]; // adjacent sites in the forward direction
    su3_site* backward[D_MAX]; // adjacent sites in the backward direction
    int thread_id; // thread assigned to this site
    bool is_locked;
    
    // methods
    void lock();
    void unlock();
    void reset_links(bool cold);
    void copy_links(su3_site* site);
    void set_link(int d, su3_link value);
    su2_link create_su2(bool random);
    su3_link create_link(bool random);
    double plaq();
    double wilson_loop(int a, int b);
    double polyakov_loop(int r);
    double correlator(int T);
    double four_point(int T, int R);
    std::complex<double> contract_3(su3_link A, su3_link B, su3_link C);
    su3_link contract_2(su3_link A, su3_link B);
    double contract_3x1_L();
    double contract_3x1();
    double contract_3x2_L();
    double contract_3x2();
    double contract_3x2x2();
    double contract_3x2x2_T();
    su3_link overrelax(su3_link g);
    su3_link make_unitary(su3_link g);
    void relax(bool coulomb);
    double sum_landau();
    double sum_coulomb();
    su3_link sum_G(bool coulomb);
    double sweep();
    double sweep_link(int d1);
};

class su3_lattice {
public:
    // variables
    int N; // array size (spacial)
    int T; // array size (time)
    int D; // number of dimensions
    double beta; // coupling factor
    su3_site* site; // site values
    int n_sites; // number of sites
    int i[D_MAX]; // current position
    std::mt19937 gen; // random number generator
    double z; // metropolis factor
    su2_link sigma[4]; // pauli spin matrices
    int n_threads; // number of parallel threads
    
    // methods
    su3_lattice(int N, int T, int D, double beta, bool cold = false);
    su3_lattice(su3_lattice* lattice);
    ~su3_lattice();
    double rand(double min = 0.0, double max = 1.0);
    void update_i(int s);
    int update_s();
    void move(int d, int n);
    double plaq();
    double wilson_loop(int a, int b);
    double polyakov_loop(int r);
    double correlator(int T);
    double four_point(int T, int R);
    double contract_3x1_L();
    double contract_3x1();
    double contract_3x2_L();
    double contract_3x2();
    double contract_3x2x2();
    double contract_3x2x2_T();
    void sweep(int n_sweeps = 0);
    int thermalize(int n_min = 0, int n_max = 0);
    double ave_link_t();
    double theta(bool coulomb);
    long double relax(long double error_target, bool coulomb);
    void write_four_point(const char* filename, bool coulomb);
    void write_correlator(const char* filename, bool coulomb);
};

#endif /* su3_hpp */
