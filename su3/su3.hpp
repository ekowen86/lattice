//
//  su3.hpp
//  Lattice
//
//  Created by Evan Owen on 10/19/18.
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
#define Q(s) s->link[0]
#define D_MAX 4
#define QUARK_MAX 2

class su3_site {
public:
    // variables
    su3_lattice* lattice; // parent lattice
    su3_link link[D_MAX]; // link values in each direction
    su3_link link_inverse[D_MAX]; // link inverse values
    su3_site* forward[D_MAX]; // adjacent sites in the forward direction
    su3_site* backward[D_MAX]; // adjacent sites in the backward direction
    bool is_locked;
    
    // methods
    bool lock();
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
    double mag_U();
    double abs_U();
    double four_point(int T, int R);
    double contract_24(su3_link* A);
    double contract_12(su3_link* A);
    su3_link contract_8(su3_link A, su3_link B, su3_link C, su3_link D, su3_link E, su3_link F, su3_link G, su3_link H);
    double contract_6(su3_link A, su3_link B, su3_link C, su3_link D, su3_link E, su3_link F);
    su3_link contract_4(su3_link A, su3_link B, su3_link C, su3_link D);
    double contract_3(su3_link A, su3_link B, su3_link C);
    su3_link contract_2(su3_link A, su3_link B);
    double contract_1(su3_link A, su3_link B);
    void contract_hadron(su3_link* A, int size, int n);
    double contract_baryon(su3_link* A, int size, int n);
    double separate_baryons(su3_link* A, int size, int n);
    double contract_rand(int D1, int D2, int D3, int n, int type);
    double contract_distance(int x1, int y1, int x2, int y2);
    void quark_distance(double* U);
    void quark_rand(double* U, int type);
    void quark_density(double* U, int type);
    void quark_0(double* U);
    void quark_1(double* U);
    void quark_2(double* U);
    void quark_3(double* U);
    void quark_4(double* U);
    void quark_5(double* U);
    void quark_special(double* U);
    void quark_gas_1(double* U, bool baryon);
    void quark_gas_2(double* U, bool baryon);
    void quark_gas_3(double* U, bool baryon);
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
    bool parallel;
    su3_site* site; // site values
    int n_sites; // number of sites
    int n_slice; // number of sites in each time slice
    int i[D_MAX]; // current position
    std::mt19937 gen; // random number generator
    double z; // metropolis factor
    su2_link sigma[4]; // pauli spin matrices
    
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
    void quark_distance(double* U);
    void quark_rand(double* U);
    void quark(double* U);
    void sweep(int n_sweeps = 0);
    int thermalize(int n_min = 0, int n_max = 0);
    double ave_link_t();
    double theta(bool coulomb);
    long double relax(long double error_target, bool coulomb);
    void write_four_point(const char* filename, bool coulomb);
    void write_correlator(const char* filename, bool coulomb);
};

#endif /* su3_hpp */
