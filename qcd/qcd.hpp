//
//  qcd.hpp
//  Lattice
//
//  Created by Evan Owen on 2/28/20.
//  Copyright © 2020 Evan Owen. All rights reserved.
//

#ifndef qcd_hpp
#define qcd_hpp

#include <random>
#include <complex>
#include "Eigen/Dense"

using namespace std;

class qcd_lattice;
class qcd_site;

#define u1 complex<double>
#define I u1(0,1)

typedef Eigen::Matrix<u1, 2, 2> su2;
typedef Eigen::Matrix<u1, 2, 1> su2_vector;
#define su2_identity su2::Identity()
#define su2_zero su2::Zero()

typedef Eigen::Matrix<u1, 3, 3> su3;
typedef Eigen::Matrix<u1, 3, 1> su3_vector;
#define su3_identity su3::Identity()
#define su3_zero su3::Zero()

#define D_MAX 4 // max number of dimensions
#define F_MAX 4 // max number of fermion flavors

class qcd_spinor {
    
public:
    double component[24];

    qcd_spinor();
    qcd_spinor(qcd_spinor& psi);
    
    su3_vector psi(int s);
    su3_vector psi_conj(int s);
    u1 contract(qcd_spinor psi);
    u1 contract(su3 U, qcd_spinor psi);
    u1 contract(int gamma, qcd_spinor psi);
    u1 contract(int gamma, su3 U, qcd_spinor psi);
};

class qcd_site {
public:
    // variables
    qcd_lattice* lattice; // parent lattice
    qcd_site* forward[D_MAX]; // adjacent sites in the forward direction
    qcd_site* backward[D_MAX]; // adjacent sites in the backward direction
    su3 link[D_MAX]; // link values in each direction
    su3 link_inverse[D_MAX]; // link inverse values
    qcd_spinor psi[F_MAX]; // fermion spinor (psi)
    bool is_locked;
    u1 forward_edge[D_MAX]; // -1 if site is at the forward edge of the lattice, 1 otherwise
    u1 backward_edge[D_MAX]; // -1 if site is at the backward edge of the lattice, 1 otherwise

    // methods
    bool lock();
    void unlock();
    void reset_site(bool cold);
    void copy_site(qcd_site* site);
    double random_mag(double z = 1.0);
    su2 create_su2(bool random, double z = 1.0);
    su3 create_link(bool random, double z = 1.0);
    u1 create_u1(bool random);
    su3_vector random_vector();
    qcd_spinor create_spinor(bool random, double z = 1.0);
    void set_link(int d, su3 value);
    double action();
    double gauge_action();
    u1 psi_bar_psi(int f);
    u1 fermion_action(int f);
    double site_fermion_action(int f, qcd_spinor t_psi);
    double sweep_gauge();
    double sweep_link(int d1);
    double sweep_fermion(int f);
};

class qcd_lattice {
public:
    // variables
    int N; // array size (spacial)
    int T; // array size (time)
    int D; // number of dimensions
    int F; // number of fermion flavors
//    double a; // lattice spacing
    double beta; // gauge coupling
    double kappa[F_MAX]; // hopping factor
    double m[F_MAX]; // fermion mass fractions
    bool parallel;
    qcd_site* site; // site values
    int n_sites; // number of sites
    int n_slice; // number of sites in each time slice
    int i[D_MAX]; // current position
    std::mt19937 gen; // random number generator
    double z_G; // metropolis factor for gauge action
    double z_F[F_MAX]; // metropolis factor for fermion action

    // methods
    qcd_lattice(int N, int T, int D, int F, double kappa, double beta, double* m, bool cold = false);
    qcd_lattice(qcd_lattice* lattice);
    ~qcd_lattice();
    double rand(double min = 0.0, double max = 1.0);
    void update_i(int s);
    int update_s();
    void move(int d, int n);
    double plaquette();
    double action();
    double gauge_action();
    double psi_bar_psi(int f);
    double fermion_action();
    double xi_avg();
    void sweep(int n_sweeps = 1, bool gauge_only = false);
    void sweep_gauge(int n_sweeps = 1);
    double sweep_f(int n_sweeps = 1);
    int thermalize(int n_min = 0, int n_max = 0);
};

#endif /* qcd_hpp */
