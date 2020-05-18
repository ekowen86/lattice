//
//  qed.hpp
//  Lattice
//
//  Created by Evan Owen on 2/16/20.
//  Copyright © 2020 Evan Owen. All rights reserved.
//

#ifndef qed_hpp
#define qed_hpp

#include <random>
#include <complex>
#include "Eigen/Dense"

using namespace std;

class qed_lattice;
class qed_site;

typedef complex<double> qed_link;
typedef complex<double> qed_vector;
#define qed_identity qed_link(1, 0)
#define qed_zero qed_vector(0, 0)

#define I qed_vector(0,1)
#define D_MAX 4 // max number of dimensions
#define F_MAX 4 // max number of fermion flavors

class qed_psi {
    
public:
    qed_vector v1;
    qed_vector v2;
    
    qed_psi() {
        
    }
    
    qed_psi(qed_vector v1, qed_vector v2) {
        this->v1 = v1;
        this->v2 = v2;
    }
    
    qed_vector operator* (qed_psi g) {
        return (v1 * g.v2 - g.v1 * v2) / 2.0;
    }

    qed_psi operator- () {
        return qed_psi(-v1, -v2);
    }

    qed_psi operator* (complex<double> d) {
        return qed_psi(v1 * d, v2 * d);
    }

    qed_psi conj() {
        return qed_psi(std::conj(v2), std::conj(v1));
    }
};

class qed_spinor {
public:
    // variables
    qed_psi psi[4]; // fermion spinor values
    qed_psi psi_bar[4]; // inverse fermion spinor values

    // methods
    qed_spinor();
    qed_spinor(qed_psi* psi);
    void set_psi(qed_psi* psi);
//    qed_spinor operator+ (qed_spinor spinor); // add two spinors
    qed_spinor operator* (qed_link u); // multiplication by a U(1) element
    qed_spinor operator* (int gamma); // multiplication by a gamma matrix
    qed_link operator* (qed_spinor spinor); // contraction with another spinor
};

class qed_site {
public:
    // variables
    qed_lattice* lattice; // parent lattice
    qed_site* forward[D_MAX]; // adjacent sites in the forward direction
    qed_site* backward[D_MAX]; // adjacent sites in the backward direction
    qed_link link[D_MAX]; // link values in each direction
    qed_link link_inverse[D_MAX]; // link inverse values
    qed_spinor spinor[F_MAX]; // fermion spinor (psi)
    bool is_locked;
    qed_link forward_edge[D_MAX]; // -1 if site is at the forward edge of the lattice, 1 otherwise
    qed_link backward_edge[D_MAX]; // -1 if site is at the backward edge of the lattice, 1 otherwise

    // methods
    bool lock();
    void unlock();
    void reset_site(bool cold);
    void copy_site(qed_site* site);
    qed_link create_link(bool random, double z = 1.0);
//    qed_vector create_vector(bool random, double z = 1.0);
    qed_spinor create_spinor(bool random, double z = 1.0);
    void set_link(int d, qed_link value);
    qed_link action();
    qed_link gauge_action();
    qed_link fermion_action(int f);
    double sweep();
    double sweep_link(int d1);
    double sweep_fermion(int f);
};

class qed_lattice {
public:
    // variables
    int N; // array size (spacial)
    int T; // array size (time)
    int D; // number of dimensions
    int F; // number of fermion flavors
    double a; // lattice spacing
    double beta; // gauge coupling
    double m[F_MAX]; // fermion masses
    bool parallel;
    qed_site* site; // site values
    int n_sites; // number of sites
    int n_slice; // number of sites in each time slice
    int i[D_MAX]; // current position
    std::mt19937 gen; // random number generator
    double z_G; // metropolis factor for gauge action
    double z_F[F_MAX]; // metropolis factor for fermion action

    // methods
    qed_lattice(int N, int T, int D, int F, double a, double beta, double* m, bool cold = false);
    qed_lattice(qed_lattice* lattice);
    ~qed_lattice();
    double rand(double min = 0.0, double max = 1.0);
    void update_i(int s);
    int update_s();
    void move(int d, int n);
    double action();
    void sweep(int n_sweeps = 0);
    int thermalize(int n_min = 0, int n_max = 0);
};

#endif /* qed_hpp */
