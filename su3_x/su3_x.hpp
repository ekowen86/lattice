//
//  su3_x.hpp
//  Lattice
//
//  Created by Evan Owen on 10/19/18.
//  Copyright © 2018 Evan Owen. All rights reserved.
//

#ifndef su3_x_hpp
#define su3_x_hpp

#include <random>
#include <complex>
#include <Eigen/Dense>

class su3_x_lattice;
typedef Eigen::Matrix<std::complex<double>, 3, 3> su3_link;
typedef Eigen::Matrix<std::complex<double>, 1, 3> su3_vector;
#define su3_identity su3_link::Identity()
#define su3_zero su3_link::Zero()

typedef Eigen::Matrix<std::complex<double>, 2, 2> su2_link;
#define su2_identity su2_link::Identity()
#define su2_zero su2_link::Zero()

#define I complex<double>(0,1)
#define Q(s) s->link[0]
#define D_MAX 6
#define SQRT3   1.73205080756887729352744634151  // sqrt(3)
#define SQRT1_3 0.577350269189625764509148780502 // 1 / sqrt(3)

class su3_x_site {
public:
    // variables
    su3_x_lattice* lattice; // parent lattice
    su3_link link[D_MAX]; // link values in each direction
    su3_link link_inverse[D_MAX]; // link inverse values
    su3_link p_link[D_MAX]; // conjugate momenta in each direction
    su3_x_site* forward[D_MAX]; // adjacent sites in the forward direction
    su3_x_site* backward[D_MAX]; // adjacent sites in the backward direction
    std::mt19937* gen; // random number generator for this time-slice
    bool forward_edge; // site is at the forward edge of the extra dimension
    bool backward_edge; // site is at the backward edge of the extra dimension
    bool is_locked;
    double eps;
    su3_link wf_z[D_MAX]; // exponent for wilson flow

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    
    // methods
    void init(su3_x_lattice* lattice, su3_x_site* lattice_sites, int s);
    su3_link make_unitary(const su3_link& g);
    su3_link cayley_ham(const su3_link& Q);
    bool lock();
    void unlock();
    double rand_double(double min = 0.0, double max = 1.0);
    int rand_int(int min, int max);
    double rand_normal(double mean = 0.0, double stdev = 1.0);
    void reset_links(bool cold);
    void copy_links(su3_x_site* site);
    void read_links(std::ifstream& ckptFile, bool bigEndian);
    void set_link(int d, const su3_link& value);
    void init_momenta();
    su2_link create_su2(bool random);
    su3_link create_link(bool random);
    double action();
    double hamiltonian();
    void hmc_step_p(double frac);
    void hmc_step_link();
    su3_link p_link_dot(int d);
    double link_trace();
    double plaq();
    su3_link plaquette(int d1, int d2);
    su3_link staple(int d1);
    su3_link reverse_staple(int d1);
    su3_link staple_x(int d1);
    su3_link reverse_staple_x(int d1);
    su3_link cloverleaf(int d1, int d2);
    double wilson_loop(int a, int b);
    std::complex<double> polyakov_loop(int r);
    double correlator(int T);
    double field_strength();
//    double field_strength_x();
    double topological_charge();
    double mag_U();
    double abs_U();
    double four_point(int T, int R);
    su3_link overrelax(su3_link g);
    void relax(bool coulomb);
    double sum_landau();
    double sum_coulomb();
    su3_link sum_G(bool coulomb);
    double heat_bath();
    double heat_bath_link(int d1);
    void cool();
    void cool_link(int d1);
    void wilson_flow(su3_x_site* target, int step);
    void wilson_flow_link(su3_x_site* target, int step, int d);
    void stout_smear(su3_x_site* target, double rho);
    void stout_smear_link(su3_x_site* target, double rho, int d1);
};

class su3_x_lattice {
public:
    // variables
    int N; // array size (spacial)
    int T; // array size (time)
    int N5; // 5th dimension size
    int D; // number of dimensions
    double beta; // coupling factor
    double eps5; // ratio of coupling in extra dimension to normal coupling
    int n_sites; // number of sites in the entire lattice
    int n_sites_5; // number of sites in a 5d sub-lattice
    int n_slice; // number of sites in each time slice of the entire lattice
    int n5_center; // center lattice
    std::vector<su3_x_site> site; // site values
    std::vector<std::mt19937> gen; // array of random number generators (one for each time-slice)
    std::vector<su2_link> sigma; // pauli spin matrices
    std::vector<su3_link> lambda; // gell-mann matrices
    double z; // heat bath metropolis factor
    int verbose;
    bool parallel;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    // hmc
    std::vector<su3_x_site> site_1; // site values for hmc
    double dt; // hmc step size
    int n_steps; // hmc step count
    int hmc_accept; // number of accepted hmc configurations
    int hmc_count; // number of attempted hmc configurations

    // time step for wilson flow algorithm
    double wf_dt;
    
    // methods
    su3_x_lattice(int N, int T, int N5, int D, double beta, double eps5, bool cold = false);
    su3_x_lattice(su3_x_lattice* lattice);
    su3_x_lattice(int N, int T, int N5, int D, double beta, double eps5, std::ifstream& ckptFile, bool isNersc = false);
    ~su3_x_lattice();
    void init();
    double rand_double(double min = 0.0, double max = 1.0);
    int rand_int(int min, int max);
    double rand_normal(double mean = 0.0, double stdev = 1.0);
    int get_site(int s, int d, int n);
    double link_trace();
    double plaq(int n5);
    double action(int n5);
    double wilson_loop(int a, int b, int n5);
    double polyakov_loop(int R, int n5);
    double correlator(int T, int n5);
    double four_point(int T, int R, int n5);
    double hamiltonian();
    void hmc(int n_sweeps = 0, bool update_dt = false, bool no_metropolis = false);
    void heat_bath(int n_sweeps = 0);
    void cool(int n5, int n_sweeps = 0);
//    void wilson_flow(int n_sweeps = 0);
    void wilson_flow(int n5, int n_sweeps = 0);
    void stout_smear(int n5, double rho = 0.1, int n_sweeps = 0);
    double field_strength(int n5);
//    double field_strength_x(int n5);
    double topological_charge(int n5);
    int thermalize(int n_min = 0, int n_max = 0);
    double ave_link_t(int n5);
    double theta(bool coulomb, int n5);
    long double relax(long double error_target, bool coulomb, int n5);
    void write_four_point(const char* filename, bool coulomb, int n5);
    void write_correlator(const char* filename, bool coulomb, int n5);
};

#endif /* su3_x_hpp */
