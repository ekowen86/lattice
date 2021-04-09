//
//  su3_x_grid.cpp
//  su3_x_grid
//
//  Created by Evan Owen on 12/9/20.
//

#include <iostream>
#include <fstream>
#include <iomanip>

#include "su3_x.hpp"

using namespace std;

void write_wilson_loop(const char* path, su3_x_lattice* lattice, double t);
void write_polyakov_loop(const char* path, su3_x_lattice* lattice, double t);
double write_field_strength(const char* path, su3_x_lattice* lattice, double t);
void write_topocharge(const char* path, su3_x_lattice* lattice, double n_c);

int main(int argc, const char * argv[]) {
    
    const char* path = argv[1]; printf("File path: %s\n", path);
    int id = atoi(argv[2]); printf("Configuration ID: %d\n", id);
    int N = atoi(argv[3]); printf("N: %d\n", N);
    int T = atoi(argv[4]); printf("T: %d\n", T);
    int N5 = atoi(argv[5]); printf("N5: %d\n", N5);
    int D = atoi(argv[6]); printf("D: %d\n", D);
    double beta = stod(argv[7]); printf("beta: %lf\n", beta);
    double eps5 = stod(argv[8]); printf("epsilon5: %lf\n", eps5);
    double t_wf = stod(argv[9]); printf("Wilson Flow time: %lf\n", t_wf);
    int n_cool = atoi(argv[10]); printf("Number of Cooling Sweeps: %d\n", n_cool);

    char ckptPath[50];
    sprintf(ckptPath, "%s/ckpoint_lat.%d", path, id);
    printf("Opening file: %s\n", ckptPath);
    
    // create a lattice from a grid data file
    ifstream ckptFile(ckptPath, ifstream::in);
    if (!ckptFile.is_open()) {
        printf("Error opening lattice data file\n");
        return 1;
    }
    su3_x_lattice lattice(N, T, N5, D, beta, eps5, ckptFile, false);
    ckptFile.close();
    
    lattice.parallel = true;
    
//    // write the average plaquette and id
//    char plaqPath[50];
//    sprintf(plaqPath, "%s/plaq.dat", path);
//    ofstream of(plaqPath, ofstream::out | ofstream::app);
//    of << id;
//    of << setprecision(8) << scientific;
//    for (int n5 = 0; n5 < lattice.N5; n5++) {
//        of << " " << lattice.plaq(n5) / 3.0;
//    }
//    of << endl;
//    of.close();
    
//    printf("Writing Wilson and Polyakov Loops with Wilson Flow...\n");
//
//    char wloopPath[50];
//    sprintf(wloopPath, "%s/wilson_loop.%d", path, id);
//    char ploopPath[50];
//    sprintf(ploopPath, "%s/polyakov_loop.%d", path, id);
//
//    // wilson flow and measure wilson loops
//    su3_x_lattice wf_lattice = su3_x_lattice(&lattice);
//    write_wilson_loop(wloopPath, &wf_lattice, 0.0);
//    write_polyakov_loop(ploopPath, &wf_lattice, 0.0);
//    for (double t = 0.02; t < (t_wf + 0.01); t += 0.02) {
//        wf_lattice.wilson_flow(lattice.n5_center);
//        write_wilson_loop(wloopPath, &wf_lattice, t);
//        write_polyakov_loop(ploopPath, &wf_lattice, t);
//    }

    printf("Measuring t0 with Wilson flow...\n");
    char t0Path[50];
    sprintf(t0Path, "%s/t0.dat", path);
    char fieldStrengthPath[50];
    sprintf(fieldStrengthPath, "%s/field_strength.%d", path, id);

    // wilson flow and measure t0
    su3_x_lattice t0_lattice = su3_x_lattice(&lattice);
    double E_t2 = write_field_strength(fieldStrengthPath, &t0_lattice, 0.0);
    double dt = t0_lattice.wf_dt;
    bool reached_t0 = false;
    for (double t = dt; E_t2 <= 0.5; t += dt) {
        t0_lattice.wilson_flow();
//        t0_lattice.wilson_flow(lattice.n5_center);
        E_t2 = write_field_strength(fieldStrengthPath, &t0_lattice, t);
        if (!reached_t0 && E_t2 > 0.3) {
            reached_t0 = true;
            ofstream of(t0Path, ofstream::out | ofstream::app);
            of << id;
            of << setprecision(8) << scientific;
            of << " " << t << endl;
            of.close();
        }
    }

//    printf("Writing Topological Charge with Link Cooling...\n");
//    char topochargePath[50];
//    sprintf(topochargePath, "%s/topo_charge.%d", path, id);
//
//    // cool and measure topological charge
//    su3_x_lattice cool_lattice = su3_x_lattice(&lattice);
//    write_topocharge(topochargePath, &cool_lattice, 0);
//    for (int n = 100; n <= n_cool; n += 100) {
//        cool_lattice.cool(lattice.n5_center, 100);
//        write_topocharge(topochargePath, &cool_lattice, n);
//    }

    printf("Done!\n\n");
}

void write_wilson_loop(const char* path, su3_x_lattice* lattice, double t) {
    // calculate wilson loops
    int n5_center = lattice->n5_center;

    // average plaquette (1x1 wilson loop)
    double plaq = lattice->plaq(n5_center);

    int N_half = lattice->N >> 1;
    double w[N_half * 2 - 1];
    for (int R = 1; R < N_half; R++) {
        w[R * 2 - 2] = lattice->wilson_loop(R, R, n5_center);
        w[R * 2 - 1] = (lattice->wilson_loop(R, R + 1, n5_center) + lattice->wilson_loop(R + 1, R, n5_center)) * 0.5;
    }
    w[N_half * 2 - 2] = lattice->wilson_loop(N_half, N_half, n5_center);

    // write to output file
    ofstream of(path, ofstream::out | ofstream::app);
    of << setprecision(4) << fixed;
    of << t;
    of << setprecision(8) << scientific;
    of << " " << plaq;
    for (int n = 0; n < (2 * N_half - 2); n++) {
        of << " " << w[n];
    }
    of << " " << w[2 * N_half - 2] << endl;
    of.close();
}

void write_polyakov_loop(const char* path, su3_x_lattice* lattice, double t) {
    // calculate observables
    int n5_center = lattice->n5_center;

    int N_half = lattice->N >> 1;
    double p[N_half];
    for (int R = 0; R < N_half; R++) {
        p[R] = lattice->polyakov_loop(R + 1, n5_center);
    }

    // write to output file
    ofstream of(path, ofstream::out | ofstream::app);
    of << setprecision(4) << fixed;
    of << t;
    of << setprecision(8) << scientific;
    for (int n = 0; n < N_half; n++) {
        of << " " << p[n];
    }
    of << endl;
    of.close();
}

double write_field_strength(const char* path, su3_x_lattice* lattice, double t) {

    // calculate field strength
    double E = lattice->field_strength(lattice->n5_center);
    double E_t2 = E * t * t;
    
    // write to output file
    ofstream of(path, ofstream::out | ofstream::app);
    of << setprecision(4) << fixed;
    of << t;
    of << setprecision(8) << scientific;
    of << " " << E << " " << E_t2 << endl;
    of.close();
    return E_t2;
}

void write_topocharge(const char* path, su3_x_lattice* lattice, double n_c) {

    // calculate topological charge
    double topological_charge = lattice->topological_charge(lattice->n5_center);

    // write to output file
    ofstream of(path, ofstream::out | ofstream::app);
    of << n_c;
    of << setprecision(8) << scientific;
    of << " " << topological_charge << endl;
    of.close();
}
