//
//  su3_x_topocharge.cpp
//  su3_x_topocharge
//
//  Created by Evan Owen on 5/26/20.
//  Copyright © 2020 Evan Owen. All rights reserved.
//
//  usage: su3_topological_charge_5d output N T N5 D beta eps5 n_therm n_sweeps n_data t_wf n_cool
//         output: output data filename
//         N: lattice size (spacial dimensions)
//         T: lattice size (time dimension)
//         N5: number of coupled parallel lattices in extra dimension
//         D: number of dimensions (excluding extra dimension)
//         beta: gauge coupling
//         eps5: ratio of coupling in extra dimension to normal coupling
//         n_therm: minimum number of thermalization sweeps
//         n_sweeps: number of sweeps per data collection
//         n_data: number of data values to collect
//         t_wf: wilson flow time
//         n_cool: number of cooling sweeps

#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <thread>
#include <sqlite3.h>

#include "su3_x.hpp"

using namespace std;

const char* filename;
int N, T, N5, D, n_therm, n_sweeps, n_data, n_cool;
double beta, eps5, t_wf;
void init_db();
void write_parameters(su3_x_lattice* lattice);
void write_action(su3_x_lattice*, double t);
void write_wilson_loop(su3_x_lattice*, double t);
void write_topocharge(su3_x_lattice*, int n_c);
void write_polyakov_loop(su3_x_lattice*, double t);

int main(int argc, const char* argv[]) {

    filename = argv[1]; cout << "output = " << filename << endl;
    N = atoi(argv[2]); cout << "N = " << N << ", ";
    T = atoi(argv[3]); cout << "T = " << T << ", ";
    N5 = atoi(argv[4]); cout << "N5 = " << N5 << ", ";
    D = atoi(argv[5]); cout << "D = " << D << endl;
    cout << setprecision(4) << fixed;
    beta = stod(argv[6]); cout << "beta = " << beta << ", ";
    eps5 = stod(argv[7]); cout << "eps5 = " << eps5 << endl;
    n_therm = atoi(argv[8]); cout << "n_therm = " << n_therm << ", ";
    n_sweeps = atoi(argv[9]); cout << "n_sweeps = " << n_sweeps << ", ";
    n_data = atoi(argv[10]); cout << "n_data = " << n_data << ", ";
    t_wf = stod(argv[11]); cout << "t_wf = " << t_wf << ", ";
    n_cool = atoi(argv[12]); cout << "n_cool = " << n_cool << endl;
    cout << thread::hardware_concurrency() << " parallel cores available" << endl;

    cout << setprecision(6) << fixed;

    // create a lattice and thermalize it
    su3_x_lattice lattice = su3_x_lattice(N, T, N5, D, beta, eps5);
    lattice.parallel = true;
    cout << "thermalizing...";
    n_therm = lattice.thermalize(n_therm);
    cout << " done" << endl;
    cout << "hmc time per config. = " << setprecision(4) << fixed << lattice.dt * lattice.n_steps * n_data << endl;

    init_db();
    write_parameters(&lattice);

    int n5_center = lattice.n5_center;
    for (int n = 0; n < n_data; n++) {

//        if (n != 0) lattice.heat_bath(n_sweeps);
        if (n != 0) lattice.hmc(n_sweeps);

        cout << "writing configuration " << n << "...";

        su3_x_lattice wf_lattice = su3_x_lattice(&lattice);
//        write_action(&wf_lattice, 0.0);
        write_wilson_loop(&wf_lattice, 0.0);
        write_polyakov_loop(&wf_lattice, 0.0);
        for (double t = 0.02; t < (t_wf + 0.01); t += 0.02) {
            wf_lattice.wilson_flow(n5_center, 0.02);
//            write_action(&wf_lattice, t);
            write_wilson_loop(&wf_lattice, t);
            write_polyakov_loop(&wf_lattice, t);
        }

        su3_x_lattice cool_lattice = su3_x_lattice(&lattice);
        write_topocharge(&cool_lattice, 0);
        for (int n = 100; n <= n_cool; n += 100) {
            cool_lattice.cool(n5_center, 100);
            write_topocharge(&cool_lattice, n);
        }
        cout << " done" << endl;
    }

    return 0;
}

void init_db() {
    // open sqlite database
    sqlite3* db;
    char* zErrMsg;
    int ret_code;
    stringstream ss;
    int N_half = N >> 1;

    // open the sqlite database file
    ret_code = sqlite3_open(filename, &db);

    if (ret_code) {
        cerr << "can't open database: " << sqlite3_errmsg(db) << endl;
        return;
    }

    // create the parameter table if it doesn't exist
    ss << "CREATE TABLE IF NOT EXISTS parameters (N INTEGER, T INTEGER, N_5 INTEGER, D INTEGER, beta VARCHAR)";

    ret_code = sqlite3_exec(db, ss.str().c_str(), NULL, 0, &zErrMsg);

    if (ret_code != SQLITE_OK) {
        cerr << "sql error: " << zErrMsg << endl;
        sqlite3_free(zErrMsg);
    }

//    // create the action table if it doesn't exist
//    ss.str("");
//    ss << "CREATE TABLE IF NOT EXISTS action (flow_t INTEGER, S VARCHAR)";
//
//    ret_code = sqlite3_exec(db, ss.str().c_str(), NULL, 0, &zErrMsg);
//
//    if (ret_code != SQLITE_OK) {
//        cerr << "sql error: " << zErrMsg << endl;
//        sqlite3_free(zErrMsg);
//    }

    // create the wilson loop table if it doesn't exist
    ss.str("");
    ss << "CREATE TABLE IF NOT EXISTS wilson_loop (flow_t VARCHAR, plaq VARCHAR";
    for (int n = 1; n < N_half; n++) {
        ss << ", W_" << n << n << " VARCHAR";
        ss << ", W_" << n << (n + 1) << " VARCHAR";
    }
    ss << ", W_" << N_half << N_half << " VARCHAR)";

    ret_code = sqlite3_exec(db, ss.str().c_str(), NULL, 0, &zErrMsg);

    if(ret_code != SQLITE_OK) {
        fprintf(stderr, "SQL error: %s\n", zErrMsg);
        sqlite3_free(zErrMsg);
    }

    // create the topological charge table if it doesn't exist
    ss.str("");
    ss << "CREATE TABLE IF NOT EXISTS topocharge (n_c INTEGER, Q VARCHAR)";

    ret_code = sqlite3_exec(db, ss.str().c_str(), NULL, 0, &zErrMsg);

    if (ret_code != SQLITE_OK) {
        cerr << "sql error: " << zErrMsg << endl;
        sqlite3_free(zErrMsg);
    }

    // create the polyakov loop table if it doesn't exist
    ss.str("");
    ss << "CREATE TABLE IF NOT EXISTS polyakov_loop (flow_t VARCHAR";
    for (int n = 1; n <= N_half; n++) {
        ss << ", P_" << n << " VARCHAR";
    }
    ss << ")";

    ret_code = sqlite3_exec(db, ss.str().c_str(), NULL, 0, &zErrMsg);

    if(ret_code != SQLITE_OK) {
        cerr << "sql error: " << zErrMsg << endl;
        sqlite3_free(zErrMsg);
    }

    // close sqlite database
    sqlite3_close(db);
}

void write_parameters(su3_x_lattice* lattice) {

    // open sqlite database
    sqlite3* db;
    char *zErrMsg = 0;
    int ret_code;

    ret_code = sqlite3_open(filename, &db);

    if (ret_code) {
        cerr << "can't open database: " << sqlite3_errmsg(db) << endl;
        return;
    }

    // write parameters to database
    stringstream ss;
    ss << "INSERT INTO parameters (N, T, N_5, D, beta) VALUES (";
    ss << N << ", ";
    ss << T << ", ";
    ss << N5 << ", ";
    ss << D << ", ";
    ss << setprecision(4) << fixed;
    ss << "'" << beta << "')";
    ret_code = sqlite3_exec(db, ss.str().c_str(), NULL, 0, &zErrMsg);

    if(ret_code != SQLITE_OK) {
        cerr << "sql error: " << zErrMsg << endl;
        sqlite3_free(zErrMsg);
    }

    // close sqlite database
    sqlite3_close(db);
}

void write_action(su3_x_lattice* lattice, double t) {

    // calculate action
    int n5_center = lattice->n5_center;

    double action = lattice->action(n5_center);

    // open sqlite database
    sqlite3* db;
    char *zErrMsg = 0;
    int ret_code;

    ret_code = sqlite3_open(filename, &db);

    if (ret_code) {
        cerr << "can't open sqlite database: " << sqlite3_errmsg(db) << endl;
        return;
    }

    // write observables to database
    stringstream ss;
    ss << "INSERT INTO action (flow_t, S) VALUES (";
    ss << setprecision(4) << fixed;
    ss << "'" << t << "'";
    ss << setprecision(8) << scientific;
    ss << ", '" << action << "')";
    ret_code = sqlite3_exec(db, ss.str().c_str(), NULL, 0, &zErrMsg);

    if(ret_code != SQLITE_OK) {
        cerr << "sql error: " << zErrMsg << endl;
        sqlite3_free(zErrMsg);
    }

    // close sqlite database
    sqlite3_close(db);
}

void write_wilson_loop(su3_x_lattice* lattice, double t) {

    // calculate observables
    int n5_center = lattice->n5_center;

    // average plaquette (1x1 wilson loop)
    double plaq = lattice->plaq(n5_center);

    int N_half = N >> 1;
    double w[N_half * 2 - 1];
    for (int R = 1; R < N_half; R++) {
        w[R * 2 - 2] = lattice->wilson_loop(R, R, n5_center);
        w[R * 2 - 1] = (lattice->wilson_loop(R, R + 1, n5_center) + lattice->wilson_loop(R + 1, R, n5_center)) * 0.5;
    }
    w[N_half * 2 - 2] = lattice->wilson_loop(N_half, N_half, n5_center);

    // open sqlite database
    sqlite3* db;
    char *zErrMsg = 0;
    int ret_code;

    ret_code = sqlite3_open(filename, &db);

    if (ret_code) {
        cerr << "can't open database: " << sqlite3_errmsg(db) << endl;
        return;
    }

    // write observables to database
    stringstream ss;
    ss << "INSERT INTO wilson_loop (flow_t, plaq";
    for (int n = 1; n < N_half; n++) {
        ss << ", W_" << n << n;
        ss << ", W_" << n << (n + 1);
    }
    ss << ", W_" << N_half << N_half << ") VALUES (";
    ss << setprecision(4) << fixed;
    ss << "'" << t << "', ";
    ss << setprecision(8) << scientific;
    ss << "'" << plaq << "', ";
    for (int n = 0; n < (2 * N_half - 2); n++) {
        ss << "'" << w[n] << "', ";
    }
    ss << "'" << w[2 * N_half - 2] << "')";
    ret_code = sqlite3_exec(db, ss.str().c_str(), NULL, 0, &zErrMsg);

    if(ret_code != SQLITE_OK) {
        fprintf(stderr, "sql error: %s\n", zErrMsg);
        sqlite3_free(zErrMsg);
    }

    // close sqlite database
    sqlite3_close(db);

}

void write_topocharge(su3_x_lattice* lattice, int n_c) {

    // calculate topological charge
    double topological_charge = lattice->topological_charge(lattice->n5_center);

    // open sqlite database
    sqlite3* db;
    char *zErrMsg = 0;
    int ret_code;

    ret_code = sqlite3_open(filename, &db);

    if (ret_code) {
        cerr << "can't open database: " << sqlite3_errmsg(db) << endl;
        return;
    }

    // write observables to database
    stringstream ss;
    ss << "INSERT INTO topocharge (n_c, Q) VALUES (";
    ss << n_c << ", ";
    ss << setprecision(8) << scientific;
    ss << "'" << topological_charge << "')";
    ret_code = sqlite3_exec(db, ss.str().c_str(), NULL, 0, &zErrMsg);

    if(ret_code != SQLITE_OK) {
        cerr << "sql error: " << zErrMsg << endl;
        sqlite3_free(zErrMsg);
    }

    // close sqlite database
    sqlite3_close(db);

}

void write_polyakov_loop(su3_x_lattice* lattice, double t) {

    // calculate observables
    int n5_center = lattice->n5_center;

    int N_half = N >> 1;
    double p[N_half];
    for (int R = 0; R < N_half; R++) {
        p[R] = lattice->polyakov_loop(R + 1, n5_center);
    }

    // open sqlite database
    sqlite3* db;
    char *zErrMsg = 0;
    int ret_code;

    ret_code = sqlite3_open(filename, &db);

    if (ret_code) {
        cerr << "can't open sqlite database: " << sqlite3_errmsg(db) << endl;
        return;
    }

    // write observables to database
    stringstream ss;
    ss << "INSERT INTO polyakov_loop (flow_t";
    for (int n = 1; n <= N_half; n++) {
        ss << ", P_" << n;
    }
    ss << ") VALUES (";
    ss << setprecision(4) << fixed;
    ss << "'" << t << "'";
    ss << setprecision(8) << scientific;
    for (int n = 0; n < N_half; n++) {
        ss << ", '" << p[n] << "'";
    }
    ss << ")";
    ret_code = sqlite3_exec(db, ss.str().c_str(), NULL, 0, &zErrMsg);

    if(ret_code != SQLITE_OK) {
        cerr << "sql error: " << zErrMsg << endl;
        sqlite3_free(zErrMsg);
    }

    // close sqlite database
    sqlite3_close(db);
}
