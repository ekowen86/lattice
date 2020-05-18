//
//  jackknife.cpp
//
//  Created by Evan Owen on 9/1/18.
//  Copyright © 2018 Evan Owen. All rights reserved.
//

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cmath>
#include <cstdlib>

#define MAX_VALUES 10000
#define MAX_DATA 500
//#define CL_FOURPOINT

using namespace std;

class jk_data {
public:
    double val[MAX_VALUES];
    int count = 0;
    
    int add_value(double value) {
        val[count++] = value;
        return count;
    }
    
    double calc_mean(int start = 0, int end = 0) {
        if (end == 0) end = count;
        double sum = 0.0;
        for (int i = start; i < end; i++) {
            sum += val[i];
        }
        return sum / double(end - start);
    }
    
    double calc_error(int n_sub) {
        // split the data set up into n_sub equally sized subsets and compute the average of each
        double a[n_sub];
        int start = 0;
        int size = count / n_sub;
        for (int n = 0; n < n_sub; n++) {
            int end = min(start + size, count);
            a[n] = calc_mean(start, end);
            start = end;
        }
        
        // calculate the average excluding each data set
        double mean = calc_mean();
        double error = 0.0;
        for (int n = 0; n < n_sub; n++) {
            double sum = 0.0;
            for (int i = 0; i < n_sub; i++) {
                if (i == n) continue;
                sum += a[i];
            }
            double a1 = (sum / double(n_sub - 1)) - mean;
            error += a1 * a1;
        }
        
        // calculate the jackknife error
        return sqrt((double(n_sub) - 1) / double(n_sub) * error);
    }
};

void parse_line(string line, jk_data* data);

#ifdef CL_FOURPOINT

// this is the main function for coulomb/landau four-point functions
//
//  usage: jackknife N beta gauge n_sub
//         N: lattice size
//         Nt: lattice size (time direction)
//         n_sub: number of data subsets

// arrays for the data
jk_data four_point[MAX_DATA];
jk_data correlators[MAX_DATA];

int main(int argc, const char * argv[]) {
    
    int r_max = atoi(argv[1]);
    int t_max = (argc == 4) ? atoi(argv[2]) : r_max;

    int n_sub = atoi(argv[3]);
    int n_sub_corr = 10;
    
    bool coulomb = strcmp(argv[3], "C") ? false : true;
    char filename[100];
    sprintf(filename, "data_%s_%s_%s.txt", argv[1], argv[2], argv[3]);
    
    // open input file
    ifstream input;
    input.open(filename, ifstream::in);
    
    // parse each line
    string line;
    while (!input.eof()) {
        getline(input, line);
        if (line.compare("") == 0) continue; // skip blank lines
        parse_line(line, four_point);
    }
    
    if (coulomb) {
        // coulomb gauge
        cout << "#T\tR\t<Gc>\td<Gc>\n";
        for (int i = 0; i < MAX_DATA; i++) {
            int R = i % r_max;
            int T = i / r_max;
            if (R >= r_max) break;
            if (T >= t_max) break;
            
            int x = R * t_max + T;
            if (four_point[x].count == 0) break;
            double mean = four_point[x].calc_mean();
            double error = four_point[x].calc_error(n_sub);
            
            cout << setprecision(12);
            cout << T + 1 << "\t";
            cout << R + 1 << "\t";
            cout << mean << "\t" << error << "\n";
        }
    } else {
        // landau gauge
        
        // open correlator file
        char filename_corr[100];
        sprintf(filename_corr, "corr_%s_%s_%s.txt", argv[1], argv[2], argv[3]);
        ifstream input_corr;
        input_corr.open(filename_corr, ifstream::in);
        
        // read correlator data
        while (!input_corr.eof()) {
            getline(input_corr, line);
            if (line.compare("") == 0) continue; // skip blank lines
            parse_line(line, correlators);
        }
        
        // close the correlator file
        input_corr.close();
        
        cout << "#T\tR\t<Gc>\td<Gc>\t<G>\td<G>\t<L>^2\td(<L>^2)\n";
        for (int i = 0; i < MAX_DATA; i++) {
            int R = i % r_max;
            int T = i / r_max;
            int x = R * t_max + T;
            
            if (correlators[T].count == 0) break;
            
            double mean = four_point[x].calc_mean();
            double mean_error = four_point[x].calc_error(n_sub);
            
            double corr = correlators[T].calc_mean();
            double corr_error = correlators[T].calc_error(n_sub_corr);
            
            // square the correlator value and its error
            corr_error = 2 * corr * corr_error;
            corr = corr * corr;
            
            // subtract out the correlator value
            double error = sqrt(mean_error * mean_error + corr_error * corr_error);
            cout << setprecision(12);
            cout << T + 1 << "\t";
            cout << R + 1 << "\t";
            cout << (mean - corr) << "\t" << error << "\t";
            cout << mean << "\t" << mean_error << "\t";
            cout << corr << "\t" << corr_error << "\n";
        }
    }
    
    // close the file and exit
    input.close();
    return 0;
}

void parse_line(string line, jk_data* data) {
    int v = 0;
    char line_1[line.length()];
    line.copy(line_1, line.length());
    char* token;
    token = strtok(line_1, "\t"); // skip lattice size
    token = strtok(NULL, "\t"); // skip beta value
    token = strtok(NULL, "\t"); // skip gauge
    token = strtok(NULL, "\t"); // this is the first data value
    while (token) {
        data[v++].add_value(strtod(token, NULL));
        token = strtok(NULL, "\t"); // get the next data value
    }
}

#else

// this is the main function for a generic dataset
//
//  usage: jackknife input n_sub
//         input: input filename
//         n_sub: number of data subsets

jk_data data[MAX_DATA];

int main(int argc, const char * argv[]) {
    const char* filename = argv[1];
    int n_sub = atoi(argv[2]);

    // open input file
    ifstream input;
    input.open(filename, ifstream::in);
    if (!input.good()) return 1;

    // parse each line
    string line;
    while (!input.eof()) {
        getline(input, line);
        if (line.compare("") == 0) continue; // skip blank lines
        parse_line(line, data);
    }

    cout << setprecision(12);
    for (int i = 0; i < MAX_DATA; i++) {
        if (data[i].count == 0) break;
        cout << i << "\t";
        cout << data[i].calc_mean() << "\t";
        cout << data[i].calc_error(n_sub) << "\n";
    }

    // close the file and exit
    input.close();
    return 0;
}

void parse_line(string line, jk_data* data) {
    int v = 0;
    char line_1[line.length()];
    line.copy(line_1, line.length());
    char* token;
    token = strtok(line_1, "\t"); // skip lattice size
    while (token) {
        data[v++].add_value(strtod(token, NULL));
        token = strtok(NULL, "\t"); // get the next data value
    }
}

#endif

