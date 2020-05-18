//
//  random_box_distance.cpp
//  random_box_distance
//
//  Created by Evan Owen on 2/10/19.
//  Copyright © 2019 Evan Owen. All rights reserved.
//

#include <iostream>
#include <cmath>

int main(int argc, const char * argv[]) {

    int a = 12;
    int b = 8;
    int c = 4;
    long double sum = 0.0L;
    double d = 0.0;
    double error = 1.0;
    int count = 0;
    
    // seed the rng
    srand((unsigned int)time(0));
    
    while (error > 1e-7) {
        
        // pick two points at random
        int x1 = rand() % a;
        int y1 = rand() % b;
        int z1 = rand() % c;
        int x2 = rand() % a;
        int y2 = rand() % b;
        int z2 = rand() % c;

        // make sure they're not the same point
        if ((x1 == x2) && (y1 == y2) && (z1 == z2)) continue;
        
        // find the distance between them
        double r2 = 0.0;
        r2 += double(x2 - x1) * double(x2 - x1);
        r2 += double(y2 - y1) * double(y2 - y1);
        r2 += double(z2 - z1) * double(z2 - z1);
        double r = sqrt(r2);
        
        if (count < 100) {
            sum += r;
            count++;
            continue;
        }
        
        // add it to the sum
        count++;
        sum += r;
        
        if ((count % 10000000) != 0) continue;

        // recalculate the average and the error
        double d1 = sum / count;
        error = abs(d - d1);
        d = d1;
        printf("%.10f\n", d);
    }
    printf("%.10Lf\n", sum / count);

    return 0;
}
