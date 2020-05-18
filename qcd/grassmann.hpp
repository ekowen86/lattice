//
//  grassmann.hpp
//  Lattice
//
//  Created by Evan Owen on 2/15/20.
//  Copyright © 2020 Evan Owen. All rights reserved.
//

#ifndef grassmann_hpp
#define grassmann_hpp

using namespace std;

template <class T>
class grassmann {
    
public:
    T v1;
    T v2;
    
    grassmann() {
        
    }
    
    grassmann(T v1, T v2) {
        this->v1 = v1;
        this->v2 = v2;
    }
    
    T operator* (grassmann<T> g) {
        return (v1 * g.v2 - g.v1 * v2) / 2.0;
    }

//    grassmann<T> operator+ (grassmann<T> g) {
//        return grassmann<T>(v1 + g.v1, v2 + g.v2);
//    }

    grassmann<T> operator- () {
        return grassmann<T>(-v1, -v2);
    }

    grassmann<T> operator* (complex<double> d) {
        return grassmann<T>(v1 * d, v2 * d);
    }

    grassmann<T> conj() {
        return grassmann<T>(std::conj(v2), std::conj(v1));
    }

    grassmann<T> conjugate() {
        return grassmann<T>(v2.conjugate(), v1.conjugate());
    }

    grassmann<T> transpose() {
        return grassmann<T>(v1.transpose(), v2.transpose());
    }
};

#endif /* grassmann_hpp */
