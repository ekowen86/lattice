//
// fft.ccp
// fast fourier transform
//

#include <iostream>
#include <complex>
#include <vector>
#include <chrono>
#include <math.h>

using namespace std;

typedef vector<double> Array;
typedef vector<complex<double>> CArray;

CArray slowFFT(CArray& x) {
    auto N = x.size();
    CArray X(N);
    fill(X.begin(), X.end(), 0.0);

    for (auto k = 0; k < N; k++) {
        for (auto n = 0; n < N; n++) {
            X[k] += x[n] * polar(1.0, -2.0 * M_PI * double(k * n) / double(N));
        }
    }
    return X;
}

CArray recursiveFFT(CArray& x) {
    auto N = x.size();
    if (N <= 1) return CArray(x);

    // split the array into even and odd parts
    CArray xEven;
    CArray xOdd;

    for (auto n = 0; n < N; n++) {
        if ((n % 2) == 0) {
            xEven.push_back(x[n]);
        } else {
            xOdd.push_back(x[n]);
        }
    }

    // recursively compute FT of both parts
    CArray XEven = recursiveFFT(xEven);
    CArray XOdd = recursiveFFT(xOdd);

    // combine results
    CArray X(N);
    for (auto k = 0; k < N / 2; k++) {
        X[k] = XEven[k] + XOdd[k] * polar(1.0, -2.0 * M_PI * k / double(N));
        X[k + N / 2] = XEven[k] - XOdd[k] * polar(1.0, -2.0 * M_PI * k / double(N));
    }
    return X;
}

void bitReversalFFT(CArray &x)
{
	// DFT
	unsigned int N = x.size();
    unsigned int k = N;
    unsigned int n;
	double thetaT = M_PI / N;
	complex<double> phiT(cos(thetaT), -sin(thetaT));
    complex<double> T;
	while (k > 1) {
		n = k;
		k >>= 1;
		phiT = phiT * phiT;
		T = 1.0L;
		for (unsigned int l = 0; l < k; l++)
		{
			for (unsigned int a = l; a < N; a += n)
			{
				unsigned int b = a + k;
				complex<double> t = x[a] - x[b];
				x[a] += x[b];
				x[b] = t * T;
			}
			T *= phiT;
		}
	}
	// Decimate
	unsigned int m = (unsigned int)log2(N);
	for (unsigned int a = 0; a < N; a++)
	{
		unsigned int b = a;
		// Reverse bits
		b = (((b & 0xaaaaaaaa) >> 1) | ((b & 0x55555555) << 1));
		b = (((b & 0xcccccccc) >> 2) | ((b & 0x33333333) << 2));
		b = (((b & 0xf0f0f0f0) >> 4) | ((b & 0x0f0f0f0f) << 4));
		b = (((b & 0xff00ff00) >> 8) | ((b & 0x00ff00ff) << 8));
		b = ((b >> 16) | (b << 16)) >> (32 - m);
		if (b > a)
		{
			complex<double> t = x[a];
			x[a] = x[b];
			x[b] = t;
		}
	}
}

double constantFn(double x) {
    return 1.0;
}

double gaussianFn(double x) {
    return exp(-x * x);
}

double sinFn(double x) {
    return sin(M_PI * x);
}

double cosFn(double x) {
    return cos(M_PI * x);
}

CArray discretizeFn(double f(double), double a, double b, int N) {
    CArray x(N);
    double step = (b - a) / double(N);
    for (auto n = 0; n < N; n++) {
        x[n] = complex<double>(f(a + n * step), 0.0);
    }
    return x;
}

void printResult(CArray X) {
    cout << fixed;
    cout.precision(12);

    for (int k = 0; k < X.size(); k++) {
        cout << k << ": " << X[k] << endl;
    }
}

int main(int argc, const char* argv[]) {

    CArray x = discretizeFn(gaussianFn, -12.0, 12.0, pow(2,14));

    {
        auto start = chrono::high_resolution_clock::now();
        slowFFT(x);
        auto end = chrono::high_resolution_clock::now();
        cout << "slowFFT: " << chrono::duration_cast<chrono::milliseconds>(end - start).count() << "ms" << endl;
    }

    {
        auto start = chrono::high_resolution_clock::now();
        recursiveFFT(x);
        auto end = chrono::high_resolution_clock::now();
        cout << "recursiveFFT: " << chrono::duration_cast<chrono::milliseconds>(end - start).count() << "ms" << endl;
    }

    CArray xBitReverse = CArray(x);
    {
        auto start = chrono::high_resolution_clock::now();
        bitReversalFFT(xBitReverse);
        auto end = chrono::high_resolution_clock::now();
        cout << "bitReversalFFT: " << chrono::duration_cast<chrono::milliseconds>(end - start).count() << "ms" << endl;
    }

    return 0;
}
