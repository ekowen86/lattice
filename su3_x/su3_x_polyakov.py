import sqlite3
import numpy as np
import matplotlib.pyplot as plt
import bz2
import scipy.optimize as opt
import cmath
import sys
import warnings

warnings.simplefilter(action='ignore', category=FutureWarning)

def complex_mean(a):
	return np.abs(np.mean(a, dtype=a[0].dtype))

def complex_stdev(a):
	return np.abs(np.std(a, dtype=a[0].dtype))

def fix_autocorr(data, n_corr):
	n = len(data)
	size = n / n_corr
	a = np.zeros(size, dtype=data[0].dtype)
	start = 0
	for n in range(0, size):
		end = np.minimum(start + n_corr, n)
		a[n] = complex_mean(data[start:end])
		start = end
	return a

def jackknife(data, f):
	n = len(data)
	a = np.zeros(n)
	da = np.zeros(n)
	a_bar = 0.0
	for i in range(0, n):
		# evaluate the function leaving out the ith data point
		a_data = np.delete(data, i)
		a[i] = f(a_data)
		da[i] = complex_stdev(a_data) / np.sqrt(n - 1)
		a_bar += a[i]
	
	fn = float(n)
	a_bar /= fn
	a_err = 0.0
	for i in range(0, n):
		a_err += (a[i] - a_bar)**2
	a_err = np.sqrt((fn - 1) / fn * a_err)
	
	return a_bar, a_err, a, da

def jackknife_deleteD(data, n_sub, n_corr):
	data = fix_autocorr(data, n_corr)
	f_sub = float(n_sub)
	a = np.zeros(n_sub, dtype=data[0].dtype)
	start = 0
	count = len(data)
	size = count / n_sub
	for n in range(0, n_sub):
		end = np.minimum(start + size, count)
		a[n] = complex_mean(data[start:end])
		start = end

	mean = complex_mean(data)
	error = 0.0
	for n in range(0, n_sub):
		sum = 0.0
		for i in range(0, n_sub):
			if (i == n):
				continue
			sum += a[i]
		a1 = (sum / (f_sub - 1)) - mean
		error += a1 * a1

	return np.sqrt((f_sub - 1) / f_sub * error, dtype=data[0].dtype)

print("\nOpening data file...")
con = sqlite3.connect("data_wp/su3_x_16_3_6000_100.db")
cursor = con.cursor()

L = int(16) # lattice size in spatial directions
L_T = int(16) # lattice size in time direction
R_half = int(L / 2) # spatial distance halfway across the lattice
p = np.empty(R_half) # polyakov measurements at each value of r
d_p = np.empty(R_half) # error in polyakov measurements
T = 2.5 # wilson flow time
R_min = 5; # smallest value of r to use
R_max = 7; # largest value of r to use

def collapse_row(row):
	arr = np.empty(len(row))
	for i in range(0, len(row)):
		arr[i] = float(row[i][0])
	return arr

print("\nReading Polyakov loop values...")
for r in range(0, R_half):
	cursor.execute("SELECT P_%d FROM polyakov_loop WHERE flow_t='%.04f'" % (r+1, T))
	row = collapse_row(cursor.fetchall())
	p[r] = np.mean(row)
	_, d_p[r], _, _ = jackknife(row, complex_mean)

cursor.close()
con.close()

for r in range(0, R_half):
	print("p(%d) = %.12f (%.12f)" % (r+1, p[r], d_p[r]))

# published beta

#   beta        sqrt(sigma)
#   5.61           0.4908
#   5.62           0.4778
#   5.63           0.4654
#   5.64           0.4535
#   5.65           0.4421
#   5.66           0.4312
#   5.67           0.4208
#   5.68           0.4108
#   5.69           0.4011
#   5.72           0.3744
#   5.78           0.3290
#   5.79           0.3223
#   5.80           0.3159
#   5.83           0.2977
#   5.85           0.2865
#   5.87           0.2759
#   5.89           0.2659
#   5.90           0.2611
#   5.92           0.2520
#   5.94           0.2433
#   5.96           0.2350
#   5.97           0.2311
#   5.98           0.2272
#   5.99           0.2234
#   6.00           0.2197
#   6.01           0.2161
#   6.02           0.2126
#   6.03           0.2092
#   6.04           0.2058
#   6.05           0.2025
#   6.06           0.1993
#   6.08           0.1931
#   6.12           0.1815
#   6.16           0.1709
#   6.20           0.1610
#   6.22           0.1564
#   6.24           0.1519
#   6.25           0.1497
#   6.26           0.1476
#   6.27           0.1455
#   6.28           0.1435
#   6.29           0.1415
#   6.30           0.1395
#   6.31           0.1375
#   6.32           0.1356
#   6.33           0.1338


# 6 * 16^3 * 1

#   beta        sqrt(sigma)                          <chi>           <DQ>
#   5.70          0.38                             0.1072(46)      1.750(56)
#   5.80          0.29                             0.0984(23)      1.190(26)
#   5.82          0.27                             0.0967(27)      1.048(29)
#   5.84          0.24                             0.0934(27)      0.864(27)
#   5.86          0.21                             0.0836(22)      0.690(26)
#   5.88          0.17                             0.0814(23)      0.495(23)
#   5.90          0.12                             0.0546(20)      0.153(19)

# 6 * 16^3 * 3 (eps5 = 0.5)

#   beta        sqrt(sigma)    beta_effective        <chi>          <DQ>
#   5.50          0.52 ?                           0.1150(50)     2.355(69)
#   5.60          0.31            ~ 5.75           0.1135(35)     2.013(41)
#   5.62          0.40                             0.1095(46)     1.822(54)
#   5.64          0.39            ~ 5.70           0.1076(47)     1.717(51)
#   5.80          0.03

# 6 * 16^3 * 5 (eps5 = 0.5)

#   beta        sqrt(sigma)    beta_effective        <chi>          <DQ>
#   5.50          0.34                             0.1155(73)     2.346(94)
#   5.60          0.37                             0.1160(78)     2.045(83)
#   5.80          0.04                             

# best fit for string tension
print("\nFitting string tension and plotting...")

R = np.zeros(R_half)
for r in range(0, R_half):
	R[r] = (r + 1.0)

def p_fit(r, sqrt_sigma, a):
	a1 = np.exp(-L_T * (a + sqrt_sigma**2 * r))
	a2 = np.exp(-L_T * (a + sqrt_sigma**2 * (L - r)))
	return a1 + a2

popt, pcov = opt.curve_fit(p_fit, R[R_min-1:R_max], p[R_min-1:R_max], [0.2, 0.05], sigma=d_p[R_min-1:R_max])
sqrt_sigma = popt[0]
d_sqrt_sigma = np.sqrt(pcov[0][0])
a = popt[1]
d_a = np.sqrt(pcov[1][1])
print("\nP(R) = exp(-L_T * (a + sigma * R)) + exp(-L_T * (a + sigma * (L - R)))")
print("a = %.12f (%.12f)" % (a, d_a))
print("sigma = %.12f (%.12f)" % (sqrt_sigma**2, 2 * np.abs(sqrt_sigma) * d_sqrt_sigma))
print("sqrt(sigma) = %.12f (%.12f)" % (np.abs(sqrt_sigma), d_sqrt_sigma))

# calculate best fit curves
R_A = np.linspace(0.0, R_half + 1, 1000)
p_A = np.zeros(len(R_A))
for r in range(0, len(R_A)):
	p_A[r] = p_fit(R_A[r], sqrt_sigma, a)

plt.rcParams.update({
    "text.usetex": False,
    "font.family": "sans-serif"})

# plot sigma vs R
plt.figure()
plt.xlim(0.0, R_half + 1.0)
plt.yscale("log")
# plt.ylim(0.0, p[0] * 10.0)
plt.errorbar(R, p, yerr=d_p, color="blue", marker='o', ms=5, mew=0.5, \
mfc='none', linestyle='none', linewidth=0.5, capsize=2.5, capthick=0.5)
plt.plot(R_A, p_A, color="blue", linewidth=0.5)
plt.xlabel("$R$")
plt.ylabel("$P(R)$")
plt.savefig("plot_P.pdf")
plt.close()
