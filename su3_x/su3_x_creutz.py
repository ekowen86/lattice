import sqlite3
import numpy as np
import matplotlib.pyplot as plt
import bz2
import scipy.optimize as opt
import cmath
import sys
import warnings

warnings.simplefilter(action='ignore', category=FutureWarning)

id = "4_1_5900"
L = int(4) # lattice size
R_half = int(L / 2) # half lattice size
R_min = 1 # min R value for fit
R_max = 1 # max R value for fit

def creutz_ratio(w00, w01, w11):
	w00_bar = np.mean(w00)
	w01_bar = np.mean(w01)
	w11_bar = np.mean(w11)
	return -np.log(w00_bar * w11_bar / w01_bar**2.0)


def jackknife_creutz(w00, w01, w11, n_sub):
	f_sub = float(n_sub) # number of data subsets
	start = 0 # start of current block of data
	count = len(w00) # number of data points
	size = count / n_sub # size of each block
	
	chi_bar = creutz_ratio(w00, w01, w11)
	d_chi = 0.0
	
	for n in range(0, n_sub):
		end = np.minimum(start + size, count)
		r = range(start, end)
		
		# copy each array and delete the current subset
		w00_d = np.delete(w00, r)
		w01_d = np.delete(w01, r)
		w11_d = np.delete(w11, r)
		
		# calculate the creutz ratio and add to the error
		d_chi += (creutz_ratio(w00_d, w01_d, w11_d) - chi_bar)**2.0

		start = end

	d_chi = np.sqrt((f_sub - 1) / f_sub * d_chi)
	return (chi_bar, d_chi)


def collapse_row(row):
	arr = np.empty(len(row))
	for i in range(0, len(row)):
		arr[i] = float(row[i][0])
	return arr


print("\nOpening data file...")
con = sqlite3.connect("data/su3_x_" + id + ".db")
cursor = con.cursor()

chi = np.empty(R_half - 1)
d_chi = np.empty(R_half - 1)
R = np.zeros(R_half - 1) # R values
_R2 = np.zeros(R_half - 1) # 1/R^2 values

print("\nCalculating Creutz ratios...")
for r in range(0, R_half - 1):

	# use a wilson flow time with a smearing radius of R - 1 when
	# calculating a Creutz ratio chi(R,R), rounded down to the nearest
	# multiple of 0.02. Note that r = R - 1
	# R_smear = sqrt(8 * T)
 	T = np.floor(r**2.0 / 8.0 * 50.0) / 50.0

	# get R and 1/R^2
	R[r] = (r + 1.0)
	_R2[r] = 1 / R[r]**2.0

	# get W(R,R)
	cursor.execute("SELECT W_%d%d FROM wilson_loop WHERE flow_t='%.04f'" % (r+1, r+1, T))
	w00 = collapse_row(cursor.fetchall())
	###
# 	w00 = w00[100:200]
	###
	
	# get W(R+1,R)
	cursor.execute("SELECT W_%d%d FROM wilson_loop WHERE flow_t='%.04f'" % (r+1, r+2, T))
	w01 = collapse_row(cursor.fetchall())
	###
# 	w01 = w01[100:200]
	###
	
	# get W(R+1,R+1)
	cursor.execute("SELECT W_%d%d FROM wilson_loop WHERE flow_t='%.04f'" % (r+2, r+2, T))
	w11 = collapse_row(cursor.fetchall())
	###
# 	w11 = w11[100:200]
	###

	# calculate Creutz ratio
	chi[r], d_chi[r] = jackknife_creutz(w00, w01, w11, 20)
	print("chi(%d,%d) = %.12f (%.12f)" % (r+1, r+1, chi[r], d_chi[r]))
	
cursor.close()
con.close()

# best fit for string tension
print("\nFitting string tension and plotting...")

def sigma_fit(R, sqrt_sigma, B):
	return sqrt_sigma**2.0 + B / R**2.0
# 	return sqrt_sigma**2 + B * (1 / R**2 + 1 / (L - R)**2)
# 	return sqrt_sigma**2 + B * (L / (R * (L - R)))**2

popt, pcov = opt.curve_fit(sigma_fit, R[R_min-1:R_max], chi[R_min-1:R_max], [0.05, 1.0], sigma=d_chi[R_min-1:R_max])
sqrt_sigma = popt[0]
d_sqrt_sigma = np.sqrt(pcov[0][0])
B = popt[1]
d_B = np.sqrt(pcov[1][1])
print("\nchi(R,R) = sigma + B / R^2")
print("B = %.12f (%.12f)" % (B, d_B))
print("sigma = %.12f (%.12f)" % (sqrt_sigma**2, 2 * np.abs(sqrt_sigma) * d_sqrt_sigma))
print("sqrt(sigma) = %.12f (%.12f)" % (np.abs(sqrt_sigma), d_sqrt_sigma))

# calculate best fit curves
R_A = np.linspace(0.001, R_half, 1000)
_R2_A = np.zeros(len(R_A))
chi_A = np.zeros(len(R_A))
for r in range(0, len(R_A)):
	chi_A[r] = sigma_fit(R_A[r], sqrt_sigma, B)
	_R2_A[r] = 1 / R_A[r]**2
# 	_R2_A[r] = 1 / R_A[r]**2 + 1 / (L - R_A[r])**2

plt.rcParams.update({
    "text.usetex": False,
    "font.family": "sans-serif"})

# plot chi vs R
# plt.figure()
# plt.xlim(0.8, R_half + 0.2)
# plt.yscale("log")
# plt.ylim(0.0, chi[1] + 0.1)
# plt.errorbar(R[1:], chi[1:], yerr=d_chi[1:], color="blue", marker='o', ms=5, mew=0.5, mfc='none', linestyle='none', linewidth=0.5, capsize=2.5, capthick=0.5)
# plt.plot(R_A, chi_A, color="blue", linewidth=0.5)
# plt.xlabel("$R$")
# plt.ylabel("$\chi(R,R)$")
# plt.savefig("plot_sigma_r2.pdf")
# plt.close()

# plot chi vs 1/R^2
plt.figure()
plt.xlim(0, 0.3)
# plt.yscale("log")
plt.ylim(0.0, chi[1] + 0.1)
plt.errorbar(_R2[1:], chi[1:], yerr=d_chi[1:], color="blue", marker='o', ms=5, mew=0.5, mfc='none', linestyle='none', linewidth=0.5, capsize=2.5, capthick=0.5)
plt.plot(_R2_A, chi_A, color="blue", linewidth=0.5)
plt.xlabel("$1/R^2$")
plt.ylabel("$\chi(R,R)$")
plt.savefig("plots/W_" + id + ".pdf")
plt.close()

# actual beta

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


# effective beta (16^4 * 5)

#   beta        sqrt(sigma)       beta_effective        <chi>
#   3.50           1.026              < 5.0             0.146
#   4.00           0.981              < 5.0             0.137
#   4.50           0.827              < 5.0             0.138

# strong

#   4.56           0.137
#   4.58           0.160
#   4.60           0.167
#   4.62           0.174
#   4.64           0.163
#   4.66           0.167
#   4.68           0.197
#   4.70           0.188                6.10            0.159
#   4.72           0.205                6.04
#   4.74           0.183                6.12
#	4.76           0.190                6.09
#	4.78           0.191                6.09
#   4.80           0.174                6.14
#   4.82
#   4.84
#	4.86
#	4.88
#   4.90           0.103   
#   5.00         < 0.1                > 6.3             0.134

# weak

# effective beta (20^4 * 3)

#   beta        sqrt(sigma)       beta_effective
#   3.00                              < 5.0
#   4.00          ~1.1                < 5.0
#   5.00           0.187                6.10
