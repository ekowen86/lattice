import sqlite3
import numpy as np
import matplotlib.pyplot as plt

id = "16_1_6300"
L = int(16) # lattice size
T = int(16) # lattice size (time direction)

def jackknife_mean(data, n_sub):
	f_sub = float(n_sub) # number of data subsets
	start = 0 # start of current block of data
	count = len(data) # number of data points
	size = count / n_sub # size of each block
	
	a_bar = np.mean(data)
	d_a = 0.0
	
	for n in range(0, n_sub):
		end = np.minimum(start + size, count)
		r = range(start, end)
		
		# copy the array and delete the current subset
		data_d = np.delete(data, r)
		
		# calculate the creutz ratio and add to the error
		d_a += (np.mean(data_d) - a_bar)**2.0

		start = end

	d_a = np.sqrt((f_sub - 1) / f_sub * d_a)
	return (a_bar, d_a)


print("Opening sqlite file...")
con = sqlite3.connect("data/su3_x_" + id + ".db")
cursor = con.cursor()
cursor.execute("SELECT Q FROM topocharge WHERE n_c='1000'")
rows = cursor.fetchall()

print("Converting SQL data to numpy array...")
def row2float(row):
	return float(row[0])

Q = np.array(map(row2float, rows))
absQ = np.abs(Q)
Q2 = Q * Q

V = L**3 * T

Q_bar, d_Q = jackknife_mean(Q, 20)
absQ_bar, d_absQ = jackknife_mean(absQ, 20)
Q2_bar, d_Q2 = jackknife_mean(Q2, 20)
chi = (Q2_bar / V)**(0.25)
d_chi = 0.5 * (1 / Q2_bar / V)**(0.25) * d_Q2

print("<Q> = %.12f (%.12f)" % (Q_bar, d_Q))
print("<|Q|> = %.12f (%.12f)" % (absQ_bar, d_absQ))
print("<Q^2> = %.12f (%.12f)" % (Q2_bar, d_Q2))
print("<chi> = %.12f (%.12f)" % (chi, d_chi))

N_best = 1.08
Q_N = np.round(N_best * Q)

# calculate average change in Q per configuration
DQ_N = np.empty(len(Q_N) - 1)
for i in range(0, len(DQ_N)):
	DQ_N[i] = np.abs(Q_N[i] - Q_N[i+1])

DQ_bar, d_DQ = jackknife_mean(DQ_N, 20)
DQnorm = DQ_bar / absQ_bar
d_DQnorm = d_DQ / absQ_bar + DQ_bar * d_absQ / absQ_bar**2

print("<DQ> = %.12f (%.12f)" % (DQ_bar, d_DQ))
print("<DQ>/<|Q|> = %.12f (%.12f)" % (DQnorm, d_DQnorm))


begin = 0.0
end = 2.0
N_array = np.linspace(begin,end,100)
N_sum = np.zeros(len(N_array))

for i in range(0, len(N_array)):
	for j in range(0, len(Q)):
		NQ = N_array[i] * Q[j]
		N_sum[i] += (NQ - np.round(NQ)) ** 2.0
		
print("Plotting integer function...")
plt.figure()
plt.plot(N_array, N_sum)
plt.xlabel("N")
plt.ylabel("f(N)")
plt.savefig("plots/q_int.pdf")

begin = 1.00
end = 1.15
N_array = np.linspace(begin,end,100)
N_sum = np.zeros(len(N_array))
for i in range(0, len(N_array)):
	for j in range(0, len(Q)):
		NQ = N_array[i] * Q[j]
		N_sum[i] += (NQ - np.round(NQ)) ** 2.0
plt.figure()
plt.plot(N_array, N_sum)
# plt.xlim()
plt.xlabel("N")
plt.ylabel("f(N)")
plt.savefig("plots/q_int_zoom.pdf")

# print(Q * N_best)
print(Q_N)

print("Plotting histogram...")
plt.figure()
plt.xlim(-6.5, 6.5)
plt.hist(Q_N, 9, (-4.5,4.5), color="orange")
plt.xlabel("Q")
plt.ylabel("Count")
plt.savefig("plots/Q_" + id + ".pdf")

