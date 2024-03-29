# su2_x

SU(2) pure gauge theory lattice model with an extra dimension ("x" for extra). The extra dimension uses open boundary conditions and all links are set to the identity. So essentially we have an array of lattice "slices" which coupled in the extra dimenension. Observables are measured on each slice separately. The extra dimension is only taken into account when doing a lattice update. Hybrid Monte Carlo is normally used for updates, but there is also a heatbath algorithm that can be used. The coupling for plaquettes in the extra dimension is the normal coupling times a factor epsilon_5 (e.g. with beta = 5.8 and epsilon_5 = 0.5, the coupling in the extra dimension will be beta_5 = 2.9).

The hybrid Monte Carlo algorithm uses an adaptive step size during thermalization to target an Metropolis acceptance rate between 70% and 90%. The step size gets locked in after thermalization. Also, the Metropolis acceptance is ignored for the first half of thermalization to speed up the initial thermalization updates.

The lattice code uses the C++ thread library to parallelize the following:

- Update methods (HMC and heat bath)
- Smearing (Wilson flow, cooling, stout)
- Observable measurements (Wilson loop, topological charge)

Most of these algorithms run on each time slice in parallel, so the optimimum number of threads to use is the number of time slices. However, the heat bath and lattice cooling algorithms run the even and odd time slices separately in order to avoid multiple threads modifying a link at the same time.

## su2_x_test

The su2_x_test.cpp program is a test I used to measure the average plaquette at different values of beta. The output is written to the console.

usage: su2_x_test N T N5 D beta eps5 n_data

- N: lattice size in spatial directions
- T: lattice size in time direction
- N5: number of slices in extra dimension
- D: number of dimensions (doesn't include extra dimension)
- beta_start: first beta value
- beta_stop: last beta value
- beta_inc: beta increment
- eps5: coupling factor in extra dimension
- n_data: number of data points to measure at each value of beta
