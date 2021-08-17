# Requirements

> Libraries used:

- ctypes
- os
- numpy
- matplotlib
- progressbar
- scipy
- osqp

> OS

- Linux only (due to the use of CDLL to call .so file, which Windows does not appear to accept, Mac may work but is untested)

# Description

The F16 class is contained in 'env.py', this contains the functions:

- step
- reset
- get_obs
- linearise
- trim
- calc_MPC_action
- _calc_xdot
- _calc_xdot_na
- _get_obs_na
 
This forms the heart of the simulation.

calc_MPC_action uses a number of functions currently found in 'utils.py', for matrix operations and reordering for input into the OSQP optimizer. The integrator used is a simple euler fixed time step.

The F16 class is instantiated in main.py as f16 as of time of writing. -> run main.py -> dir(f16), in a python console to examine the object.

All parameters are stored in 'parameters.py', stateVector, inputVector, and simulationParameter dataclasses are also instantiated here and used to instantiate the F16 class. These dataclasses contain things like the units, limits, names, initial conditions, of the state vector and input vector.

There are examples from notes found here: https://markcannon.github.io/assets/downloads/teaching/C21_Model_Predictive_Control/mpc_notes.pdf, that have been implemented with the same algorithm in 'notes_examples/example_2_x.py'.

# Credits

This simulation is based upon two works:

1. Nguyen "Simulator study of stall/post-stall characteristics of a fighter airplane with relaxed longitudinal static stability" - https://core.ac.uk/download/pdf/42866809.pdf

2. Stevens and Lewis "Aircraft Control and Simulation: Dynamics, Controls Design, and Autonomous Systems, 3rd Edition" - https://www.wiley.com/en-gb/Aircraft+Control+and+Simulation%3A+Dynamics%2C+Controls+Design%2C+and+Autonomous+Systems%2C+3rd+Edition-p-9781118870983

The former is the higher fidelity model, with accurate angles of attack between -20 -> +90 degrees, thanks to some crazy pilots and a lot of empirical testing. The latter is a reduced version of this model correct between angles of attacked between -10 -> +45 degrees. Both are accurate in a sideslip range between +-30 degrees.

The C implementation of both of these codes for MATLAB found here "https://dept.aem.umn.edu/~balas/darpa_sec/SEC.Software.html#F16Manual". This was originally designed to be a MEX file for MATLAB compilation and use, but I have cut those bits and am now using it as a shared object (.so) file to be called from Python.

This C implementation produces the instantaneous time derivatives of the 12 primary aircraft states (p5 of this document describes the IO: https://dept.aem.umn.edu/~balas/darpa_sec/software/F16Manual.pdf) these are:

x = {npos,epos,h,phi,theta,psi,Vt,alpha,beta,p,q,r,T,dh,da,dr,lf2,lf1}

# Technicalities

The C .so file has been compiled twice for two primary configurations of the F16 -> stable and unstable. These occur when xcg = 25, and 35 respectively. Should other parameters be changed the flags for gcc compilation of "nlplant.c" were:

gcc -fPIC -shared -lm -o nlplant.so nlplant.c

These two files are selected in the Python parameters.py using the "stab_flag" which is 1 for unstable and 0 for stable.

main.py calls the .so file using ctypes.CDLL. This allows access for individual functions in nlplant.so, two of which are called by main.py -> Nlplant and atmos. Nlplant is the main aforementioned one which generates the state derivates, and atmos is used for calculating the leading edge flap (lef) deflection.

