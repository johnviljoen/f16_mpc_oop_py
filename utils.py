#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  7 14:45:08 2021

@author: johnviljoen
"""


# import time for tic toc functions
import time

# import matplotlib for visualisation
import matplotlib.pyplot as plt

import numpy as np
from numpy import pi

import scipy

# In[ All expect u in vertical vector format, hzn as a scalar ]

def gen_rate_lim_constr_mat(u, hzn):
    
    rate_lim_constr_mat = np.eye(len(u)*hzn)
    
    for i in range(len(u)*hzn):
        if i >= len(u):
            rate_lim_constr_mat[i,i-len(u)] = -1
            
    return rate_lim_constr_mat

def gen_rate_lim_constr_upper_lower(u, hzn, lower_limits, upper_limits):
    
    rlcl = np.zeros([len(u)*hzn,1])
    rlcu = np.zeros([len(u)*hzn,1])
    
    rlcl[0:len(u),0] = -np.infty
    rlcu[0:len(u),0] = np.infty
    
    for i in range(hzn):
        if i >= 1:
            for j in range(len(u)):
                rlcl[len(u)*i+j,0] = lower_limits[j]
                rlcu[len(u)*i+j,0] = upper_limits[j]
            
    return rlcl, rlcu
    
def gen_cmd_sat_constr_mat(u, hzn):
    
    return dmom(np.eye(len(u)), hzn)

def gen_cmd_sat_constr_upper_lower(u, hzn, lower_limits, upper_limits):
    
    cscl = np.zeros([len(u)*hzn,1])
    cscu = np.zeros([len(u)*hzn,1])
    
    for i in range(hzn):
        for j in range(len(u)):
            cscl[len(u)*i + j,0] = lower_limits[j]
            cscu[len(u)*i + j,0] = upper_limits[j]
            
    return cscl, cscu

def gen_OSQP_A(CC, cscm, rlcm):
    return np.concatenate((CC, cscm, rlcm), axis=0)
    

# In[]
def calc_MC(hzn, A, B, dt):

    # hzn is the horizon
    nstates = A.shape[0]
    ninputs = B.shape[1]
    
    # x0 is the initial state vector of shape (nstates, 1)
    # u is the matrix of input vectors over the course of the prediction of shape (ninputs,horizon)
    
    # initialise CC, MM, Bz
    CC = np.zeros([nstates*hzn, ninputs*hzn])
    MM = np.zeros([nstates*hzn, nstates])
    Bz = np.zeros([nstates, ninputs])
    
    for i in range(hzn):
        MM[nstates*i:nstates*(i+1),:] = np.linalg.matrix_power(A,i+1) * dt ** (i+1)
        for j in range(hzn):
            if i-j >= 0:
                CC[nstates*i:nstates*(i+1),ninputs*j:ninputs*(j+1)] = np.matmul(np.linalg.matrix_power(A,(i-j)),B) * dt ** (i-j+1)
            else:
                CC[nstates*i:nstates*(i+1),ninputs*j:ninputs*(j+1)] = Bz

    return MM, CC

# In[]
# from https://github.com/python-control/python-control/issues/359:
def dlqr(A,B,Q,R):
    """
    Solve the discrete time lqr controller.
    x[k+1] = A x[k] + B u[k]
    cost = sum x[k].T*Q*x[k] + u[k].T*R*u[k]
    
    
    Discrete-time Linear Quadratic Regulator calculation.
    State-feedback control  u[k] = -K*(x_ref[k] - x[k])
    select the states that you want considered and make x[k] the difference
    between the current x and the desired x.
      
    How to apply the function:    
        K = dlqr(A_d,B_d,Q,R)
      
    Inputs:
      A_d, B_d, Q, R  -> all numpy arrays  (simple float number not allowed)
      
    Returns:
      K: state feedback gain
    
    """
    # first, solve the ricatti equation
    P = np.matrix(scipy.linalg.solve_discrete_are(A, B, Q, R))
    # compute the LQR gain
    K = np.matrix(scipy.linalg.inv(B.T @ P @ B+R) @ (B.T @ P @ A))
    return K

# In[]

def square_mat_degen_2d(mat, degen_idx):
    
    degen_mat = np.zeros([len(degen_idx),len(degen_idx)])
    
    for i in range(len(degen_idx)):
        
        degen_mat[:,i] = mat[degen_idx, [degen_idx[i] for x in range(len(degen_idx))]]
        
    return degen_mat

# In[]

# def gmim(mat, submat_dims):
#     # get matrix in matrix -> gmim
    
#     mat[]

def dmom(mat, num_mats):
    # diagonal matrix of matrices -> dmom
    
    # dimension extraction
    nrows = mat.shape[0]
    ncols = mat.shape[1]
    
    # matrix of matrices matomats -> I thought it sounded cool
    matomats = np.zeros((nrows*num_mats,ncols*num_mats))
    
    for i in range(num_mats):
        for j in range(num_mats):
            if i == j:
                matomats[nrows*i:nrows*(i+1),ncols*j:ncols*(j+1)] = mat
                
    return matomats

# In[]

def TicTocGenerator():
    # Generator that returns time differences
    ti = 0           # initial time
    tf = time.time() # final time
    while True:
        ti = tf
        tf = time.time()
        yield tf-ti # returns the time difference

TicToc = TicTocGenerator() # create an instance of the TicTocGen generator

# This will be the main function through which we define both tic() and toc()
def toc(tempBool=True):
    # Prints the time difference yielded by generator instance TicToc
    tempTimeInterval = next(TicToc)
    if tempBool:
        print( "Elapsed time: %f seconds.\n" %tempTimeInterval )

def tic():
    # Records a time in TicToc, marks the beginning of a time interval
    toc(False)
    
# In[]

def vis(x_storage, rng):

    fig, axs = plt.subplots(12, 1)
    #fig.suptitle('Vertically stacked subplots')
    axs[0].plot(rng, x_storage[:,0])
    axs[0].set_ylabel('npos (ft)')
    
    axs[1].plot(rng, x_storage[:,1])
    axs[1].set_ylabel('epos (ft)')
    
    axs[2].plot(rng, x_storage[:,2])
    axs[2].set_ylabel('h (ft)')
    
    axs[3].plot(rng, x_storage[:,3])
    axs[3].set_ylabel('$\phi$ (rad)')
    
    axs[4].plot(rng, x_storage[:,4])
    axs[4].set_ylabel('$\theta$ (rad)')
    
    axs[5].plot(rng, x_storage[:,5])
    axs[5].set_ylabel('$\psi$ (rad)')
    
    axs[6].plot(rng, x_storage[:,6])
    axs[6].set_ylabel("V_t (ft/s)")
    
    axs[7].plot(rng, x_storage[:,7]*180/pi)
    axs[7].set_ylabel('alpha (deg)')
    
    axs[8].plot(rng, x_storage[:,8]*180/pi)
    axs[8].set_ylabel('beta (deg)')
    
    axs[9].plot(rng, x_storage[:,9]*180/pi)
    axs[9].set_ylabel('p (deg/s)')
    
    axs[10].plot(rng, x_storage[:,10]*180/pi)
    axs[10].set_ylabel('q (deg/s)')
    
    axs[11].plot(rng, x_storage[:,11]*180/pi)
    axs[11].set_ylabel('r (deg/s)')
    axs[11].set_xlabel('time (s)')
    
    fig2, axs2 = plt.subplots(5,1)
    
    axs2[0].plot(rng, x_storage[:,12])
    axs2[0].set_ylabel('P3')
    
    axs2[1].plot(rng, x_storage[:,13])
    axs2[1].set_ylabel('dh')
    
    axs2[2].plot(rng, x_storage[:,14])
    axs2[2].set_ylabel('da')
    
    axs2[3].plot(rng, x_storage[:,15])
    axs2[3].set_ylabel('dr')
    
    axs2[4].plot(rng, x_storage[:,16])
    axs2[4].set_ylabel('lef')