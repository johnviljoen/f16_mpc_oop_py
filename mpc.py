#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  7 14:55:50 2021

@author: johnviljoen
"""

import numpy as np

from sim import calc_xdot, calc_out

# In[]

def linearise(x, u, output_vars, fi_flag, nlplant):
    
    eps = 1e-06
    
    A = np.zeros([len(x),len(x)])
    B = np.zeros([len(x),len(u)])
    C = np.zeros([len(output_vars),len(x)])
    D = np.zeros([len(output_vars),len(u)])
    
    # Perturb each of the state variables and compute linearization
    for i in range(len(x)):
        
        dx = np.zeros((len(x),))
        dx[i] = eps
        
        A[:, i] = (calc_xdot(x + dx, u, fi_flag, nlplant) - calc_xdot(x, u, fi_flag, nlplant)) / eps
        C[:, i] = (calc_out(x + dx, u, output_vars) - calc_out(x, u, output_vars)) / eps
        
    # Perturb each of the input variables and compute linearization
    for i in range(len(u)):
        
        du = np.zeros((len(u),))
        du[i] = eps
                
        B[:, i] = (calc_xdot(x, u + du, fi_flag, nlplant) - calc_xdot(x, u, fi_flag, nlplant)) / eps
        D[:, i] = (calc_out(x, u + du, output_vars) - calc_out(x, u, output_vars)) / eps
    
    return A, B, C, D

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

def calc_HFG(A, B, C, hzn, Q, R):
    
    MM, CC = calc_MC(hzn, A, B, 1)

    Q = np.matmul(C.T,C)
    
    Q_full = dmom(Q, hzn)
    # Q_full = np.eye(hzn)
    
    R_full = np.eye(hzn) * 0.01
    
    H = np.matmul(np.matmul(CC.T, Q_full),CC) + R_full
    
    F = np.matmul(np.matmul(CC.T, Q_full), MM)
    
    G = np.matmul(np.matmul(MM.T, Q_full), MM)
    
    return H, F, G

# In[]

# dual mode predicted HFG
def calc_dm_HFG(A, B, C, K, hzn, Q, R):
    
    MM, CC = calc_MC(hzn, A, B, 1)

    Q = np.matmul(C.T,C)
    
    Q_full = dmom(Q, hzn)
    # Q_full = np.eye(hzn)
    
    rhs = Q + np.matmul(np.matmul(K.T,R), K)
    
    Qbar = np.array([])
    
    R_full = np.eye(hzn) * 0.01
    
    H = np.matmul(np.matmul(CC.T, Q_full),CC) + R_full
    
    F = np.matmul(np.matmul(CC.T, Q_full), MM)
    
    G = np.matmul(np.matmul(MM.T, Q_full), MM)
    
    return H, F, G