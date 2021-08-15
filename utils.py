#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  7 14:45:08 2021

@author: johnviljoen
"""
import time
import matplotlib.pyplot as plt
import numpy as np
from numpy import pi
import scipy
import ctypes

from parameters import u_lb, u_ub, udot_lb, udot_ub, x_lb, x_ub


# In[]

def setup_OSQP(A, B, Q, R, hzn, dt, x, act_states, x_lb, x_ub, u_lb, u_ub, udot_lb, udot_ub):
    
    """
    args:
        x:
            1D numpy array -> horizontal vector
        x_lb:
            2D numpy array -> vertical vector
        x_ub:
            2D numpy array -> vertical vector
        hzn:
            int
        MM:
            2D numpy array
        CC:
            2D numpy array
        u_lb:
            2D numpy array -> vertical vector
        u_ub:
            2D numpy array -> vertical vector
        udot_lb:
            2D numpy array -> vertical vector
        udot_ub:
            2D numpy array -> vertical vector
        u_len:
            int
            
    """
    
    # number of states, inputs = m, n
    m = len(x)
    n = len(act_states)
    
    MM, CC = calc_MC(A, B, dt, hzn)
    
    K = - dlqr(A, B, Q, R)
    Q_bar = scipy.linalg.solve_discrete_lyapunov((A + B @ K).T, Q + K.T @ R @ K)
    QQ = dmom(Q, hzn)
    QQ[-m:,-m:] = Q_bar
    
    RR = dmom(R, hzn)
    
    H = CC.T @ QQ @ CC + RR
    F = CC.T @ QQ @ MM
    G = MM.T @ QQ @ MM
    
    OSQP_P = 2*H
    OSQP_q = (2 * x[np.newaxis] @ F.T).T
    
    """
    There are three constraints to be enforced on the system:
        
        state constraints:
            x(n+1) = Ax(n) + Bu(n)
            
        input command limit constraints:
            u_min <= u <= u_max
            
        input command rate limit constraints:
            udot_min <= udot <= udot_max
    """
    
    # calculate state constraint limits vector
        
    x_lb = np.tile(x_lb,(hzn,1))    
    x_ub = np.tile(x_ub,(hzn,1))
    
    state_constr_lower = x_lb - MM @ x[np.newaxis].T
    state_constr_upper = x_ub - MM @ x[np.newaxis].T
    
    # the state constraint input sequence matrix is just CC
    
    # calculate the command saturation limits vector
    
    cmd_constr_lower = np.tile(u_lb,(hzn,1))
    cmd_constr_upper = np.tile(u_ub,(hzn,1))
    
    # calculate the command saturation input sequence matrix -> just eye
    
    cmd_constr_mat = np.eye(n*hzn)
    
    # calculate the command rate saturation limits vector
    
    u0_rate_constr_lower = act_states[np.newaxis].T + udot_lb
    u0_rate_constr_upper = act_states[np.newaxis].T + udot_ub
    
    cmd_rate_constr_lower = np.concatenate((u0_rate_constr_lower,np.tile(udot_lb,(hzn-1,1))))
    cmd_rate_constr_upper = np.concatenate((u0_rate_constr_upper,np.tile(udot_ub,(hzn-1,1))))
    
    # calculate the command rate saturation input sequence matrix
    
    cmd_rate_constr_mat = np.eye(n*hzn)
    for i in range(n*hzn):
        if i >= n:
            cmd_rate_constr_mat[i,i-n] = -1
            
    # assemble the complete matrices to send to OSQP
            
    OSQP_A = np.concatenate((CC, cmd_constr_mat, cmd_rate_constr_mat), axis=0)
            
    OSQP_l = np.concatenate((state_constr_lower, cmd_constr_lower, cmd_rate_constr_lower))
    
    OSQP_u = np.concatenate((state_constr_upper, cmd_constr_upper, cmd_rate_constr_upper))
    
    return OSQP_P, OSQP_q, OSQP_A, OSQP_l, OSQP_u

    

# In[ Calculate the MM, CC matrices -> squiggly M and squiggly C in the notes ]

def calc_MC(A, B, dt, hzn, includeFirstRow=False):
    
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
        MM[nstates*i:nstates*(i+1),:] = np.linalg.matrix_power(A,i+1) 
        for j in range(hzn):
            if i-j >= 0:
                CC[nstates*i:nstates*(i+1),ninputs*j:ninputs*(j+1)] = np.matmul(np.linalg.matrix_power(A,(i-j)),B)
            else:
                CC[nstates*i:nstates*(i+1),ninputs*j:ninputs*(j+1)] = Bz
                
    if includeFirstRow:
        CC = np.concatenate((np.zeros([nstates,ninputs*hzn]), CC), axis=0)
        MM = np.concatenate((np.eye(nstates), MM), axis=0)

    return MM, CC

# In[ OSQP setup functions ]

# def gen_rate_lim_constr_mat(u_len, hzn):
    
#     rate_lim_constr_mat = np.eye(u_len*hzn)
    
#     for i in range(u_len*hzn):
#         if i >= u_len:
#             rate_lim_constr_mat[i,i-u_len] = -1
            
#     return rate_lim_constr_mat

# def gen_rate_lim_constr_upper_lower(u_len, hzn, lower_limits, upper_limits):
    
#     rlcl = np.zeros([u_len*hzn,1])
#     rlcu = np.zeros([u_len*hzn,1])
    
#     rlcl[0:u_len,0] = -np.infty
#     rlcu[0:u_len,0] = np.infty
    
#     for i in range(hzn):
#         if i >= 1:
#             for j in range(u_len):
#                 rlcl[u_len*i+j,0] = lower_limits[j]
#                 rlcu[u_len*i+j,0] = upper_limits[j]
            
#     return rlcl, rlcu
    
# def gen_cmd_sat_constr_mat(u_len, hzn):
    
#     return dmom(np.eye(u_len), hzn)

# def gen_cmd_sat_constr_upper_lower(u_len, hzn, lower_limits, upper_limits):
    
#     cscl = np.zeros([u_len*hzn,1])
#     cscu = np.zeros([u_len*hzn,1])
    
#     for i in range(hzn):
#         for j in range(u_len):
#             cscl[u_len*i + j,0] = lower_limits[j]
#             cscu[u_len*i + j,0] = upper_limits[j]
            
#     return cscl, cscu

# def gen_OSQP_A(CC, cscm, rlcm):
#     return np.concatenate((CC, cscm, rlcm), axis=0)

# def setup_OSQP_paras(CC, A, x, hzn, ninputs, x_ub, x_lb, u_ub, u_lb, udot_ub, udot_lb):
    
#     """ 
#     args:
#         CC - numpy 2D array
#         A - numpy 2D array
#         x - numpy 2D array (vertical vector)
#         hzn - int
#         ninputs - int
#         x_ub - list
#         x_lb - list
#         u_ub - list
#         u_lb - list
#         udot_ub - list
#         udot_lb - list
        
#     returns:
        
#     """
    
#     cscm = gen_cmd_sat_constr_mat(ninputs, hzn)
#     rlcm = gen_rate_lim_constr_mat(ninputs, hzn)
#     OSQP_A = gen_OSQP_A(CC, cscm, rlcm)

#     x_ub = np.array(x_ub)[np.newaxis].T
#     x_lb = np.array(x_lb)[np.newaxis].T
    
#     u1 = np.concatenate(([x_ub - A @ x] * hzn))
#     l1 = np.concatenate(([x_lb - A @ x] * hzn))
    
#     cscl, cscu = gen_cmd_sat_constr_upper_lower(ninputs, hzn, u_lb, u_ub)
#     rlcl, rlcu = gen_rate_lim_constr_upper_lower(ninputs, hzn, udot_lb, udot_ub)
#     OSQP_l = np.concatenate((l1, cscl, rlcl))
#     OSQP_u = np.concatenate((u1, cscu, rlcu))
    
#     return OSQP_A, OSQP_l, OSQP_u
    



# In[]

def is_pos_def(x):
    return np.all(np.linalg.eigvals(x) > 0)

# In[ discrete linear quadratic regulator ]
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

# In[ degenerate 2D square matrix ]

def square_mat_degen_2d(mat, degen_idx):
    
    degen_mat = np.zeros([len(degen_idx),len(degen_idx)])
    
    for i in range(len(degen_idx)):
        
        degen_mat[:,i] = mat[degen_idx, [degen_idx[i] for x in range(len(degen_idx))]]
        
    return degen_mat

# In[ diagonal matrix of matrices ]

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

# In[ calculate actuator model time derivatives ]

def upd_lef(h, V, coeff, alpha, lef_state_1, lef_state_2, nlplant):
    
    nlplant.atmos(ctypes.c_double(h),ctypes.c_double(V),ctypes.c_void_p(coeff.ctypes.data))
    atmos_out = coeff[1]/coeff[2] * 9.05
    alpha_deg = alpha*180/pi
    
    LF_err = alpha_deg - (lef_state_1 + (2 * alpha_deg))
    #lef_state_1 += LF_err*7.25*time_step
    LF_out = (lef_state_1 + (2 * alpha_deg)) * 1.38
    
    lef_cmd = LF_out + 1.45 - atmos_out
    
    # command saturation
    lef_cmd = np.clip(lef_cmd,x_lb[16],x_ub[16])
    # rate saturation
    lef_err = np.clip((1/0.136) * (lef_cmd - lef_state_2),-25,25)
    
    return LF_err*7.25, lef_err

def upd_thrust(T_cmd, T_state):
    # command saturation
    T_cmd = np.clip(T_cmd,u_lb[0],u_ub[0])
    # rate saturation
    return np.clip(T_cmd - T_state, -10000, 10000)

def upd_dstab(dstab_cmd, dstab_state):
    # command saturation
    dstab_cmd = np.clip(dstab_cmd,u_lb[1],u_ub[1])
    # rate saturation
    return np.clip(20.2*(dstab_cmd - dstab_state), -60, 60)

def upd_ail(ail_cmd, ail_state):
    # command saturation
    ail_cmd = np.clip(ail_cmd,u_lb[2],u_ub[2])
    # rate saturation
    return np.clip(20.2*(ail_cmd - ail_state), -80, 80)

def upd_rud(rud_cmd, rud_state):
    # command saturation
    rud_cmd = np.clip(rud_cmd,u_lb[3],u_ub[3])
    # rate saturation
    return np.clip(20.2*(rud_cmd - rud_state), -120, 120)

# In[ replicate MATLAB tic toc functions ]

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

def bmatrix(a):
    """Returns a LaTeX bmatrix

    :a: numpy array
    :returns: LaTeX bmatrix as a string
    """
    if len(a.shape) > 2:
        raise ValueError('bmatrix can at most display two dimensions')
    lines = str(a).replace('[', '').replace(']', '').splitlines()
    rv = [r'\begin{bmatrix}']
    rv += ['  ' + ' & '.join(l.split()) + r'\\' for l in lines]
    rv +=  [r'\end{bmatrix}']
    return '\n'.join(rv)
    
# In[ visualise full 18 DoF system time history ]

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