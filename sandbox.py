#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  6 08:06:11 2021

@author: johnviljoen
"""

from scipy.optimize import minimize
from scipy.signal import cont2discrete
import numpy as np
import scipy

from parameters import initial_state_vector_ft_rad as x0
from parameters import simulation_parameters as paras_sim
from parameters import model_predictive_control_parameters as paras_mpc

from env import F16

from scipy.signal import cont2discrete

from parameters import simulation_parameters as paras_sim
from parameters import initial_state_vector_ft_rad as x0

from utils import tic, toc, vis, dmom, square_mat_degen_2d, gen_rate_lim_constr_mat, \
    gen_rate_lim_constr_upper_lower, gen_cmd_sat_constr_mat, gen_cmd_sat_constr_upper_lower, \
        gen_OSQP_A, dlqr, calc_MC

from scipy.sparse import csc_matrix

import osqp

# In[ make OSQP_A ]

f16 = F16(x0, x0[12:16], paras_sim)

A,B,C,D = f16.linearise(f16.x, f16.u)
A,B,C,D = cont2discrete((A,B,C,D), f16.dt)[0:4]

A_degen = square_mat_degen_2d(A, f16.x_degen_idx)
B_degen = B[f16.x_degen_idx,1:4]

Q = square_mat_degen_2d(C.T @ C, f16.x_degen_idx)  

R = np.eye(len(f16.u_degen))*1000

lim = np.copy(f16.lim)
x_degen_idx = np.copy(f16.x_degen_idx)
x_degen = np.copy(f16.x_degen)

dt = 0.001

u = np.copy(f16.u_degen)

hzn = 10

def calc_MPC_action(u, hzn, A, B, dt, Q, R, x_idx, x, lim):

    cscm = gen_cmd_sat_constr_mat(u, hzn)
    
    rlcm = gen_rate_lim_constr_mat(u,hzn)
    
    #################################
    
    
    MM, CC = calc_MC(hzn, A, B, dt)
    
    # CCi = CC[]
    
    OSQP_A = gen_OSQP_A(CC, cscm, rlcm)
    
    # In[ make OSQP_l and OSQP_u ]
    
    K = dlqr(A, B, Q, R)
    
    Q_bar = scipy.linalg.solve_discrete_lyapunov((A + np.matmul(B, K)).T, Q + np.matmul(np.matmul(K.T,R), K))
    
    QQ = dmom(Q,hzn)
    RR = dmom(R,hzn)
    
    QQ[-A.shape[0]:,-A.shape[0]:] = Q_bar
    
    H = CC.T @ QQ @ CC + RR
    F = CC.T @ QQ @ MM
    G = MM.T @ QQ @ MM
    
    # u_seq = np.concatenate((u,u,u,u,u,u,u,u,u,u),axis=0)
    
    x = np.copy(x)
    
    x_u = np.array(list(list(lim.flatten(order='F'))[i] for i in x_idx))[np.newaxis].T
    x_l = np.array(list(list(lim.flatten(order='F'))[i] for i in [i + 17 for i in x_idx]))[np.newaxis].T
    
    u1 = np.concatenate(([x_u - A @ x] * hzn))
    l1 = np.concatenate(([x_l - A @ x] * hzn))
    
    cscl, cscu = gen_cmd_sat_constr_upper_lower(u, hzn, lim[13:16,1], lim[13:16,0])
    
    rlcl, rlcu = gen_rate_lim_constr_upper_lower(u, hzn, [-60, -80, -120], [60, 80, 120])
    
    OSQP_l = np.concatenate((l1, cscl, rlcl))
    OSQP_u = np.concatenate((u1, cscu, rlcu))
    
    # In[]
    
    
    P = 2*H
    q = (2 * x.T @ F.T).T
    
    m = osqp.OSQP()
    
    m.setup(P=csc_matrix(P), q=q, A=csc_matrix(OSQP_A), l=OSQP_l, u=OSQP_u, max_iter=40000, verbose=False)
    
    res = m.solve()
    
    u_star = res.x
    
    return u_star, H, F, G

# In[]
u_star, H, F, G = calc_MPC_action(u, hzn, A_degen, B_degen, dt, Q, R, x_degen_idx, x_degen, lim)
