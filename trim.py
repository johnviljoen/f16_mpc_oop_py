#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  7 14:49:02 2021

@author: johnviljoen
"""

import numpy as np
from numpy import pi

# import scipy fmin for trim function
from scipy.optimize import minimize

from sim import calc_xdot

from parameters import act_lim, x_lim

# calculate objective function for trimming
def obj_func(UX0, h_t, v_t, fi_flag, nlplant):
    
    V = v_t
    h = h_t
        
    P3, dh, da, dr, alpha = UX0
    
    npos = 0
    epos = 0
    #h
    phi = 0
    #theta = alpha in straight level flight
    psi = 0
    #V
    #alpha
    beta = 0
    p = 0
    q = 0
    r = 0
    #P3
    #dh
    #da
    #dr
    #dlef1
    #dlef2
    
    rho0 = 2.377e-3
    tfac = 1 - 0.703e-5*h
    
    temp = 519*tfac
    if h >= 35000:
        temp = 390
        
    rho = rho0*tfac**4.14
    qbar = 0.5*rho*V**2
    ps = 1715*rho*temp
    
    dlef = 1.38*alpha*180/pi - 9.05*qbar/ps + 1.45
    
    x = np.array([npos, epos, h, phi, alpha, psi, V, alpha, beta, p, q, r, P3, dh, da, dr, dlef, -alpha*180/pi])
    
    # set thrust limits
    if x[12] > act_lim[0][0]:
        x[12] = act_lim[0][0]
    elif x[12] < act_lim[1][0]:
        x[12] = act_lim[1][0]
       
    # set elevator limits
    if x[13] > act_lim[0][1]:
        x[13] = act_lim[0][1]
    elif x[13] < act_lim[1][1]:
        x[13] = act_lim[1][1]
        
    # set aileron limits
    if x[14] > act_lim[0][2]:
        x[14] = act_lim[0][2]
    elif x[14] < act_lim[1][2]:
        x[14] = act_lim[1][2]
        
    # set rudder limits
    if x[15] > act_lim[0][3]:
        x[15] = act_lim[0][3]
    elif x[15] < act_lim[1][3]:
        x[15] = act_lim[1][3]
        
    # set alpha limits
    if x[7] > x_lim[0][7]*pi/180:
        x[7] = x_lim[0][7]*pi/180
    elif x[7] < x_lim[1][7]*pi/180:
        x[7] = x_lim[1][7]*pi/180
        
    u = np.array([x[12],x[13],x[14],x[15]])
    
    xdot = calc_xdot(x, u, fi_flag, nlplant)
    
    phi_w = 10
    theta_w = 10
    psi_w = 10
    
    weight = np.array([0, 0, 5, phi_w, theta_w, psi_w, 2, 10, 10, 10, 10, 10]).transpose()
    
    cost = np.matmul(weight,xdot[0:12]**2)
    
    return cost

def trim(h_t, v_t, fi_flag, nlplant):
    
    # initial guesses
    thrust = 5000           # thrust, lbs
    elevator = -0.09        # elevator, degrees
    alpha = 8.49            # AOA, degrees
    rudder = -0.01          # rudder angle, degrees
    aileron = 0.01          # aileron, degrees
    
    UX0 = [thrust, elevator, alpha, rudder, aileron]
            
    opt = minimize(obj_func, UX0, args=((h_t, v_t, fi_flag, nlplant)), method='Nelder-Mead',tol=1e-10,options={'maxiter':5e+04})
    
    P3_t, dstab_t, da_t, dr_t, alpha_t  = opt.x
    
    rho0 = 2.377e-3
    tfac = 1 - 0.703e-5*h_t
    
    temp = 519*tfac
    if h_t >= 35000:
        temp = 390
        
    rho = rho0*tfac**4.14
    qbar = 0.5*rho*v_t**2
    ps = 1715*rho*temp
    
    dlef = 1.38*alpha_t*180/pi - 9.05*qbar/ps + 1.45
    
    x_trim = np.array([0, 0, h_t, 0, alpha_t, 0, v_t, alpha_t, 0, 0, 0, 0, P3_t, dstab_t, da_t, dr_t, dlef, -alpha_t*180/pi])
    
    return x_trim, opt