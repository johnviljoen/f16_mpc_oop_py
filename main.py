#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 31 14:41:39 2021

@author: johnviljoen
"""

# dependencies
import numpy as np
from numpy import pi
import gym
from gym import spaces
from scipy.optimize import minimize
import ctypes
from ctypes import CDLL
import os
import progressbar

# custom files
from parameters import act_lim, x_lim
from parameters import initial_state_vector_ft_rad as x0
from parameters import simulation_parameters as paras_sim
from utils import tic, toc, vis

class F16(gym.Env):
    
    def __init__(self, x0, u0, paras_sim):
        
        super().__init__()
        
        # system state
        self.x = np.copy(x0[np.newaxis].T)
        self.x0 = np.copy(x0[np.newaxis].T)
        # input demand
        self.u = np.copy(u0[np.newaxis].T)
        self.u0 = np.copy(u0[np.newaxis].T)
        # output state indices
        self.y_vars = [6,7,8,9,10,11]
        # measured state indices
        self.z_vars = [6,7,8,9]
        # fidelity flGag
        self.fi_flag = paras_sim[4]
        # time step
        self.dt = paras_sim[0]
        
        self.xdot = np.zeros([x0.shape[0]])
        
        # create interface with c shared library .so file in folder "C"
        if paras_sim[3] == 1:
            so_file = os.getcwd() + "/C/nlplant_xcg35.so"
        elif paras_sim[3] == 0:
            so_file = os.getcwd() + "/C/nlplant_xcg25.so"
        nlplant = CDLL(so_file)
        self.nlplant = nlplant
        
        self.action_space = spaces.Box(low=np.array(act_lim[1])[0:4], high=np.array(act_lim[0])[0:4], dtype=np.float32)
        self.observation_space = spaces.Box(low=np.array(x_lim[1] + act_lim[1]), high=np.array(x_lim[0] + act_lim[0]), shape=(17,), dtype=np.float32)
        
    def calc_xdot(self, x, u):
        
        def upd_thrust(T_cmd, T_state):
            # command saturation
            T_cmd = np.clip(T_cmd,act_lim[1][0],act_lim[0][0])
            # rate saturation
            return np.clip(T_cmd - T_state, -10000, 10000)
        
        def upd_dstab(dstab_cmd, dstab_state):
            # command saturation
            dstab_cmd = np.clip(dstab_cmd,act_lim[1][1],act_lim[0][1])
            # rate saturation
            return np.clip(20.2*(dstab_cmd - dstab_state), -60, 60)
        
        def upd_ail(ail_cmd, ail_state):
            # command saturation
            ail_cmd = np.clip(ail_cmd,act_lim[1][2],act_lim[0][2])
            # rate saturation
            return np.clip(20.2*(ail_cmd - ail_state), -80, 80)
        
        def upd_rud(rud_cmd, rud_state):
            # command saturation
            rud_cmd = np.clip(rud_cmd,act_lim[1][3],act_lim[0][3])
            # rate saturation
            return np.clip(20.2*(rud_cmd - rud_state), -120, 120)
        
        def upd_lef(h, V, coeff, alpha, lef_state_1, lef_state_2, nlplant):
            
            nlplant.atmos(ctypes.c_double(h),ctypes.c_double(V),ctypes.c_void_p(coeff.ctypes.data))
            atmos_out = coeff[1]/coeff[2] * 9.05
            alpha_deg = alpha*180/pi
            
            LF_err = alpha_deg - (lef_state_1 + (2 * alpha_deg))
            #lef_state_1 += LF_err*7.25*time_step
            LF_out = (lef_state_1 + (2 * alpha_deg)) * 1.38
            
            lef_cmd = LF_out + 1.45 - atmos_out
            
            # command saturation
            lef_cmd = np.clip(lef_cmd,act_lim[1][4],act_lim[0][4])
            # rate saturation
            lef_err = np.clip((1/0.136) * (lef_cmd - lef_state_2),-25,25)
            
            return LF_err*7.25, lef_err
        
        # initialise variables
        xdot = np.zeros([18,1])
        temp = np.zeros(6)
        coeff = np.zeros(3)
        #--------------Thrust Model--------------#
        temp[0] = upd_thrust(u[0], x[12])
        #--------------Dstab Model---------------#
        temp[1] = upd_dstab(u[1], x[13])
        #-------------aileron model--------------#
        temp[2] = upd_ail(u[2], x[14])
        #--------------rudder model--------------#
        temp[3] = upd_rud(u[3], x[15])
        #--------leading edge flap model---------#
        temp[5], temp[4] = upd_lef(x[2], x[6], coeff, x[7], x[17], x[16], self.nlplant)
        #----------run nlplant for xdot----------#
        self.nlplant.Nlplant(ctypes.c_void_p(x.ctypes.data), ctypes.c_void_p(xdot.ctypes.data), ctypes.c_int(self.fi_flag))    
        #----------assign actuator xdots---------#
        xdot[12:18,0] = temp
        return xdot
        
    def step(self, action):
        self.x += self.calc_xdot(self.x, self.u)*self.dt
        return self.x
    
    def reset(self):
        self.x = np.copy(self.x0)
        self.u = np.copy(self.u0)
        
    def get_obs(self, x, u):
        return x[self.y_vars]
        
    def trim(self, h_t, v_t):
        
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
            
            xdot = self.calc_xdot(x, u)
            
            phi_w = 10
            theta_w = 10
            psi_w = 10
            
            weight = np.array([0, 0, 5, phi_w, theta_w, psi_w, 2, 10, 10, 10, 10, 10]).transpose()
            
            cost = np.matmul(weight,xdot[0:12]**2)
            
            return cost
        
        # initial guesses
        thrust = 5000           # thrust, lbs
        elevator = -0.09        # elevator, degrees
        alpha = 8.49            # AOA, degrees
        rudder = -0.01          # rudder angle, degrees
        aileron = 0.01          # aileron, degrees
        
        UX0 = [thrust, elevator, alpha, rudder, aileron]
                
        opt = minimize(obj_func, UX0, args=((h_t, v_t, self.fi_flag, self.nlplant)), method='Nelder-Mead',tol=1e-10,options={'maxiter':5e+04})
        
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
        
    def linearise(self, x, u):
        
        eps = 1e-06
        
        A = np.zeros([len(x),len(x)])
        B = np.zeros([len(x),len(u)])
        C = np.zeros([len(self.y_vars),len(x)])
        D = np.zeros([len(self.y_vars),len(u)])
        
        # Perturb each of the state variables and compute linearization
        for i in range(len(x)):
            
            dx = np.zeros([len(x),1])
            dx[i] = eps
            
            A[:, i] = (self.calc_xdot(x + dx, u)[:,0] - self.calc_xdot(x, u)[:,0]) / eps
            C[:, i] = (self.get_obs(x + dx, u)[:,0] - self.get_obs(x, u)[:,0]) / eps
            
        # Perturb each of the input variables and compute linearization
        for i in range(len(u)):
            
            du = np.zeros([len(u),1])
            du[i] = eps
                    
            B[:, i] = (self.calc_xdot(x, u + du)[:,0] - self.calc_xdot(x, u)[:,0]) / eps
            D[:, i] = (self.get_obs(x, u + du)[:,0] - self.get_obs(x, u)[:,0]) / eps
        
        return A, B, C, D        
        
# make starting array immutable to cause error if used innapropriately
x0.flags.writeable = False

# instantiate the object
f16 = F16(x0, x0[12:16], paras_sim)

# trim the aircraft at 10000ft, 700 ft/s
f16.x = f16.trim(10000,700)[0][np.newaxis].T
f16.u = f16.x[12:16]


rng = np.linspace(paras_sim[1], paras_sim[2], int((paras_sim[2]-paras_sim[1])/paras_sim[0]))

# create storage
x_storage = np.zeros([len(rng),len(f16.x)])
A = np.zeros([len(f16.x),len(f16.x),len(rng)])
B = np.zeros([len(f16.x),len(f16.u),len(rng)])
C = np.zeros([len(f16.y_vars),len(f16.x),len(rng)])
D = np.zeros([len(f16.y_vars),len(f16.u),len(rng)])



bar = progressbar.ProgressBar(maxval=len(rng)).start()

tic()

for idx, val in enumerate(rng):
    
    #----------------------------------------#
    #------------linearise model-------------#
    #----------------------------------------#
    
    [A[:,:,idx], B[:,:,idx], C[:,:,idx], D[:,:,idx]] = f16.linearise(f16.x, f16.u)
    
    #----------------------------------------#
    #--------------Take Action---------------#
    #----------------------------------------#
    
    # MPC prediction using squiggly C and M matrices
    #CC, MM = calc_MC(paras_mpc[0], A[:,:,idx], B[:,:,idx], time_step)
    
    
    #----------------------------------------#
    #--------------Integrator----------------#
    #----------------------------------------#    
    
    x = f16.step(f16.u)
    
    #----------------------------------------#
    #------------Store History---------------#
    #----------------------------------------#
    
    x_storage[idx,:] = x[:,0]
    
    bar.update(idx)

toc()

# In[]

#----------------------------------------------------------------------------#
#---------------------------------Visualise----------------------------------#
#----------------------------------------------------------------------------#

#%matplotlib qt

vis(x_storage, rng)

