#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 18 12:41:54 2021

@author: johnviljoen
"""
import numpy as np

import F16_sim
import gym
from gym import spaces

from parameters import act_lim, x_lim
from parameters import initial_state_vector_ft_rad as x0
from parameters import simulation_parameters as paras_sim

import progressbar
from utils import tic, toc, vis

# import scipy fmin for trim function
from scipy.optimize import minimize

class F16_env(gym.Env, F16_sim.F16_sim):
    
    def __init__(self, x0, paras_sim):
        super(F16_env, self).__init__(x0, paras_sim[4], paras_sim[3], paras_sim[0])
        
        # ignoring lef actuator limits as they are not directly commanded
        self.action_space = spaces.Box(low=np.array(act_lim[1])[0:4], high=np.array(act_lim[0])[0:4], dtype=np.float32)
        self.observation_space = spaces.Box(low=np.array(x_lim[1]), high=np.array(x_lim[0]), dtype=np.float32)
        self.measured_states = [2,6,7,8,9,10,11]
        self.t = 0
        
    def step(self, action):
        
        self.u = action
        
        # find xdot
        self.xdot = self.calc_xdot(self.x, self.u)
        
        # update x
        self.x += self.xdot*self.dt
        self.t += self.dt
        
        reward = 1
        isdone = 0
        
        return self.x, reward, isdone, {}
    
    def reset(self):
        
        self.x = x0
        self.u = x0[12:16]
        self.t = 0
    
    def _trim_obj_func(self, UX0, h, V):
              
        P3, dh, da, dr, alpha = UX0
        
        rho0 = 2.377e-3
        tfac = 1 - 0.703e-5*h
        
        temp = 519*tfac
        if h >= 35000:
            temp = 390
            
        rho = rho0*tfac**4.14
        qbar = 0.5*rho*V**2
        ps = 1715*rho*temp
        
        dlef = 1.38*alpha*180/np.pi - 9.05*qbar/ps + 1.45
        
        self.x = np.array([0, 0, h, 0, alpha, 0, V, alpha, 0, 0, 0, 0, P3, dh, da, dr, dlef, -alpha*180/np.pi])
        
        # set thrust limits
        if self.x[12] > self.action_space.high[0]:
            self.x[12] = self.action_space.high[0]
        elif self.x[12] < self.action_space.low[0]:
            self.x[12] = self.action_space.low[0]
           
        # set elevator limits
        if self.x[13] > self.action_space.high[1]:
            self.x[13] = self.action_space.high[1]
        elif self.x[13] < self.action_space.low[1]:
            self.x[13] = self.action_space.low[1]
            
        # set aileron limits
        if self.x[14] > self.action_space.high[2]:
            self.x[14] = self.action_space.high[2]
        elif self.x[14] < self.action_space.low[2]:
            self.x[14] = self.action_space.low[2]
            
        # set rudder limits
        if self.x[15] > self.action_space.high[3]:
            self.x[15] = self.action_space.high[3]
        elif self.x[15] < self.action_space.low[3]:
            self.x[15] = self.action_space.low[3]
            
        # set alpha limits
        if self.x[7] > self.observation_space.high[7]*np.pi/180:
            self.x[7] = self.observation_space.high[7]*np.pi/180
        elif self.x[7] < self.observation_space.low[7]*np.pi/180:
            self.x[7] = self.observation_space.low[7]*np.pi/180
            
        self.u = self.x[12:16]
        
        # calc xdot
        self.xdot = self.calc_xdot(self.x, self.u)
        
        phi_w = 10
        theta_w = 10
        psi_w = 10
        
        weight = np.array([0, 0, 5, phi_w, theta_w, psi_w, 2, 10, 10, 10, 10, 10]).transpose()
        
        cost = np.matmul(weight,self.xdot[0:12]**2)
        
        return cost
        
    def trim(self, h, V):
        
        # initial guesses
        thrust = 5000           # thrust, lbs
        elevator = -0.09        # elevator, degrees
        alpha = 8.49            # AOA, degrees
        rudder = -0.01          # rudder angle, degrees
        aileron = 0.01          # aileron, degrees
        
        UX0 = [thrust, elevator, alpha, rudder, aileron]
                
        opt = minimize(self._trim_obj_func, UX0, args=((h, V)), method='Nelder-Mead',tol=1e-10,options={'maxiter':5e+04})
        
        P3_t, dh_t, da_t, dr_t, alpha_t  = opt.x
        
        rho0 = 2.377e-3
        tfac = 1 - 0.703e-5*h
        
        temp = 519*tfac
        if h >= 35000:
            temp = 390
            
        rho = rho0*tfac**4.14
        qbar = 0.5*rho*V**2
        ps = 1715*rho*temp
        
        dlef = 1.38*alpha_t*180/np.pi - 9.05*qbar/ps + 1.45
        
        self.x = np.array([0, 0, h, 0, alpha_t, 0, V, alpha_t, 0, 0, 0, 0, P3_t, dh_t, da_t, dr_t, dlef, -alpha_t*180/np.pi])
        self.u = np.array([P3_t, dh_t, da_t, dr_t])
        
        return opt
    
    def linearise(self):
    
        eps = 1e-06
        
        A = np.zeros([len(self.x),len(self.x)])
        B = np.zeros([len(self.x),len(self.u)])
        C = np.zeros([len(self.measured_states),len(self.x)])
        D = np.zeros([len(self.measured_states),len(self.u)])
        
        # Perturb each of the state variables and compute linearization,
        for i in range(len(self.x)):
            
            dx = np.zeros((len(self.x),))
            dx[i] = eps
                        
            A[:, i] = (self.calc_xdot(self.x + dx, self.u) - self.calc_xdot(self.x, self.u)) / eps
            C[:, i] = np.array(list((self.calc_xdot(self.x + dx, self.u) - 
                                     self.calc_xdot(self.x, self.u))[i] for i in self.measured_states)) / eps
            
        # Perturb each of the input variables and compute linearization
        for i in range(len(self.u)):
            
            du = np.zeros((len(self.u),))
            du[i] = eps
                    
            B[:, i] = (self.calc_xdot(self.x, self.u + du) - self.calc_xdot(self.x, self.u)) / eps
            D[:, i] = np.array(list((self.calc_xdot(self.x, self.u + du) - 
                                     self.calc_xdot(self.x, self.u))[i] for i in self.measured_states))  / eps
        
        return A, B, C, D
    
# exit()
    
# model = F16_env(x0,paras_sim)
# A,B,C,D = model.linearise()

# model.trim(10000,700)

# A,B,C,D = model.linearise()

# #exit()

# rng = np.linspace(paras_sim[1], paras_sim[2], int((paras_sim[2]-paras_sim[1])/paras_sim[0]))
# output_vars = [6,7,8,9,10,11]

# # create storage
# x_storage = np.zeros([len(rng),len(model.x)])
# A = np.zeros([len(model.x),len(model.x),len(rng)])
# B = np.zeros([len(model.x),len(model.u),len(rng)])
# C = np.zeros([len(model.measured_states),len(model.x),len(rng)])
# D = np.zeros([len(model.measured_states),len(model.u),len(rng)])

# bar = progressbar.ProgressBar(maxval=len(rng)).start()

# tic()

# for idx, val in enumerate(rng):
    
#     #----------------------------------------#
#     #------------linearise model-------------#
#     #----------------------------------------#
    
#     # the object oriented linearisation function is about 10x slower than the functional programming one
#     [A[:,:,idx], B[:,:,idx], C[:,:,idx], D[:,:,idx]] = model.linearise()
#     #[A[:,:,idx], B[:,:,idx], C[:,:,idx], D[:,:,idx]] = linearise(model.x, model.u, model.measured_states, model.paras_sim[4], model._nlplant)
    
#     #----------------------------------------#
#     #--------------Take Action---------------#
#     #----------------------------------------#
    
#     # MPC prediction using squiggly C and M matrices
#     #CC, MM = calc_MC(paras_mpc[0], A[:,:,idx], B[:,:,idx], time_step)
    
    
#     #----------------------------------------#
#     #--------------Integrator----------------#
#     #----------------------------------------#    
    
#     x = model.step(model.u)
    
#     #----------------------------------------#
#     #------------Store History---------------#
#     #----------------------------------------#
    
#     x_storage[idx,:] = x
    
#     bar.update(idx)
    
# toc()
    
# vis(x_storage, rng)
