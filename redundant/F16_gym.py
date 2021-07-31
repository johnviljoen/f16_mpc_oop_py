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

import os
from ctypes import CDLL
import ctypes

from stable_baselines3 import A2C
from stable_baselines3.common.vec_env import DummyVecEnv

from stable_baselines3.common.env_checker import check_env

from numba import jit

class F16_env(gym.Env):
    
    metadata = {'render.modes': ['human']}
    
    def __init__(self, x0, paras_sim):
        super(F16_env, self).__init__()
        
        self.x = x0
        self.u = x0[12:16]
        self.fi_flag = paras_sim[4]
        self.dt = paras_sim[0]
        self.xdot = np.zeros([x0.shape[0]])
        
        # create interface with c shared library .so file in folder "C"
        if paras_sim[3] == 1:
            so_file = os.getcwd() + "/C/nlplant_xcg35.so"
        elif paras_sim[3] == 0:
            so_file = os.getcwd() + "/C/nlplant_xcg25.so"
        nlplant = CDLL(so_file)
        self.nlplant = nlplant
        
        # ignoring lef actuator limits as they are not directly commanded
        self.action_space = spaces.Box(low=np.array(act_lim[1])[0:4], high=np.array(act_lim[0])[0:4], dtype=np.float32)
        self.observation_space = spaces.Box(low=np.array(x_lim[1] + act_lim[1]), high=np.array(x_lim[0] + act_lim[0]), shape=(17,), dtype=np.float32)
        self.measured_states = [2,6,7,8,9,10,11]
        self.t = 0
        
    def upd_thrust(self, T_cmd, T_state):
        # command saturation
        T_cmd = np.clip(T_cmd,act_lim[1][0],act_lim[0][0])
        # rate saturation
        T_err = np.clip(T_cmd - T_state, -10000, 10000)
        return T_err
    
    def upd_dstab(self, dstab_cmd, dstab_state):
        # command saturation
        dstab_cmd = np.clip(dstab_cmd,act_lim[1][1],act_lim[0][1])
        # rate saturation
        dstab_err = np.clip(20.2*(dstab_cmd - dstab_state), -60, 60)
        return dstab_err
    
    def upd_ail(self, ail_cmd, ail_state):
        # command saturation
        ail_cmd = np.clip(ail_cmd,act_lim[1][2],act_lim[0][2])
        # rate saturation
        ail_err = np.clip(20.2*(ail_cmd - ail_state), -80, 80)
        return ail_err
    
    def upd_rud(self, rud_cmd, rud_state):
        # command saturation
        rud_cmd = np.clip(rud_cmd,act_lim[1][3],act_lim[0][3])
        # rate saturation
        rud_err = np.clip(20.2*(rud_cmd - rud_state), -120, 120)
        return rud_err
    
    def upd_lef(self, h, V, coeff, alpha, lef_state_1, lef_state_2, nlplant):
        
        nlplant.atmos(ctypes.c_double(h),ctypes.c_double(V),ctypes.c_void_p(coeff.ctypes.data))
        atmos_out = coeff[1]/coeff[2] * 9.05
        alpha_deg = alpha*180/np.pi
        
        LF_err = alpha_deg - (lef_state_1 + (2 * alpha_deg))
        #lef_state_1 += LF_err*7.25*time_step
        LF_out = (lef_state_1 + (2 * alpha_deg)) * 1.38
        
        lef_cmd = LF_out + 1.45 - atmos_out
        
        # command saturation
        lef_cmd = np.clip(lef_cmd,act_lim[1][4],act_lim[0][4])
        # rate saturation
        lef_err = np.clip((1/0.136) * (lef_cmd - lef_state_2),-25,25)
        # integrate
        #lef_state_2 += lef_err*time_step
        
        return LF_err*7.25, lef_err
    
    def calc_xdot(self, x, u):
        
        # initialise variables
        xdot = np.zeros(18)
        temp = np.zeros(6)
        coeff = np.zeros(3)
        
        #--------------Thrust Model--------------#
        temp[0] = self.upd_thrust(u[0], x[12])
        #--------------Dstab Model---------------#
        temp[1] = self.upd_dstab(u[1], x[13])
        #-------------aileron model--------------#
        temp[2] = self.upd_ail(u[2], x[14])
        #--------------rudder model--------------#
        temp[3] = self.upd_rud(u[3], x[15])
        #--------leading edge flap model---------#
        temp[5], temp[4] = self.upd_lef(x[2], x[6], coeff, x[7], x[17], x[16], self.nlplant)
        
        #----------run nlplant for xdot----------#
        self.nlplant.Nlplant(ctypes.c_void_p(x.ctypes.data), ctypes.c_void_p(xdot.ctypes.data), ctypes.c_int(self.fi_flag))    
        
        xdot[12:18] = temp
        
        return xdot
        
    def step(self, action):
        
        self.u = action
        
        # find xdot
        self.xdot = self.calc_xdot(self.x, self.u)
        
        # update x
        self.x += self.xdot*self.dt
        self.t += self.dt
        
        reward = 1
        isdone = False
        info = {'fidelity':'high'}
        
        return self.get_obs(), reward, isdone, info
    
    def get_obs(self):
        
        return  np.array(list(self.x[i] for i in range(17)), dtype='float32')
    
    def reset(self):
        
        self.x = x0
        self.u = x0[12:16]
        self.t = 0
        
        return self.get_obs()
    
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
    
        eps = 1e-05
        
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
    
    
sim = F16_env(x0,paras_sim)

A,B,C,D = sim.linearise()

check_env(sim, warn=True)

# env = DummyVecEnv([lambda: F16_env(x0, paras_sim)])

model = A2C('MlpPolicy', sim, verbose = 1)
# model.learn(total_timesteps=10000)