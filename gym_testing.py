#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 10 13:54:33 2021

@author: johnviljoen
"""

from parameters import x_lim, act_lim
from parameters import initial_state_vector_ft_rad as x0
from parameters import simulation_parameters as paras_sim
from parameters import paras_mpc

from ctypes import CDLL
import ctypes
import os

import numpy as np

import gym
from gym import spaces

import progressbar
from utils import tic, toc, vis

# import scipy fmin for trim function
from scipy.optimize import minimize

class F16_env(gym.Env):
    """Custom Environment that follows gym interface"""
  
    metadata = {'render.modes': ['human']}

    def __init__(self, act_lim, x_lim, x0, paras_sim, paras_mpc):
        super(F16_env, self).__init__()
        # Define action and observation space
        # They must be gym.spaces objects
        
        # ignoring lef actuator limits as they are not directly commanded
        self.action_space = spaces.Box(low=np.array(act_lim[1])[0:4], high=np.array(act_lim[0])[0:4])
        
        self.observation_space = spaces.Box(low=np.array(x_lim[1]), high=np.array(x_lim[0]))
        
        self.ic = x0
        
        self.x = x0
        
        self.u = x0[12:16]
        
        self.dt = paras_sim[0]
        
        self.paras_sim = paras_sim
        
        self.paras_mpc = paras_mpc
        
        self.coeff = np.zeros(3)
        
        # create interface with c shared library .so file in folder "C"
        if paras_sim[3] == 1:
            so_file = os.getcwd() + "/C/nlplant_xcg35.so"
        elif paras_sim[3] == 0:
            so_file = os.getcwd() + "/C/nlplant_xcg25.so"
            
        nlplant = CDLL(so_file)
        
        self._nlplant = nlplant
        
        self.xdot = np.zeros(len(x0))
        
        
    def _upd_thrust(self):
        # command saturation
        self.u[0] = np.clip(self.u[0],self.action_space.low[0],self.action_space.high[0])
        # rate saturation
        T_grad = np.clip(self.u[0] - self.x[12], -10000, 10000)
        return T_grad
    
    def _upd_dstab(self):
        # command saturation
        self.u[1] = np.clip(self.u[1],self.action_space.low[1],self.action_space.high[1])
        # rate saturation
        dh_grad = np.clip(20.2*(self.u[1] - self.x[13]), -60, 60)
        return dh_grad
    
    def _upd_ail(self):
        # command saturation
        self.u[2] = np.clip(self.u[2],self.action_space.low[2],self.action_space.high[2])
        # rate saturation
        da_grad = np.clip(20.2*(self.u[2] - self.x[14]), -80, 80)
        return da_grad
    
    def _upd_rud(self):
        # command saturation
        self.u[3] = np.clip(self.u[3],self.action_space.low[3],self.action_space.high[3])
        # rate saturation
        dr_grad = np.clip(20.2*(self.u[3] - self.x[15]), -120, 120)
        return dr_grad
    
    def _upd_lef(self):
        
        self._nlplant.atmos(ctypes.c_double(self.x[2]),ctypes.c_double(self.x[6]),ctypes.c_void_p(self.coeff.ctypes.data))
        atmos_out = self.coeff[1]/self.coeff[2] * 9.05
        alpha_deg = self.x[7]*180/np.pi
        
        LF_err = alpha_deg - (self.x[17] + (2 * alpha_deg))
        #lef_state_1 += LF_err*7.25*time_step
        LF_out = (self.x[17] + (2 * alpha_deg)) * 1.38
        
        lef_cmd = LF_out + 1.45 - atmos_out
        
        # command saturation
        lef_cmd = np.clip(lef_cmd,0,25)
        # rate saturation
        lef_err = np.clip((1/0.136) * (lef_cmd - self.x[16]),-25,25)        
        return LF_err*7.25, lef_err
    
    def calc_xdot(self):
        
        # initialise variables
        #xdot = np.zeros(18)
        temp = np.zeros(6)
        
        #--------------Thrust Model--------------#
        temp[0] = self._upd_thrust()
        #--------------Dstab Model---------------#
        temp[1] = self._upd_dstab()
        #-------------aileron model--------------#
        temp[2] = self._upd_ail()
        #--------------rudder model--------------#
        temp[3] = self._upd_rud()
        #--------leading edge flap model---------#
        temp[5], temp[4] = self._upd_lef()
        
        #----------run nlplant for xdot----------#
        self._nlplant.Nlplant(ctypes.c_void_p(self.x.ctypes.data), ctypes.c_void_p(self.xdot.ctypes.data), ctypes.c_int(self.paras_sim[4]))    
        
        self.xdot[12:18] = temp
        
        #return xdot
        
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
        if self.x[7] > x_lim[0][7]*np.pi/180:
            self.x[7] = x_lim[0][7]*np.pi/180
        elif self.x[7] < x_lim[1][7]*np.pi/180:
            self.x[7] = x_lim[1][7]*np.pi/180
            
        self.u = self.x[12:16]
        
        # calc xdot
        self.calc_xdot()
        
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
    
    # def linearise(self, x, u, output_vars, fi_flag, nlplant):
    
    #     eps = 1e-06
        
    #     A = np.zeros([len(self.x),len(self.x)])
    #     B = np.zeros([len(self.x),len(self.u)])
    #     C = np.zeros([len(output_vars),len(x)])
    #     D = np.zeros([len(output_vars),len(u)])
        
    #     # Perturb each of the state variables and compute linearization
    #     for i in range(len(x)):
            
    #         dx = np.zeros((len(x),))
    #         dx[i] = eps
            
    #         xdot_control
            
    #         A[:, i] = (calc_xdot(x + dx, u, fi_flag, nlplant) - calc_xdot(x, u, fi_flag, nlplant)) / eps
    #         C[:, i] = (calc_out(x + dx, u, output_vars) - calc_out(x, u, output_vars)) / eps
            
    #     # Perturb each of the input variables and compute linearization
    #     for i in range(len(u)):
            
    #         du = np.zeros((len(u),))
    #         du[i] = eps
                    
    #         B[:, i] = (calc_xdot(x, u + du, fi_flag, nlplant) - calc_xdot(x, u, fi_flag, nlplant)) / eps
    #         D[:, i] = (calc_out(x, u + du, output_vars) - calc_out(x, u, output_vars)) / eps
        
    #     return A, B, C, D
        
    
    def step(self, action):
        
        # apply action
        self.u = action
        
        # find xdot
        self.calc_xdot()
        
        # update x
        self.x += self.xdot*self.dt
        
        return self.x
    
    def reset(self):
        
        self.x = self.ic
        # Randomise inputs
        
        
    def render(self, mode='human', close=False):
        pass
    


model = F16_env(act_lim, x_lim, x0, paras_sim, paras_mpc)
model.trim(10000,700)

rng = np.linspace(paras_sim[1], paras_sim[2], int((paras_sim[2]-paras_sim[1])/paras_sim[0]))
output_vars = [6,7,8,9,10,11]

# create storage
x_storage = np.zeros([len(rng),len(model.x)])
A = np.zeros([len(model.x),len(model.x),len(rng)])
B = np.zeros([len(model.x),len(model.u),len(rng)])
C = np.zeros([len(output_vars),len(model.x),len(rng)])
D = np.zeros([len(output_vars),len(model.u),len(rng)])

bar = progressbar.ProgressBar(maxval=len(rng)).start()

tic()

for idx, val in enumerate(rng):
    
    #----------------------------------------#
    #------------linearise model-------------#
    #----------------------------------------#
    
    #[A[:,:,idx], B[:,:,idx], C[:,:,idx], D[:,:,idx]] = linearise(x, u, output_vars, fi_flag, nlplant)
    
    #----------------------------------------#
    #--------------Take Action---------------#
    #----------------------------------------#
    
    # MPC prediction using squiggly C and M matrices
    #CC, MM = calc_MC(paras_mpc[0], A[:,:,idx], B[:,:,idx], time_step)
    
    
    #----------------------------------------#
    #--------------Integrator----------------#
    #----------------------------------------#    
    
    x = model.step(model.u)
    
    #----------------------------------------#
    #------------Store History---------------#
    #----------------------------------------#
    
    x_storage[idx,:] = x
    
    bar.update(idx)
    
toc()
    
vis(x_storage, rng)
