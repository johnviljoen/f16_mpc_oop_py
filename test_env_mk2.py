#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 25 16:55:15 2021

@author: johnviljoen
"""

import unittest
import progressbar
import numpy as np

from env import F16
from utils import *
from scipy.signal import cont2discrete
import matplotlib.pyplot as plt

from parameters import x_lb, x_ub, u_lb, u_ub, udot_lb, udot_ub

class test_F16(unittest.TestCase, F16):
    
    def __init__(self, state_vector, input_vector, simulation_parameters, state_space, nlplant):
        super(unittest.TestCase, self).__init__(state_vector, input_vector, simulation_parameters, state_space, nlplant)
    
    def LQR(self, linear=False):
        
        # setup LQR controlled simulation
        self.paras.time_end = 10
        
        rng = np.linspace(self.paras.time_start, self.paras.time_end, int((self.paras.time_end-self.paras.time_start)/self.paras.dt))
        
        # create storage
        x_storage = np.zeros([len(rng),self.ss.Ad.shape[0]])
        u_storage = np.zeros([len(rng),self.ss.Bd.shape[1]])
        
        # find the offline gain for cruise
        K = self._calc_LQR_gain()
        
        # pick a reference point, 0's for no movement
        p_dem = 0
        q_dem = 0
        r_dem = 0
        
        if linear: # simulating on the linearised state space system
        
            x = self.x._get_mpc_x()
            u0 = self.u.initial_condition[1:]
            u = np.copy(u0)
            
            # create storage
            x_storage = np.zeros([len(rng),len(x)])
            u_storage = np.zeros([len(rng),len(u)])
            
            for idx, val in enumerate(rng):
                
                print('idx:', idx)
                
                u = self._calc_LQR_action(p_dem, q_dem, r_dem, K, x, u0)
                
                x = self.ssr.Ad @ x + self.ssr.Bd @ u
                
                x_storage[idx,:] = x
                u_storage[idx,:] = u
                
            vis_mpc_u(u_storage, rng)
            vis_mpc_x(x_storage, rng)
        
        else: # as you might expect... simulating on the nonlinear system
            
            for idx, val in enumerate(rng):
            
                print('idx:', idx)
                
                # use the built in function for determining LQR action for consistency
                u = self._calc_LQR_action(p_dem, q_dem, r_dem, K, self.x._get_mpc_x(), self.u.initial_condition[1:])
                
                # of course we arent controlling T here so only change final 3 inputs
                self.u.values[1:] = u                
                
                print('u:',self.u.values)
                
                # use the u values calculated to use the built in step function
                self.step(self.u.values)
                
                x_storage[idx,:] = self.x.values
                
            vis_x(x_storage, rng)
            vis_u(u_storage, rng)
        
        