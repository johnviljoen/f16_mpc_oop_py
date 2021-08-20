#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug  1 20:09:42 2021

@author: johnviljoen
"""

import unittest
import progressbar
import numpy as np

from env import F16
from parameters import state_vector, input_vector, simulation_parameters, nlplant
from utils import vis_x, vis_u, dlqr, vis_mpc_u, vis_mpc_x
from scipy.signal import cont2discrete
import matplotlib.pyplot as plt


# from control.matlab import *
from parameters import x_lb, x_ub, u_lb, u_ub, udot_lb, udot_ub

class test_F16(unittest.TestCase, F16):
    
    def __init__(self, state_vector, input_vector, simulation_parameters, nlplant):
        super(unittest.TestCase, self).__init__(state_vector, input_vector, simulation_parameters, nlplant)
    
    def test_upd_thrust(self):
        # self.assertAlmostEqual()
        """ Function to test that the thrust model is behaviing as expected. This
        is done by creating a seperate control system modelled in a different way
        and then verifying the outputs are the same
        
        The model is first order and so it is also simple to test the command
        and rate saturations to ensure accuracy"""
        pass
        
        
        
    def test_act_cmd_lims(self):
        
        """ Function to test that the actuator command limits are implemented correctly by
        testing the boundaries of their operating envelope. It should be noted that
        this test is conducted simulatenously on all actuators, and therefore 
        it assumes that the actuator dynamics are not coupled, which they shouldnt
        be... it is fairly simple to see they arent in the code written in "env.py". 
                
        The engine is tested by commanding it to go 1000 below and 1000 above
        its maximum value whilst having the current engine state be 1 below and 
        1 above its maximum value respectively. This test is successful if  and 
        only if this command is seen to be ignored by the generation of an engine 
        state time derivative in the opposite direction. 
        
        The other actuators will crash the C code if they go above or below their 
        maximum and minimum values respectively. Therefore their states are set
        at their maximums and minimums and they are commanded to go beyond them.
        The test is successful if and only if their rate of change is found to be
        zero"""
        
        self.reset()
        
        # set command 1000 below the minimum
        self.u.values[0] = u_lb[0] - 1000
        self.u.values[1] = u_lb[1] - 1000
        self.u.values[2] = u_lb[2] - 1000
        self.u.values[3] = u_lb[3] - 1000
        
        # set the current engine state 1 below the minimum
        self.x.values[12] = u_lb[0] - 1
        
        # the C code crashes if ordered to lookup data beyond its tables for the
        # other actuators, therefore dh, da, dr are set to their minimum exactly
        # and it is tested if a command can decrease this
        self.x.values[13] = u_lb[1] 
        self.x.values[14] = u_lb[2]
        self.x.values[15] = u_lb[3]
        
        # calculate the state vector time derivative
        xdot = self._calc_xdot(self.x.values, self.u.values)
        
        # check that the rate of change of the engine state is indeed positive
        self.assertGreater(xdot[12], 0)
        
        # check the rate of change of dh, da, dr is indeed zero
        self.assertAlmostEqual(xdot[13], 0)
        self.assertAlmostEqual(xdot[14], 0)
        self.assertAlmostEqual(xdot[15], 0)

        # check the maximums now, but first reset the simulation.
        self.reset()
        
        # set command 1000 above the maximum
        self.u.values[0] = u_ub[0] + 1000
        self.u.values[1] = u_ub[1] + 1000
        self.u.values[2] = u_ub[2] + 1000
        self.u.values[3] = u_ub[3] + 1000
        
        # set the current engine state 1 above the maximum
        self.x.values[12] = u_ub[0] + 1
        
        # set dh, da, dr to their maximums
        self.x.values[13] = u_ub[1] 
        self.x.values[14] = u_ub[2]
        self.x.values[15] = u_ub[3]
        
        # calculate the state vector time derivative
        xdot = self._calc_xdot(self.x.values, self.u.values)
        
        # check that the rate of change of the engine state is indeed negative
        self.assertLess(xdot[12], 0)
        
        # check that the rate of change of dh, da, dr is indeed zero
        self.assertAlmostEqual(xdot[13], 0)
        self.assertAlmostEqual(xdot[14], 0)
        self.assertAlmostEqual(xdot[15], 0)
        
    def test_act_rate_lims(self):
        
        """ Function to test the rate limits of the 1st order actuators is behaving as 
        expected. Note this does not include the engine as it is a more complex system."""
        
        # begin from rough trim
        self.reset()
        
        # command maximums on all actuators
        self.u.values[1] = u_ub[1]
        self.u.values[2] = u_ub[2]
        self.u.values[3] = u_ub[3]
        
        xdot = self._calc_xdot(self.x.values, self.u.values)
        
        self.assertAlmostEqual(xdot[13], 60)
        self.assertAlmostEqual(xdot[14], 80)
        self.assertAlmostEqual(xdot[15], 120)
        
        self.reset()
                        
        # now test the inverse
        self.u.values[1] = u_lb[1]
        self.u.values[2] = u_lb[2]
        self.u.values[3] = u_lb[3]
                
        xdot = self._calc_xdot(self.x.values, self.u.values)

        self.assertAlmostEqual(xdot[13], -60)
        self.assertAlmostEqual(xdot[14], -80)
        self.assertAlmostEqual(xdot[15], -120)
        
    def test_aerodynamics(self):
        
        pass
    
    def test_LQR_lin(self, f16=True):
        
        if f16:
        
            # trim the simulation and set it as the initial condition
            self.x.initial_condition, _ = self.trim(10000,700)
            self.u.initial_condition = self.x.initial_condition[12:16]
            self.reset()
            
            A,B,C,D = self.linearise(self.x._get_mpc_x(), self.u._get_mpc_u(), _calc_xdot=self._calc_xdot_na, get_obs=self._get_obs_na)
            A,B,C,D = cont2discrete((A,B,C,D), self.paras.dt)[0:4]    
            
            x = self.x._get_mpc_x()[:,None]
            u = self.u._get_mpc_u()[:,None]
        
            # reference x is trim, i.e. the current state, therefore it should stay there
            x_ref = np.copy(x)
        
        else:
            
            # example
            
            self.paras.dt = 0.1
            self.paras.time_end = 3
            
            x = np.array([3, 1])[np.newaxis].T
            u = np.array([0])
            x_ref = np.array([-3, 0])[np.newaxis].T
            
            A = np.array([[1, 1.0], [0, 1]])
            B = np.array([0.0, 1])[:,None]
            Q = np.array([[1.0, 0.0], [0.0, 0.0]])
            R = np.array([[1.0]])
        
        Q = np.eye(len(x))
        R = np.eye(len(u))*0.01
        
        K = dlqr(A,B,Q,R)
                
        rng = np.linspace(self.paras.time_start, self.paras.time_end, int((self.paras.time_end-self.paras.time_start)/self.paras.dt))
        # bar = progressbar.ProgressBar(maxval=len(rng)).start()
        
        # create storage
        x_storage = np.zeros([len(rng),len(x)])
        u_storage = np.zeros([len(rng),len(u)])
        
        for i, val in enumerate(rng):
            
            u = - K @ (x - x_ref)
            x = A @ x + B @ u
                        
            x_storage[i,:] = x[:,0]
            u_storage[i,:] = u[:,0]
            
        if f16:
            
            vis_mpc_u(u_storage, rng)
            vis_mpc_x(x_storage, rng)
            
        else:
            
            plt.plot(rng, u_storage, "-r", label="input")
            plt.plot(rng, x_storage[:,0], "-b", label="x1")
            plt.plot(rng, x_storage[:,1], "-g", label="x2")
            plt.grid(True)
            plt.xlim([0, self.paras.time_end])
            plt.title("LQR Regulator")
            plt.legend()
            plt.show()
            
        return x_storage, u_storage, rng
        
        
    def test_LQR(self):
        
        self.x.initial_condition, _ = self.trim(10000,700)
        self.u.initial_condition = self.x.initial_condition[12:16]
        self.reset()
        
        print('x:',self.x._get_mpc_x())
        print('u:',self.u._get_mpc_u())
        
        # now we have the simulation p r i m e d for timehistory from trim
        # if the LQR deviates from where it is right now theres a problem
        # assuming LQR is implemented as u = -K @ (xref - x)
        
        rng = np.linspace(self.paras.time_start, self.paras.time_end, int((self.paras.time_end-self.paras.time_start)/self.paras.dt))
        # bar = progressbar.ProgressBar(maxval=len(rng)).start()
        
        # create storage
        x_storage = np.zeros([len(rng),len(self.x.values)])
        u_storage = np.zeros([len(rng),len(self.u._get_mpc_u())])
        xdot_storage = np.zeros([len(rng),len(self.x.values)])
        
        x_ref = np.copy(self.x._get_mpc_x())
        
        p_dem = 0
        q_dem = 0
        r_dem = 0
        
        x_ref[5] = p_dem
        x_ref[6] = q_dem
        x_ref[7] = r_dem
        
        K = self._calc_LQR_gain()
        print('K:',K)
        print('x:',self.x._get_mpc_x())
        
        # return K, self.x._get_mpc_x(), x_ref
        
        # exit()
        
        for idx, val in enumerate(rng):
            
            # print('idx:', idx)
            
            
            cmd = (- K @ (self.x._get_mpc_x() - x_ref)) * np.pi/180
            u_storage[idx,:] = cmd
            self.u.values[1:] = cmd
            # self.u.values[0] = self.u.initial_condition[0]
            # print('u:',self.u.values)
            # self.u.values = np.copy(self.u.initial_condition)
            
            print('u:',self.u.values)
            # print('x:',self.x.values)
            
            self.step(self.u.values)
            
            x_storage[idx,:] = self.x.values
            # bar.update(idx)
            
        vis_x(x_storage, rng)
        vis_u(u_storage, rng)
        
    
    def test_MPC(self):
        
        """ Function to simulate the MPC controlled F16 to test it is behaving correctly
        
        exact methods are TBD """
        
        rng = np.linspace(self.paras.time_start, self.paras.time_end, int((self.paras.time_end-self.paras.time_start)/self.paras.dt))
        # bar = progressbar.ProgressBar(maxval=len(rng)).start()
        
        # create storage
        x_storage = np.zeros([len(rng),len(self.x.values)])
        u_storage = np.zeros([len(rng),len(self.u._get_mpc_u())])
        xdot_storage = np.zeros([len(rng),len(self.x.values)])
        
        for idx, val in enumerate(rng):
            
            p_dem = 0 * np.pi/180 # rad
            q_dem = 0 * np.pi/180  # rad
            r_dem = 0   # rad
            
            print('idx:', idx)
            
            cmd = self._calc_MPC_action(p_dem, q_dem, r_dem,10)
            u_storage[idx,:] = cmd
            self.u.values[1:] = cmd
            # self.u.values[0] = self.u.initial_condition[0]
            print('u:',self.u.values)
            # self.u.values = np.copy(self.u.initial_condition)
            
            # print('u:',self.u.values)
            # print('x:',self.x.values)
            
            self.step(self.u.values)
            
            x_storage[idx,:] = self.x.values
            #bar.update(idx)
            
        vis_x(x_storage, rng)
        vis_u(u_storage, rng)
        

test_f16 = test_F16(state_vector, input_vector, simulation_parameters, nlplant)

x_storage, u_storage, rng = test_f16.test_LQR_lin(f16=True)
# x = test_f16.test_LQR_lin()