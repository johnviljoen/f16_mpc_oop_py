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
# from parameters import state_vector, input_vector, simulation_parameters, state_space, nlplant
from utils import *
from scipy.signal import cont2discrete
import matplotlib.pyplot as plt
from control.matlab import *

# from control.matlab import *
from parameters import x_lb, x_ub, u_lb, u_ub, udot_lb, udot_ub

class test_F16(unittest.TestCase, F16):
    
    def __init__(self, state_vector, input_vector, simulation_parameters, state_space, nlplant):
        super(unittest.TestCase, self).__init__(state_vector, input_vector, simulation_parameters, state_space, nlplant)
    
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
    
    def test_linearisation(self):
        
        """ In this function we import both the Python and MATLAB linearisations
        conducted at 10000ft and 700fps wings level flight in order to compare the two"""
        
        # MATLAB output for steady wings level cruise at 10000ft 700fps
        # MATLAB_evals = np.array([
        # 0.0000 + 0.0000j,
        # 0.0000 + 0.0000j,
        # 0.0000 + 0.0000j,
        # -7.3529 + 0.0000j,
        # -7.2500 + 0.0000j,
        # -1.3929 + 2.7668j,
        # -1.3929 - 2.7668j,
        # -0.0067 + 0.0670j,
        # -0.0067 - 0.0670j,
        # 0.0000 + 0.0000j,
        # -0.4478 + 3.9347j,
        # -0.4478 - 3.9347j,
        # -3.7888 + 0.0000j,
        # -0.0089 + 0.0000j,
        # -1.0000 + 0.0000j,
        # -20.2000 + 0.0000j,
        # -20.2000 + 0.0000j, 
        # -20.2000 + 0.0000j])
        
        py_A,py_B,py_C,py_D = self.linearise(self.x.values, self.u.values)
        
        # Python_evals = np.linalg.eig(A)
        
        # print("MATLAB_evals:", MATLAB_evals)
        # print("Python_evals:", Python_evals)
        
        mat = scipy.io.loadmat("MATLAB_SS.mat")
        
        MAT_A = mat['A']
        MAT_B = mat['B']
        MAT_C = mat['C']
        MAT_D = mat['D']
        
        # lets simulate a timehistory of both linear systems        
        self.paras.time_end = 10
        
        rng = np.linspace(self.paras.time_start, self.paras.time_end, int((self.paras.time_end-self.paras.time_start)/self.paras.dt))
        
        # create storage
        MAT_x_storage = np.zeros([len(rng),self.ss.Ad.shape[0]])
        py_x_storage = np.zeros([len(rng),self.ss.Ad.shape[0]])
        u_storage = np.zeros([len(rng),self.ss.Bd.shape[1]])
        
        py_x = np.copy(self.x.initial_condition)
        MAT_x = np.copy(self.x.initial_condition)
        u = np.copy(self.u.initial_condition)
        
        for idx, val in enumerate(rng):
            
            print('idx:', idx)
            
            py_x = py_A @ py_x + py_B @ u
            MAT_x = MAT_A @ MAT_x + MAT_B @ u
            
            py_x_storage[idx,:] = py_x
            MAT_x_storage[idx,:] = MAT_x
            u_storage[idx,:] = u
            
        vis_x(py_x_storage, rng)
        vis_x(MAT_x_storage, rng)
        
        print('MATLAB eigenvalues: \n', np.linalg.eig(MAT_A)[0])
        print('python eigenvalues: \n', np.linalg.eig(py_A)[0])
        
    def offline_full_LQR_lin(self):
        
        self.paras.time_end = 10
        
        rng = np.linspace(self.paras.time_start, self.paras.time_end, int((self.paras.time_end-self.paras.time_start)/self.paras.dt))
        
        # create storage
        x_storage = np.zeros([len(rng),self.ss.Ad.shape[0]])
        u_storage = np.zeros([len(rng),self.ss.Bd.shape[1]])
        
        Q = self.ss.Cd.T @ self.ss.Cd
        R = np.eye(4)
        
        K = dlqr(self.ss.Ad, self.ss.Bd, Q, R)
        
        x_ref = np.copy(self.x.values)
        
        u0 = np.copy(self.u.values)
        x = np.copy(self.x.values)
        u = np.copy(u0)
        
        for idx, val in enumerate(rng):
            
            print('idx:', idx)
            print('u:',u)
            
            u = u0 - K @ (x - x_ref)
            u = u0
            
            x = self.ss.Ad @ x + self.ss.Bd @ u
            
            x_storage[idx,:] = x
            u_storage[idx,:] = u
                        
        vis_mpc_x(x_storage, rng)
        vis_mpc_u(u_storage, rng)
        
    
    def offline_LQR_nl(self):
        
        self.paras.time_end = 10
        # self.paras.dt = 1/60
        
        rng = np.linspace(self.paras.time_start, self.paras.time_end, int((self.paras.time_end-self.paras.time_start)/self.paras.dt))
        
        # create storage
        x_storage = np.zeros([len(rng),self.ss.Ad.shape[0]])
        u_storage = np.zeros([len(rng),self.ss.Bd.shape[1]])
        
        Q = self.ssr.Cd.T @ self.ssr.Cd
        R = np.eye(3)
        
        K = dlqr(self.ssr.Ad, self.ssr.Bd, Q, R)
        
        x_ref = np.copy(self.x._get_mpc_x())
        x_ref[4] = 0
        x_ref[5] = 0
        x_ref[6] = 0
        
        u0 = self.u.initial_condition[1:]
        
        for idx, val in enumerate(rng):
            
            print('idx:', idx)

            self.step(self.u.values)
            print('u:',self.u.values)
            x_storage[idx,:] = self.x.values
            u_storage[idx,:] = self.u.values
            
            u = u0 - K @ (self.x._get_mpc_x() - x_ref)
            u = self._calc_LQR_action(0,0,0,K)
            self.u.values[1:] = u
            
        vis_x(x_storage, rng)
        vis_u(u_storage, rng)
        
    def SSR_continuous_PID_lin(self):
        
        self.paras.time_end = 10
        
        rng = np.linspace(self.paras.time_start, self.paras.time_end, int((self.paras.time_end-self.paras.time_start)/self.paras.dt))
        
        # create storage
        x_storage = np.zeros([len(rng),self.ssr.Ac.shape[0]])
        u_storage = np.zeros([len(rng),self.ssr.Bc.shape[1]])
        
        q_dem = 0
        
        u0 = np.copy(self.u._get_mpc_u())
        x = np.copy(self.x._get_mpc_x())
        u = np.copy(u0)
        
        error_int = 0
        
        # with negatives here it is trying to self correct! 
        # Stability has so far failed with gains attempted
        Kp = -10
        Ki = -5
        
        for idx, val in enumerate(rng):
            
            # print('idx:', idx)
            # print('u:',u)
            
            error = q_dem-x[5]
            error_int += error
            
            dh_cmd = error * Kp + error_int * Ki + u0[0]
            
            da_cmd = 0
            dr_cmd = 0
            
            u = np.array([dh_cmd, da_cmd, dr_cmd])
            print('error:',error)
            print('error_int:',error_int)
            print('u:',u)
            
            xdot = self.ssr.Ac @ x + self.ssr.Bc @ u
            x += xdot*self.paras.dt
            
            # print(np.max(np.real(np.linalg.eig(self.ssr.Ac - self.ssr.Bc@K)[0])))
            
            if idx == 1000:
                exit()
            
            x_storage[idx,:] = x
            u_storage[idx,:] = u
                        
        vis_mpc_x(x_storage, rng)
        vis_mpc_u(u_storage, rng)
        
    def SSR_continuous_LQR_lin(self):
        
        self.paras.time_end = 10
        
        rng = np.linspace(self.paras.time_start, self.paras.time_end, int((self.paras.time_end-self.paras.time_start)/self.paras.dt))
        
        # create storage
        x_storage = np.zeros([len(rng),self.ssr.Ac.shape[0]])
        u_storage = np.zeros([len(rng),self.ssr.Bc.shape[1]])
        
        Q = self.ssr.Cc.T @ self.ssr.Cc
        R = np.eye(3)
        
        K = np.array(lqr(self.ssr.Ac, self.ssr.Bc, Q, R)[0])
        print(K)
        
        x_ref = np.copy(self.x._get_mpc_x())
        
        u0 = np.copy(self.u._get_mpc_u())
        x = np.copy(self.x._get_mpc_x())
        u = np.copy(u0)
        
        for idx, val in enumerate(rng):
            
            # print('idx:', idx)
            # print('u:',u)
            
            u = (u0 - K @ (x - x_ref))
            # u = u0
            print('u:',u)
            
            xdot = self.ssr.Ac @ x + self.ssr.Bc @ u
            x += xdot*self.paras.dt
            
            # print(np.max(np.real(np.linalg.eig(self.ssr.Ac - self.ssr.Bc@K)[0])))
            
            x_storage[idx,:] = x
            u_storage[idx,:] = u
                        
        vis_mpc_x(x_storage, rng)
        vis_mpc_u(u_storage, rng)
        
    def SSR_discrete_LQR_lin(self):
        
        self.paras.time_end = 10
        
        rng = np.linspace(self.paras.time_start, self.paras.time_end, int((self.paras.time_end-self.paras.time_start)/self.paras.dt))
        
        # create storage
        x_storage = np.zeros([len(rng),self.ssr.Ad.shape[0]])
        u_storage = np.zeros([len(rng),self.ssr.Bd.shape[1]])
        
        Q = self.ssr.Cd.T @ self.ssr.Cd
        Q[4,4] = 100
        Q[5,5] = 100
        Q[6,6] = 100
        R = np.eye(3)
        R[0,0] = 0
        R[1,1] = 0
        R[2,2] = 0
        
        K = dlqr(self.ssr.Ad, self.ssr.Bd, Q, R)
        
        x_ref = np.copy(self.x._get_mpc_x())
        
        u0 = np.copy(self.u._get_mpc_u())
        x = np.copy(self.x._get_mpc_x())
        u = np.copy(u0)
        
        for idx, val in enumerate(rng):
            
            # print('idx:', idx)
            # print('u:',u)
            
            u = u0 - K @ (x - x_ref)
            # u = u0
            print(x-x_ref)
            
            x = self.ssr.Ad @ x + self.ssr.Bd @ u
            
            x_storage[idx,:] = x
            u_storage[idx,:] = u
                        
        vis_mpc_x(x_storage, rng)
        vis_mpc_u(u_storage, rng)
        
        return x, x_ref, K, self.ssr.Ad, self.ssr.Bd, u0, u
        
    def test_control(self):
        
        """ Function to simulate the MPC controlled F16 to test it is behaving correctly
        
        exact methods are TBD """   
        
        rng = np.linspace(self.paras.time_start, self.paras.time_end, int((self.paras.time_end-self.paras.time_start)/self.paras.dt))
        
        # create storage
        x_storage = np.zeros([len(rng),len(self.x.values)])
        u_storage = np.zeros([len(rng),len(self.u._get_mpc_u())])
        
        for idx, val in enumerate(rng):
            
            print('idx:', idx)
            print('u:',self.u.values)

            self.step(self.u.values)
            x_storage[idx,:] = self.x.values
            
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
        
        for idx, val in enumerate(rng):
            
            p_dem = 0 * np.pi/180 # rad
            q_dem = 0 * np.pi/180  # rad
            r_dem = 0   # rad
            
            print('idx:', idx)
            
            cmd = self._calc_MPC_action(p_dem, q_dem, r_dem,10)
            u_storage[idx,:] = cmd
            self.u.values[1:] = cmd
            print('u:',self.u.values)
            
            self.step(self.u.values)
            
            x_storage[idx,:] = self.x.values
            #bar.update(idx)
            
        vis_x(x_storage, rng)
        vis_u(u_storage, rng)
        
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
            
            Q = np.eye(len(x))
        
            R = np.eye(len(u))
        
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
            
        return x_storage, u_storage, rng, K
        
        
    def test_LQR_static_nl(self):
        
        # now we have the simulation p r i m e d for timehistory from trim
        # if the LQR deviates from where it is right now theres a problem
        # assuming LQR is implemented as u = -K @ (xref - x)
        p_dem = 0
        q_dem = 0
        r_dem = 0
        
        x = self.x._get_mpc_x()
        u = self.u._get_mpc_u()
        
        K = self._calc_LQR_gain()
                
        rng = np.linspace(self.paras.time_start, self.paras.time_end, int((self.paras.time_end-self.paras.time_start)/self.paras.dt))
        # bar = progressbar.ProgressBar(maxval=len(rng)).start()
        
        # create storage
        x_storage = np.zeros([len(rng),len(self.x.values)])
        u_storage = np.zeros([len(rng),len(self.u._get_mpc_u())])
        
        x_ref = np.copy(self.x._get_mpc_x())
        x_ref[4] = p_dem
        x_ref[5] = q_dem
        x_ref[6] = r_dem
        
        u0 = np.copy(u)

        for idx, val in enumerate(rng):
            
            print('idx:', idx)
            
            
            u = self._calc_LQR_action(p_dem, q_dem, r_dem, K)
            
            self.u.values[1:] = u                
            
            print('u:',self.u.values)
            
            self.step(self.u.values)
            
            x_storage[idx,:] = self.x.values
            # bar.update(idx)
            
        vis_x(x_storage, rng)
        vis_u(u_storage, rng)
        
    def test_LQR_dynamic_nl(self):
        
        self.paras.time_end = 2
        
        # now we have the simulation p r i m e d for timehistory from trim
        # if the LQR deviates from where it is right now theres a problem
        # assuming LQR is implemented as u = -K @ (xref - x)
        
        A,B,C,D = self.linearise(self.x._get_mpc_x(), self.u._get_mpc_u(), _calc_xdot=self._calc_xdot_na, get_obs=self._get_obs_na)
        A,B,C,D = cont2discrete((A,B,C,D), self.paras.dt)[0:4]    
        
        x = self.x._get_mpc_x()[:,None]
        u = self.u._get_mpc_u()[:,None]
        
        Q = np.eye(len(x))
        R = np.eye(len(u))*10000
        
        # reference x is trim, i.e. the current state, therefore it should stay there
        x_ref = np.copy(x)
        
        rng = np.linspace(self.paras.time_start, self.paras.time_end, int((self.paras.time_end-self.paras.time_start)/self.paras.dt))
        # bar = progressbar.ProgressBar(maxval=len(rng)).start()
        
        # create storage
        x_storage = np.zeros([len(rng),len(self.x.values)])
        u_storage = np.zeros([len(rng),len(self.u._get_mpc_u())])
        xdot_storage = np.zeros([len(rng),len(self.x.values)])
        
        x_ref = np.copy(self.x._get_mpc_x())

        for idx, val in enumerate(rng):
            
            print('idx:', idx)
            
            Ac,Bc,Cc,Dc = self.linearise(self.x._get_mpc_x(), self.u._get_mpc_u(), _calc_xdot=self._calc_xdot_na, get_obs=self._get_obs_na)
            A,B,C,D = cont2discrete((Ac,Bc,Cc,Dc), self.paras.dt)[0:4]   
            K = dlqr(A,B,Q,R)
            
            evals, _ = np.linalg.eig(A-B@K)
            
            if np.amax(np.abs(evals)) > 1:
                print('unstable poles')
            print('max pole for this time step:', np.amax(np.abs(evals)))
            
            cmd = (- K @ (self.x._get_mpc_x() - x_ref)) #* np.pi/180
            u_storage[idx,:] = cmd
            self.u.values[1:] = cmd
            
            print('u:',self.u.values)
            
            self.step(self.u.values)
            
            Ac,Bc,Cc,Dc = self.linearise(self.x._get_mpc_x(), self.u._get_mpc_u(), _calc_xdot=self._calc_xdot_na, get_obs=self._get_obs_na)
            A,B,C,D = cont2discrete((Ac,Bc,Cc,Dc), self.paras.dt)[0:4]   
            evals, _ = np.linalg.eig(A-B@K)
            
            print('max pole for next time step with old gain:', np.amax(np.abs(evals)))
            
            x_storage[idx,:] = self.x.values
            # bar.update(idx)
            
        vis_x(x_storage, rng)
        vis_u(u_storage, rng)
        
    def test_LQR_lin_all_states(self):
        
        A,B,C,D = self.linearise(x, u)
        
        x = x[:,None]
        u = u[:,None]
        
        A,B,C,D = cont2discrete((A,B,C,D), self.paras.dt)[0:4]
        
        rng = np.linspace(self.paras.time_start, self.paras.time_end, int((self.paras.time_end-self.paras.time_start)/self.paras.dt))
        # bar = progressbar.ProgressBar(maxval=len(rng)).start()
        
        Q = np.eye(len(x)-3)
        R = np.eye(len(u)) * 0.1
        
        Q[0,0] = 0
        Q[1,1] = 0
        Q[5,5] = 0
        
        r_idx = [2,3,4,6,7,8,9,10,11,12,13,14,15,16,17]
        
        Ar = square_mat_degen_2d(A, r_idx)
        Br = np.vstack((B[2:6,:],  B[7:,:]))
        K = dlqr(Ar, Br, Q, R)
        
        # create storage
        x_storage = np.zeros([len(rng),len(self.x.values)])
        u_storage = np.zeros([len(rng),len(self.u.values)])
        
        def get_x_red(x):
            return np.array([x[i,0] for i in r_idx])[:,None]
        
        for idx, val in enumerate(rng):
            
            u = - K @ (get_x_red(x) - get_x_red(x_ref))
            x = A @ x + B @ u
                        
            x_storage[idx,:] = x[:,0]
            u_storage[idx,:] = u[:,0]
            
        vis_x(x_storage, rng)