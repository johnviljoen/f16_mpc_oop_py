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
from scipy.signal import cont2discrete
import scipy
import osqp
from scipy.sparse import csc_matrix

# custom files
from parameters import x_lim, u_cmd_lim, u_rate_lim
from utils import *

class F16(gym.Env):
    
    def __init__(self, x0, u0, paras_sim):
        
        super().__init__()
        
        self.x = np.copy(x0)                    # dynamics state
        self.u = np.copy(u0)                    # dynamic input
        
        self.x0 = np.copy(x0)                   # initial state
        self.u0 = np.copy(u0)                   # initial input
        
        self.dt = paras_sim[0]                  # simulation time step
        self.time_start = paras_sim[1]          # simulation start time
        self.time_end = paras_sim[2]            # simulation end time
        self.stab_flag = paras_sim[3]           # stability flag
        self.fi_flag = paras_sim[4]             # fidelity flag
        self.x_lim = paras_sim[5]               # state limits
        self.u_cmd_lim = paras_sim[6]           # command limit saturation
        self.u_rate_lim = paras_sim[7]          # command rate saturation
        self.observed_states = paras_sim[8]     # observed states -> generates y
        self.nlplant = paras_sim[9]             # C function interface

        self.action_space = spaces.Box(low=np.array(self.u_cmd_lim[1])[0:4], high=np.array(self.u_cmd_lim[0])[0:4], dtype=np.float32)
        self.observation_space = spaces.Box(low=np.array(list((self.x_lim[1] + self.u_cmd_lim[1])[i] for i in self.observed_states), dtype='float32'),\
                                            high=np.array(list((self.x_lim[0] + self.u_cmd_lim[0])[i] for i in self.observed_states), dtype='float32'),\
                                                shape=(len(self.observed_states),), dtype=np.float32)
        
    def _calc_xdot(self, x, u):
        
        """ calculates, and returns the rate of change of the state vector, x, using the empirical
        aerodynamic data held in folder 'C', also using equations of motion found in the
        shared library C file. This function includes all actuator models.
        
        Args:
            x:
                numpy 2D array (vertical vector) of 18 elements
                {xe,ye,h,phi,theta,psi,V,alpha,beta,p,q,r,T,dh,da,dr,lf2,lf1}
            u:
                numpy 2D array (vertical vector) of 4 elements
                {T,dh,da,dr}
    
        Returns:
            xdot:
                numpy 2D array (vertical vector) of 18 elements
                time derivatives of {xe,ye,h,phi,theta,psi,V,alpha,beta,p,q,r,T,dh,da,dr,lf2,lf1}
        """
        
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
        self.x += self._calc_xdot(self.x, self.u)*self.dt
        reward = 1
        isdone = False
        info = {'fidelity':'high'}
        return self.get_obs(self.x, self.u), reward, isdone, info
    
    def reset(self):
        self.x = np.copy(self.x0)
        self.u = np.copy(self.u0)
        return self.get_obs(self.x, self.u)
        
    def get_obs(self, x, u):
        
        """ Function for acquiring the current observation from the state space.
        
        Args:
            x -> the state vector
            of form numpy 2D array (vertical vector)
            
        Returns:
            y -> system output
            of form numpy 1D array to match gym requirements
        """
        
        return np.copy(np.array(list(x[i] for i in self.observed_states), dtype='float32').flatten())
    
    def _sim(self, x0, visualise=True):
        """ 
        Args:
            x0:
                initial state vector, 2D numpy array (vertical vector)
                
        Returns:
            x_storage:
                timehistory sequence of state vectors, 2D numpy array
        """
        rng = np.linspace(self.time_start, self.time_end, int((self.time_end-self.time_start)/self.dt)) # setup time sequence
        x_storage = np.zeros([len(rng),len(self.x)]) # initialise storage
        bar = progressbar.ProgressBar(maxval=len(rng)).start() # begin progressbar
        self.x = x0 # reset x
        tic() # begin timer
        for idx, val in enumerate(rng):
            # MPC prediction using squiggly C and M matrices
            #CC, MM = calc_MC(paras_mpc[0], A[:,:,idx], B[:,:,idx], time_step)
            # self.u[1:4] = self.calc_MPC_action_mk2(2,2,2,[10,0.001])
            self.step(self.u) # integrator
            x_storage[idx,:] = self.x[:,0] # store history
            bar.update(idx) # update progressbar
        toc() # finish timer
        if visualise:
            # run this in spyder terminal to have plots appear in standalone windows
            # %matplotlib qt
            vis(x_storage, rng)
        return x_storage
    
    # NEEDS WORK
    def _sim_na(self, visualise=True):
        rng = np.linspace(self.time_start, self.time_end, int((self.time_end-self.time_start)/self.dt)) # setup time sequence
        x_storage = np.zeros([len(rng),len(self.x_na)]) # create storage
        bar = progressbar.ProgressBar(maxval=len(rng)).start() # begin progress bar
        tic() # start timer
        for idx, val in enumerate(rng):
            self.linearise(self.x_na, self.u, _calc_xdot=self._calc_xdot_na, get_obs=self.get_obs_na)
            self.x_na += self._calc_xdot_na(self.x_na, self.u) * self.dt
            x_storage[idx,:] = self.x_na[:,0]
            bar.update(idx)
        toc() # end timer
    
    def __calc_xdot_na(self, x, u):
        """
        Args:
            x:
                {h,phi,theta,V,alpha,beta,p,q,r,lf1,lf2}
                numpy 2D array (vertical vector) of 14 elements
                {xe,ye,h,phi,theta,psi,V,alpha,beta,p,q,r,lf1,lf2}
            u:
                numpy 2D array (vertical vector) of 4 elements
                {T,dh,da,dr}
    
        Returns:
            xdot:
                numpy 2D array (vertical vector) of 14 elements
                time derivatives of {xe,ye,h,phi,theta,psi,V,alpha,beta,p,q,r,lf1,lf2}
        """
        x = np.array([0, 0, x[0,0], x[1,0], x[2,0], 0, x[3,0] ,x[4,0], x[5,0], x[6,0], x[7,0], x[8,0], x[9,0], x[10,0]])[np.newaxis].T
        # initialise variables
        xdot = np.zeros([18,1])
        coeff = np.zeros(3)
        C_input_x = np.zeros(18)
        #--------leading edge flap model---------#
        lf_state1_dot, lf_state2_dot = upd_lef(x[2], x[6], coeff, x[7], x[12], x[13], self.nlplant)
        #----------run nlplant for xdot----------#
        C_input_x = np.concatenate((x[0:12],u,x[13:14]))
        self.nlplant.Nlplant(ctypes.c_void_p(C_input_x.ctypes.data), ctypes.c_void_p(xdot.ctypes.data), ctypes.c_int(self.fi_flag))    
        #----------assign actuator xdots---------#
        return np.concatenate((xdot[2:5], xdot[6:12], np.array([lf_state1_dot, lf_state2_dot])))
    
    def _get_obs_na(self, x, u):
        return np.copy(np.array(list(x[i] for i in list(np.array(self.observed_states) - 3)), dtype='float32').flatten())
        
    def calc_MPC_action_pqr(self, p_dem, q_dem, r_dem, paras_mpc):
        
        hzn = paras_mpc[0] * 10
        # dt = self.dt
        dt = 1
        x = np.concatenate((self.x[2:5],self.x[6:12],self.x[16:17],self.x[15:16])) # select the new state vector -> removing actuators
        x_lim = [list(self.x_lim[0][i] for i in [2,3,4,6,7,8,9,10,11]) , list(self.x_lim[1][i] for i in [2,3,4,6,7,8,9,10,11])]
        u = np.copy(self.u[1:4]) # we ignore the engine actuator here as we only command p,q,r
        A,B,C,D = self.linearise(x, u, _calc_xdot=self.__calc_xdot_na, get_obs=self._get_obs_na)
        A,B,C,D = cont2discrete((A,B,C,D), dt)[0:4]
        ignore_states = [0,1,5]
        ninputs = len(u)
        nstates = len(x)
        MM, CC = calc_MC(hzn, A, B, dt)
        
        Q = C.T @ C
        R = np.eye(ninputs) * 1 # incredibly sensitive
        
        # return A, B, Q, R
        
        K = - dlqr(A, B, Q, R)
        Q_bar = scipy.linalg.solve_discrete_lyapunov((A + B @ K).T, Q + K.T @ R @ K)
        QQ = dmom(Q, hzn)
        RR = dmom(R, hzn)
        QQ[-nstates:,-nstates:] = Q_bar
        
        H = CC.T @ QQ @ CC + RR
        F = CC.T @ QQ @ MM
        G = MM.T @ QQ @ MM
        
        P = 2*H
        q = (2 * x.T @ F.T).T
        
        OSQP_A, OSQP_l, OSQP_u = setup_OSQP_paras(CC, A, x, hzn, ninputs, x_lim, self.u_cmd_lim, self.u_rate_lim)
        
        m = osqp.OSQP()
        m.setup(P=csc_matrix(P), q=q, A=csc_matrix(OSQP_A), l=OSQP_l, u=OSQP_u, max_iter=40000, verbose=True)
        res = m.solve()
        
        return res.x[0:ninputs], P, q, OSQP_A, OSQP_l, OSQP_u
        
    def trim(self, h_t, v_t):
        
        """ Function for trimming the aircraft in straight and level flight. The objective
        function is built to be the same as that of the MATLAB version of the Nguyen 
        simulation.
        
        Args:
            h_t:
                altitude above sea level in ft, float
            v_t:
                airspeed in ft/s, float
                
        Returns:
            x_trim:
                trim state vector, 1D numpy array
            opt:
                scipy.optimize.minimize output information
        """
        
        def obj_func(UX0, h_t, v_t, fi_flag, nlplant):
    
            V = v_t
            h = h_t
            P3, dh, da, dr, alpha = UX0
            npos = 0
            epos = 0
            phi = 0
            psi = 0
            beta = 0
            p = 0
            q = 0
            r = 0
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
            
            # thrust limits
            x[12] = np.clip(x[12], self.u_cmd_lim[1][0], self.u_cmd_lim[0][0])
            # elevator limits
            x[13] = np.clip(x[13], self.u_cmd_lim[1][1], self.u_cmd_lim[0][1])
            # aileron limits
            x[14] = np.clip(x[14], self.u_cmd_lim[1][2], self.u_cmd_lim[0][2])
            # rudder limits
            x[15] = np.clip(x[15], self.u_cmd_lim[1][3], self.u_cmd_lim[0][3])
            # alpha limits
            x[7] = np.clip(x[7], self.x_lim[1][7]*pi/180, self.x_lim[0][7]*pi/180)
               
            u = np.array([x[12],x[13],x[14],x[15]])
            xdot = self._calc_xdot(x, u)
            
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
        
    def linearise(self, x, u, _calc_xdot=None, get_obs=None):
        
        """ Function to linearise the aircraft at a given state vector and input demand.
        This is done by perturbing each state and measuring its effect on every other state.
        
        Args:
            x:
                state vector, 2D numpy array (vertical vector)
            u:
                input vector, 2D numpy array (vertical vector)
                
        Returns:
            4 2D numpy arrays, representing the 4 state space matrices, A,B,C,D.
        """
        
        if _calc_xdot == None:
            _calc_xdot = self._calc_xdot
        if get_obs == None:
            get_obs = self.get_obs
        
        eps = 1e-06
        
        A = np.zeros([len(x),len(x)])
        B = np.zeros([len(x),len(u)])
        C = np.zeros([len(self.observed_states),len(x)])
        D = np.zeros([len(self.observed_states),len(u)])
        
        # Perturb each of the state variables and compute linearization
        for i in range(len(x)):
            
            dx = np.zeros([len(x),1])
            dx[i] = eps
            
            A[:, i] = (_calc_xdot(x + dx, u)[:,0] - _calc_xdot(x, u)[:,0]) / eps
            C[:, i] = (get_obs(x + dx, u) - get_obs(x, u)) / eps
            
        # Perturb each of the input variables and compute linearization
        for i in range(len(u)):
            
            du = np.zeros([len(u),1])
            du[i] = eps
                    
            B[:, i] = (_calc_xdot(x, u + du)[:,0] - _calc_xdot(x, u)[:,0]) / eps
            D[:, i] = (get_obs(x, u + du) - get_obs(x, u)) / eps
        
        return A, B, C, D   