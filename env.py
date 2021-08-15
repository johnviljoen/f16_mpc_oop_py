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
import os
import progressbar
from scipy.signal import cont2discrete
import scipy
import osqp
from scipy.sparse import csc_matrix

# custom files
from utils import *

class F16(gym.Env):
    
    def __init__(self, state_vector, input_vector, simulation_parameters, nlplant):
        
        super().__init__()
        
        self.x = state_vector                   # mutable state dataclass
        self.u = input_vector                   # mutable input dataclass
        self.paras = simulation_parameters      # immutable simulation parameters dataclass
        self.nlplant = nlplant                  # C interface - the heart of the simulation
        
        # self.action_space = spaces.Box(low=self.u.lower_cmd_bound, high=self.u.upper_cmd_bound, dtype=np.float32)
        # self.observation_space = spaces.Box(low=self.x.lower_bound, high=self.x.upper_bound, shape=(len(self.x.states)), dtype=np.float32)
        
        
        
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
        xdot = np.zeros(18)
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
        self.nlplant.Nlplant(ctypes.c_void_p(x.ctypes.data), ctypes.c_void_p(xdot.ctypes.data), ctypes.c_int(self.paras.fi_flag))    
        #----------assign actuator xdots---------#
        xdot[12:18] = temp
        return xdot
    
    def step(self, action):
        self.x.values += self._calc_xdot(self.x.values, self.u.values)*self.paras.dt
        reward = 1
        isdone = False
        info = {'fidelity':'high'}
        return self.get_obs(self.x.values, self.u.values), reward, isdone, info
    
    def reset(self):
        self.x.values = np.copy(self.x.initial_condition)
        self.u.values = np.copy(self.u.initial_condition)
        return self.get_obs(self.x.values, self.u.values)
    
    def get_obs(self, x, u):
        
        """ Function for acquiring the current observation from the state space.
        
        Args:
            x -> the state vector
            of form numpy 2D array (vertical vector)
            
        Returns:
            y -> system output
            of form numpy 1D array to match gym requirements
        """
        
        return np.copy(np.array(list(x[i] for i in self.x._obs_x_idx), dtype='float32').flatten())
    
    def _calc_xdot_na(self, x, u):
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
        state_vector = np.zeros(18)
        
        for i in range(len(self.x._mpc_x_idx)):
            state_vector[self.x._mpc_x_idx[i]] = x[i]
            
        state_vector[12:16] = u
                                
        # initialise variables
        xdot = np.zeros(18)
        coeff = np.zeros(3)
        C_input_x = np.zeros(18)
        #--------leading edge flap model---------#
        lf_state1_dot, lf_state2_dot = upd_lef(state_vector[2], state_vector[6], coeff, state_vector[7], state_vector[12], state_vector[13], self.nlplant)
        #----------run nlplant for xdot----------#
        # C_input_x = np.concatenate((x[0:12],u,x[13:14]))
        self.nlplant.Nlplant(ctypes.c_void_p(state_vector.ctypes.data), ctypes.c_void_p(xdot.ctypes.data), ctypes.c_int(self.paras.fi_flag))    
        #----------assign actuator xdots---------#
        state_vector_dot = np.concatenate((xdot[0:12],u,np.array([lf_state1_dot, lf_state2_dot])))
        return np.array([xdot[i] for i in self.x._mpc_x_idx])
    
    def _get_obs_na(self, x, u):
        return np.copy(np.array(list(x[i] for i in self.x._mpc_obs_x_idx), dtype='float32').flatten())
    
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
            x[12] = np.clip(x[12], self.u.lower_cmd_bound[0], self.u.upper_cmd_bound[0])
            # elevator limits
            x[13] = np.clip(x[13], self.u.lower_cmd_bound[1], self.u.upper_cmd_bound[1])
            # aileron limits
            x[14] = np.clip(x[14], self.u.lower_cmd_bound[2], self.u.upper_cmd_bound[2])
            # rudder limits
            x[15] = np.clip(x[15], self.u.lower_cmd_bound[3], self.u.upper_cmd_bound[3])
            # alpha limits
            x[7] = np.clip(x[7], self.x.lower_bound[7]*pi/180, self.x.upper_bound[7]*pi/180)
               
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
                
        opt = minimize(obj_func, UX0, args=((h_t, v_t, self.paras.fi_flag, self.nlplant)), method='Nelder-Mead',tol=1e-10,options={'maxiter':5e+04})
        
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
        C = np.zeros([len(self.x.observed_states),len(x)])
        D = np.zeros([len(self.x.observed_states),len(u)])
        
        # Perturb each of the state variables and compute linearization
        for i in range(len(x)):
            
            dx = np.zeros([len(x)])
            dx[i] = eps
                        
            A[:, i] = (_calc_xdot(x + dx, u) - _calc_xdot(x, u)) / eps
            C[:, i] = (get_obs(x + dx, u) - get_obs(x, u)) / eps
            
        # Perturb each of the input variables and compute linearization
        for i in range(len(u)):
            
            du = np.zeros([len(u)])
            du[i] = eps
                    
            B[:, i] = (_calc_xdot(x, u + du) - _calc_xdot(x, u)) / eps
            D[:, i] = (get_obs(x, u + du) - get_obs(x, u)) / eps
        
        return A, B, C, D
    
    def _calc_MPC_action_mk2(self):
        
        hzn = 3
        dt = 0.001
        x = self.x._get_mpc_x()
        u = self.u._get_mpc_u()
        act_states = self.x._get_mpc_act_states()        
        
        A,B,C,D = self.linearise(x, u, _calc_xdot=self._calc_xdot_na, get_obs=self._get_obs_na)
        A,B,C,D = cont2discrete((A,B,C,D), dt)[0:4]
        
        Q = C.T @ C
        R = np.eye(len(u)) * 0.01 # incredibly sensitive
        
        # return A, B, Q, R, hzn, dt, x, act_states,\
        #     self.x._vec_mpc_x_lb,\
        #     self.x._vec_mpc_x_ub,\
        #     self.u._vec_mpc_u_lb,\
        #     self.u._vec_mpc_u_ub,\
        #     self.u._vec_mpc_udot_lb,\
        #     self.u._vec_mpc_udot_ub
        
        OSQP_P, OSQP_q, OSQP_A, OSQP_l, OSQP_u = setup_OSQP(
            A, B, Q, R, hzn, dt, x, u,\
            self.x._vec_mpc_x_lb,\
            self.x._vec_mpc_x_ub,\
            self.u._vec_mpc_u_lb,\
            self.u._vec_mpc_u_ub,\
            self.u._vec_mpc_udot_lb,\
            self.u._vec_mpc_udot_ub)
        
        return OSQP_P, OSQP_q, OSQP_A, OSQP_l, OSQP_u

    def _calc_MPC_action(self, paras_mpc):
        
        hzn = paras_mpc[0]
        dt = 1
        x = self.x._get_mpc_x()
        u = self.u._get_mpc_u()
                
        A,B,C,D = self.linearise(x, u, _calc_xdot=self._calc_xdot_na, get_obs=self._get_obs_na)
        # C[0,3] = 1
        # return A,B,C,D
        A,B,C,D = cont2discrete((A,B,C,D), dt)[0:4]
        
        MM, CC = calc_MC(hzn, A, B, dt)
        
        Q = C.T @ C
        R = np.eye(len(u)) * 0.01 # incredibly sensitive
                
        K = - dlqr(A, B, Q, R)
        Q_bar = scipy.linalg.solve_discrete_lyapunov((A + B @ K).T, Q + K.T @ R @ K)
        QQ = dmom(Q, hzn)
        # return QQ
        RR = dmom(R, hzn)
        QQ[-len(x):,-len(x):] = Q_bar
                
        H = CC.T @ QQ @ CC + RR
        F = CC.T @ QQ @ MM
        G = MM.T @ QQ @ MM
        
        P = 2*H
        q = (2 * x[np.newaxis] @ F.T).T
        
        OSQP_A, OSQP_l, OSQP_u = setup_OSQP_paras(CC, A, x[np.newaxis].T, hzn, len(u), self.x._mpc_x_ub, self.x._mpc_x_lb, self.u._mpc_u_cmd_ub, self.u._mpc_u_cmd_lb, self.u._mpc_u_rate_lb, self.u._mpc_u_rate_lb)
        
        m = osqp.OSQP()
        m.setup(P=csc_matrix(P), q=q, A=csc_matrix(OSQP_A), l=OSQP_l, u=OSQP_u, max_iter=40000, verbose=True)
        res = m.solve()
        
        return res.x[0:len(u)], P, q, OSQP_A, OSQP_l, OSQP_u
    
