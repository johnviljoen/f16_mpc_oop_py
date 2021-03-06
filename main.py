#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 31 19:07:44 2021

@author: johnviljoen
"""

# dependencies
# from stable_baselines3 import A2C
# from stable_baselines3.common.vec_env import DummyVecEnv
from scipy.signal import cont2discrete
import numpy as np
# from stable_baselines3.common.env_checker import check_env
from sys import exit
from scipy.sparse import csc_matrix
import osqp
import ctypes
import scipy.io

# custom files
from env import F16
from test_env_mk2 import test_F16
from parameters import state_vector, input_vector, simulation_parameters, state_space, nlplant

f16 = F16(state_vector, input_vector, simulation_parameters, state_space, nlplant)
test_f16 = test_F16(state_vector, input_vector, simulation_parameters, state_space, nlplant)

# 1 means visualise in flight sim, 0 means dont
flight_sim = False

# test_f16.test_LQR_lin(f16=False)
# test_f16.test_LQR_static_nl()

test_f16.LQR(linear=True)

# test_f16.offline_LQR_nl()
# test_f16.SSR_continuous_PID_lin()

# mat = scipy.io.loadmat("MATLAB_SS.mat")

# test_f16.test_linearisation()

# run the flight sim
if flight_sim:
    exec(open('flight_sim.py').read())

exit()


env = DummyVecEnv([lambda: F16(x0, x0[12:16], paras_sim)])

model = A2C('MlpPolicy', env, verbose=1)

# model.learn(total_timesteps=100)