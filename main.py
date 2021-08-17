#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 31 19:07:44 2021

@author: johnviljoen
"""

# dependencies
from stable_baselines3 import A2C
from stable_baselines3.common.vec_env import DummyVecEnv
from scipy.signal import cont2discrete
import numpy as np
from scipy.signal import cont2discrete
from stable_baselines3.common.env_checker import check_env
from sys import exit
from scipy.sparse import csc_matrix
import osqp

# custom files
from env import F16
from parameters import state_vector, input_vector, simulation_parameters, nlplant

f16 = F16(state_vector, input_vector, simulation_parameters, nlplant)

# res = f16._calc_MPC_action(2,0,0,10)






exit()


env = DummyVecEnv([lambda: F16(x0, x0[12:16], paras_sim)])

model = A2C('MlpPolicy', env, verbose=1)

# model.learn(total_timesteps=100)