#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 31 19:07:44 2021

@author: johnviljoen
"""

from env_mk4 import F16

from parameters_mk2 import state_vector, input_vector, simulation_parameters, nlplant

from stable_baselines3 import A2C
from stable_baselines3.common.vec_env import DummyVecEnv

from scipy.signal import cont2discrete

import numpy as np
from scipy.signal import cont2discrete

from stable_baselines3.common.env_checker import check_env

from sys import exit

f16 = F16(state_vector, input_vector, simulation_parameters, nlplant)
# u_opt = f16.calc_MPC_action_mk2(10,10,10,paras_mpc)
# A,B,C,D = f16.linearise(f16.x, f16.u)
# A,B,C,D = cont2discrete((A,B,C,D),)
# res, P, q, OSQP_A, OSQP_l, OSQP_u = f16.calc_MPC_action_pqr(0, 0, 0, paras_mpc)




exit()


env = DummyVecEnv([lambda: F16(x0, x0[12:16], paras_sim)])

model = A2C('MlpPolicy', env, verbose=1)

# model.learn(total_timesteps=100)