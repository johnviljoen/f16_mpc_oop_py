#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 31 19:07:44 2021

@author: johnviljoen
"""

from env import F16

from parameters import initial_state_vector_ft_rad as x0
from parameters import simulation_parameters as paras_sim

import gym

from stable_baselines3 import A2C
from stable_baselines3.common.vec_env import DummyVecEnv
from stable_baselines3.common.policies import ActorCriticPolicy

from stable_baselines3.common.env_checker import check_env

from sys import exit

f16 = F16(x0, x0[12:16], paras_sim)

f16.validate_sim(f16.x0)

check_env(f16, warn=True)

exit()


env = DummyVecEnv([lambda: F16(x0, x0[12:16], paras_sim)])

model = A2C('MlpPolicy', env, verbose=1)

# model.learn(total_timesteps=100)