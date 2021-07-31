#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 18 14:14:34 2021

@author: johnviljoen
"""

from parameters import initial_state_vector_ft_rad as x0
from parameters import simulation_parameters as paras_sim

import gym

from stable_baselines3 import A2C
from stable_baselines3.common.vec_env import DummyVecEnv
from stable_baselines3.common.policies import ActorCriticPolicy

from F16_gym import F16_env

env = DummyVecEnv([lambda: F16_env(x0, paras_sim)])

model = A2C('MlpPolicy', env, verbose = 1)
# model.learn(total_timesteps=10000)

# obs = env.reset()
# for i in range(10000):
#     action, _states = model.predict(obs)
#     obs, rewards, isdone, info = env.step(action)