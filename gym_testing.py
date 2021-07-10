#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 10 13:54:33 2021

@author: johnviljoen
"""

from parameters import x_lim, act_lim

import numpy as np

import gym
from gym import spaces

class F16_env(gym.Env):
  """Custom Environment that follows gym interface"""
  
  metadata = {'render.modes': ['human']}

  def __init__(self, arg1, arg2):
    super(F16_env, self).__init__()
    # Define action and observation space
    # They must be gym.spaces objects
    # Example when using discrete actions:
    self.action_space = spaces.Box(low=np.array(act_lim[1])[0:4], high=np.array(act_lim[0])[0:4])
    
    self.observation_space = spaces.Box(low=np.array(x_lim[1]), high=np.array(x_lim[0]))
    #self.action_space = spaces.Box(low=act_lim[1], high=act_lim[0], shape=(np.array([len(act_lim[0]),1,1])), dtype=np.float64)
    # Example for using image as input:
    #self.observation_space = spaces.Box(low=x_lim[1], high=x_lim[0] shape=(len(x_lim[0]), 1, 1), dtype=np.float64)

  # def step(self, action):
  #   # Execute one time step within the environment
  #   ...
  # def reset(self):
  #   # Reset the state of the environment to an initial state
  #   ...
  # def render(self, mode='human', close=False):
  #   # Render the environment to the screen
  #   ...