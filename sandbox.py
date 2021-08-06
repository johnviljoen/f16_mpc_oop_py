#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  6 08:06:11 2021

@author: johnviljoen
"""

from scipy.optimize import minimize
from scipy.signal import cont2discrete
import numpy as np

from env import F16

from parameters import simulation_parameters as paras_sim
from parameters import initial_state_vector_ft_rad as x0

from utils import tic, toc, vis, dmom, square_mat_degen_2d


f16 = F16(x0, x0[12:16], paras_sim)

A,B,C,D = f16.linearise(f16.x, f16.u)
A,B,C,D = cont2discrete((A,B,C,D), f16.dt)[0:4]

A_degen = square_mat_degen_2d(A, f16.x_degen_idx)
B_degen = B[f16.x_degen_idx,1:4]

