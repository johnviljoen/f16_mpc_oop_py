#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  9 19:40:12 2021

@author: johnviljoen
"""

# In[] imports

# from ctypes import *
from ctypes import CDLL
#import ctypes
import os

# import numpy and sin, cos for convenience
import numpy as np

# handbuilt functions for all this
from utils import tic, toc, vis
from trim import trim
from sim import upd_sim, calc_xdot
from mpc import linearise, dmom, calc_HFG, calc_MC

# import progressbar for convenience
import progressbar

# import parameters
from parameters import initial_state_vector_ft_rad, simulation_parameters, paras_mpc

# import exit() function for debugging
from sys import exit

# In[]

#----------------------------------------------------------------------------#
#-------------------------prepare data for nlplant.c-------------------------#
#----------------------------------------------------------------------------#

# unwrap simulation parameters
time_step, time_start, time_end, stab_flag, fi_flag = simulation_parameters

# create interface with c shared library .so file in folder "C"
if stab_flag == 1:
    so_file = os.getcwd() + "/C/nlplant_xcg35.so"
elif stab_flag == 0:
    so_file = os.getcwd() + "/C/nlplant_xcg25.so"
    
nlplant = CDLL(so_file)

# initialise x
x = initial_state_vector_ft_rad

output_vars = [6,7,8,9,10,11]

# trim aircraft
h_t = 10000
v_t = 700

x, opt_res = trim(h_t, v_t, fi_flag, nlplant)

u = x[12:16]

A,B,C,D = linearise(x, u, output_vars, fi_flag, nlplant)

# In[]

# Import do_mpc package:
import do_mpc

model_type = 'discrete' # either 'discrete' or 'continuous'
model = do_mpc.model.Model(model_type)

do_mpc.controller.MPC

import casadi 
casadi.casadi.Function

# doesnt seem compatible with my simulation unfortunately
