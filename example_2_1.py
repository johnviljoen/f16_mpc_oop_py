#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 10 16:23:19 2021

@author: johnviljoen
"""
import numpy as np
from utils import *

A = np.array([[1.1,2],[0,0.95]])
B = np.array([[0],[0.0787]])
C = np.array([-1,1])[np.newaxis]

ninputs = B.shape[1]
nstates = A.shape[0]

hzn = 4
dt = 1

MM, CC = calc_MC(hzn, A, B, dt)

Q = C.T @ C
R = np.eye(ninputs) * 0.01
Q_bar = Q

# turn into diagonals
QQ = dmom(Q, hzn)
RR = dmom(R, hzn)

H = CC.T @ QQ @ CC + RR
F = CC.T @ QQ @ MM
G = MM.T @ QQ @ MM

