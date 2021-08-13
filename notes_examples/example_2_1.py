#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 10 16:23:19 2021

@author: johnviljoen
"""
import numpy as np
from utils import *
from scipy.signal import cont2discrete

f16 = False

hzn = 4
dt = 1

if f16:
    
    npzfile = np.load('f16SS.npz')
    A = npzfile['A']
    B = npzfile['B']
    C = npzfile['C']
    D = npzfile['D']
    A,B,C,D = cont2discrete((A,B,C,D), dt)[0:4]
    
else:

    A = np.array([[1.1,2],[0,0.95]])
    B = np.array([[0],[0.0787]])
    C = np.array([-1,1])[np.newaxis]

ninputs = B.shape[1]
nstates = A.shape[0]



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

