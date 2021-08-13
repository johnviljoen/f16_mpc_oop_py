#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 10 16:35:35 2021

@author: johnviljoen
"""

from example_2_1 import *

out = -np.linalg.inv(H) @ F

L = out[0,:]


K = - dlqr(A, B, Q, R)

from scipy.linalg import solve_discrete_lyapunov

Q_bar = solve_discrete_lyapunov((A + np.matmul(B, K)).T, Q + np.matmul(np.matmul(K.T,R), K))
Q_bar = solve_discrete_lyapunov((A + B @ K).T, Q + K.T @ R @ K)

QQ[-A.shape[0]:,-A.shape[0]:] = Q_bar

H = CC.T @ QQ @ CC + RR
F = CC.T @ QQ @ MM
G = MM.T @ QQ @ MM

L = - np.linalg.inv(H) @ F


