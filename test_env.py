#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug  1 20:09:42 2021

@author: johnviljoen
"""

import unittest

from env import F16
from parameters import initial_state_vector_ft_rad as x0
from parameters import simulation_parameters as paras_sim

class test_F16(unittest.TestCase, F16):
    
    def __init__(self, x0, u0, paras_sim):
        super(unittest.TestCase, self).__init__(x0, u0, paras_sim)
    
    def test_upd_thrust(self):
        # self.assertAlmostEqual()
        pass
    
    
test_f16 = test_F16(x0, x0[12:16], paras_sim)