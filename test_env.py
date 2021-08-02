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

from control.matlab import *
from parameters import act_lim, x_lim

class test_F16(unittest.TestCase, F16):
    
    def __init__(self, x0, u0, paras_sim):
        super(unittest.TestCase, self).__init__(x0, u0, paras_sim)
    
    def test_upd_thrust(self):
        # self.assertAlmostEqual()
        """ Function to test that the thrust model is behaviing as expected. This
        is done by creating a seperate control system modelled in a different way
        and then verifying the outputs are the same
        
        The model is first order and so it is also simple to test the command
        and rate saturations to ensure accuracy"""
        pass
        
        
        
    def test_act_cmd_lims(self):
        
        """ Function to test that the actuator command limits are implemented correctly by
        testing the boundaries of their operating envelope. It should be noted that
        this test is conducted simulatenously on all actuators, and therefore 
        it assumes that the actuator dynamics are not coupled, which they shouldnt
        be... it is fairly simple to see they arent in the code written in "env.py". 
                
        The engine is tested by commanding it to go 1000 below and 1000 above
        its maximum value whilst having the current engine state be 1 below and 
        1 above its maximum value respectively. This test is successful if  and 
        only if this command is seen to be ignored by the generation of an engine 
        state time derivative in the opposite direction. 
        
        The other actuators will crash the C code if they go above or below their 
        maximum and minimum values respectively. Therefore their states are set
        at their maximums and minimums and they are commanded to go beyond them.
        The test is successful if and only if their rate of change is found to be
        zero"""
        
        self.reset()
        
        # set command 1000 below the minimum
        self.u[0] = act_lim[1][0] - 1000
        self.u[1] = act_lim[1][1] - 1000
        self.u[2] = act_lim[1][2] - 1000
        self.u[3] = act_lim[1][3] - 1000
        
        # set the current engine state 1 below the minimum
        self.x[12] = act_lim[1][0] - 1
        
        # the C code crashes if ordered to lookup data beyond its tables for the
        # other actuators, therefore dh, da, dr are set to their minimum exactly
        # and it is tested if a command can decrease this
        self.x[13] = act_lim[1][1] 
        self.x[14] = act_lim[1][2]
        self.x[15] = act_lim[1][3]
        
        # calculate the state vector time derivative
        xdot = self.calc_xdot(self.x, self.u)
        
        # check that the rate of change of the engine state is indeed positive
        self.assertGreater(xdot[12], 0)
        
        # check the rate of change of dh, da, dr is indeed zero
        self.assertAlmostEqual(xdot[13], 0)
        self.assertAlmostEqual(xdot[14], 0)
        self.assertAlmostEqual(xdot[15], 0)

        # check the maximums now, but first reset the simulation.
        self.reset()
        
        # set command 1000 above the maximum
        self.u[0] = act_lim[0][0] + 1000
        self.u[1] = act_lim[0][1] + 1000
        self.u[2] = act_lim[0][2] + 1000
        self.u[3] = act_lim[0][3] + 1000
        
        # set the current engine state 1 above the maximum
        self.x[12] = act_lim[0][0] + 1
        
        # set dh, da, dr to their maximums
        self.x[13] = act_lim[0][1] 
        self.x[14] = act_lim[0][2]
        self.x[15] = act_lim[0][3]
        
        # calculate the state vector time derivative
        xdot = self.calc_xdot(self.x, self.u)
        
        # check that the rate of change of the engine state is indeed negative
        self.assertLess(xdot[12], 0)
        
        # check that the rate of change of dh, da, dr is indeed zero
        self.assertAlmostEqual(xdot[13], 0)
        self.assertAlmostEqual(xdot[14], 0)
        self.assertAlmostEqual(xdot[15], 0)
        
    def test_act_rate_lims(self):
        
        """ Function to test the rate limits of the 1st order actuators is behaving as 
        expected. Note this does not include the engine as it is a more complex system."""
        
        # begin from rough trim
        self.reset()
        
        # command maximums on all actuators
        self.u[1] = act_lim[0][1]
        self.u[2] = act_lim[0][2]
        self.u[3] = act_lim[0][3]
        
        xdot = self.calc_xdot(self.x, self.u)
        
        self.assertAlmostEqual(xdot[13], 60)
        self.assertAlmostEqual(xdot[14], 80)
        self.assertAlmostEqual(xdot[15], 120)
        
        self.reset()
                        
        # now test the inverse
        self.u[1] = act_lim[1][1]
        self.u[2] = act_lim[1][2]
        self.u[3] = act_lim[1][3]
                
        xdot = self.calc_xdot(self.x, self.u)

        self.assertAlmostEqual(xdot[13], -60)
        self.assertAlmostEqual(xdot[14], -80)
        self.assertAlmostEqual(xdot[15], -120)
        
    def test_aerodynamics(self):
        
        pass
    
    
                    
    
test_f16 = test_F16(x0, x0[12:16], paras_sim)
test_f16.test_act_rate_lims()