#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  7 14:50:18 2021

@author: johnviljoen
"""

import numpy as np
from numpy import pi

from parameters import act_lim

import ctypes

def upd_thrust(T_cmd, T_state):
    # command saturation
    T_cmd = np.clip(T_cmd,act_lim[1][0],act_lim[0][0])
    # rate saturation
    T_err = np.clip(T_cmd - T_state, -10000, 10000)
    return T_err

def upd_dstab(dstab_cmd, dstab_state):
    # command saturation
    dstab_cmd = np.clip(dstab_cmd,act_lim[1][1],act_lim[0][1])
    # rate saturation
    dstab_err = np.clip(20.2*(dstab_cmd - dstab_state), -60, 60)
    return dstab_err

def upd_ail(ail_cmd, ail_state):
    # command saturation
    ail_cmd = np.clip(ail_cmd,act_lim[1][2],act_lim[0][2])
    # rate saturation
    ail_err = np.clip(20.2*(ail_cmd - ail_state), -80, 80)
    return ail_err

def upd_rud(rud_cmd, rud_state):
    # command saturation
    rud_cmd = np.clip(rud_cmd,act_lim[1][3],act_lim[0][3])
    # rate saturation
    rud_err = np.clip(20.2*(rud_cmd - rud_state), -120, 120)
    return rud_err

def upd_lef(h, V, coeff, alpha, lef_state_1, lef_state_2, nlplant):
    
    nlplant.atmos(ctypes.c_double(h),ctypes.c_double(V),ctypes.c_void_p(coeff.ctypes.data))
    atmos_out = coeff[1]/coeff[2] * 9.05
    alpha_deg = alpha*180/pi
    
    LF_err = alpha_deg - (lef_state_1 + (2 * alpha_deg))
    #lef_state_1 += LF_err*7.25*time_step
    LF_out = (lef_state_1 + (2 * alpha_deg)) * 1.38
    
    lef_cmd = LF_out + 1.45 - atmos_out
    
    # command saturation
    lef_cmd = np.clip(lef_cmd,act_lim[1][4],act_lim[0][4])
    # rate saturation
    lef_err = np.clip((1/0.136) * (lef_cmd - lef_state_2),-25,25)
    # integrate
    #lef_state_2 += lef_err*time_step
    
    return LF_err*7.25, lef_err

def calc_xdot(x, u, fi_flag, nlplant):
    
    # initialise variables
    xdot = np.zeros(18)
    temp = np.zeros(6)
    coeff = np.zeros(3)
    
    #--------------Thrust Model--------------#
    temp[0] = upd_thrust(u[0], x[12])
    #--------------Dstab Model---------------#
    temp[1] = upd_dstab(u[1], x[13])
    #-------------aileron model--------------#
    temp[2] = upd_ail(u[2], x[14])
    #--------------rudder model--------------#
    temp[3] = upd_rud(u[3], x[15])
    #--------leading edge flap model---------#
    temp[5], temp[4] = upd_lef(x[2], x[6], coeff, x[7], x[17], x[16], nlplant)
    
    #----------run nlplant for xdot----------#
    nlplant.Nlplant(ctypes.c_void_p(x.ctypes.data), ctypes.c_void_p(xdot.ctypes.data), ctypes.c_int(fi_flag))    
    
    xdot[12:18] = temp
    
    return xdot

def upd_sim(x, u, fi_flag, time_step, nlplant):
    
    # find xdot
    xdot = calc_xdot(x, u, fi_flag, nlplant)
    
    # update x
    x += xdot*time_step
    
    return x

def calc_out(x, u, output_vars):
    # return the variables    
    return x[output_vars]
