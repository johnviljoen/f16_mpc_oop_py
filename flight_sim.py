#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 22 17:24:58 2021

@author: johnviljoen
"""

from ursina import *
from ursina.prefabs.first_person_controller import FirstPersonController

from env import F16
from parameters import state_vector, input_vector, simulation_parameters, state_space, nlplant

import numpy as np

f16 = F16(state_vector, input_vector, simulation_parameters, state_space, nlplant)
txt = Text(text="")

class F16_body(Entity):
    def __init__(self):
        super().__init__(
        model = 'f16.obj',
        color = color.white,
        texture = 'white_cube',
        rotation = Vec3(0,0,0),
        scale = (0.01,0.01,0.01)
        )
        
class Floor(Entity):
    def __init__(self):
        super().__init__(
            model = 'plane',
            scale = (100,1,100),
            color = color.rgb(1,1,1),
            position = (0, -10, 0),
            collider = 'box'
            )
        
class Cam(Entity):
    def __init__(self):
        super().__init__()
        camera.orthographic = True
        camera.position = (0, 0, 0)
        camera.rotation = (10,0,0)
        camera.clip_plane_near = 0.0001
        camera.clip_plane_far = 100000
        camera.fov = 45
        camera.z = -10
        
def update():
    dh_cmd = 0
    da_cmd = 0
    dr_cmd = 0
    txt.text = ""
    if held_keys['a']:
        da_cmd = -21.5
        txt.text = 'a'
    if held_keys['d']:
        da_cmd = 21.5
    if held_keys['w']:
        dh_cmd = -25
    if held_keys['s']:
        dh_cmd = 25
    if held_keys['q']:
        dr_cmd = -30
    if held_keys['e']:
        dr_cmd = 30
    Thrust_cmd = 2000
    f16.step(np.array([Thrust_cmd, dh_cmd, da_cmd, dr_cmd]))
    if held_keys['esc']:
        application.quit()
    phi = f16.x.values[3]
    theta = f16.x.values[4]
    psi = f16.x.values[5]
    f16_body.rotation_z = phi
    f16_body.rotation_y = theta
    f16_body_rotation_x = psi
    # rotate_entity(f16_body)
    # translate_entity(camera)

app = Ursina()

def rotate_entity(entity):
    if held_keys['a']:
        entity.rotation_y -= 100 * time.dt
    if held_keys['d']:
        entity.rotation_y += 100 * time.dt
    if held_keys['w']:
        entity.rotation_x += 100 * time.dt
    if held_keys['s']:
        entity.rotation_x -= 100 * time.dt
    if held_keys['q']:
        entity.rotation_z += 100 * time.dt
    if held_keys['e']:
        entity.rotation_z -= 100 * time.dt
        
def translate_entity(entity):
    if held_keys['j']:
        entity.y -= 100 * time.dt
    if held_keys['l']:
        entity.y += 100 * time.dt
    if held_keys['i']:
        entity.x += 100 * time.dt
    if held_keys['k']:
        entity.x -= 100 * time.dt
    if held_keys['u']:
        entity.z += 100 * time.dt
    if held_keys['o']:
        entity.z -= 100 * time.dt

f16_body = F16_body()

floor = Floor()

cam = Cam()

f16.paras.dt = 1/60



app.run()

