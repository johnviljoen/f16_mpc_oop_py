#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 22 17:24:58 2021

@author: johnviljoen
"""

from ursina import *
# from ursina.prefabs.first_person_controller import FirstPersonController

class F16_body(Entity):
    def __init__(self, position = (0,0,0)):
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
            model = 'cube',
            scale = (100,3,100),
            color = color.rgb(1,1,1),
            position = (0, -20, 0),
            collider = 'box'
            )

def update():
    if held_keys['esc']:
        application.quit()
    if held_keys['a']:
        test_f16.rotation_y -= 100 * time.dt
    if held_keys['d']:
        test_f16.rotation_y += 100 * time.dt
    if held_keys['w']:
        test_f16.rotation_x += 100 * time.dt
    if held_keys['s']:
        test_f16.rotation_x -= 100 * time.dt
    if held_keys['q']:
        test_f16.rotation_z += 100 * time.dt
    if held_keys['e']:
        test_f16.rotation_z -= 100 * time.dt

app = Ursina()

xscale = 0.1
yscale = 0.1
xpos = 1
ypos = 2


test_f16 = F16_body()

floor = Floor()

camera.orthographic = True
camera.position = (-20, 0, 0)
camera.rotation = (10,0,0)

app.run()

