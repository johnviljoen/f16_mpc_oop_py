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
from sim import upd_sim
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


# In[]

#----------------------------------------------------------------------------#
#---------------------------------Simulate-----------------------------------#
#----------------------------------------------------------------------------#

def lin_x_traj(hzn, A, B, x0, dt, u_seq):
    x = np.zeros((hzn+1, len(x0)))
    x[0,:] = x0
    for i in range(hzn):
        x[i+1,:] = x[i,:] + (np.matmul(A,x[i,:]) + np.matmul(B, u_seq[i,:])) * dt
    return x[1:]

# def nl_x_traj(hzn, x0, u_seq, nlplant, dt):
#     x_temp = np.zeros([len(x0)])
#     x_seq = np.zeros((hzn,len(x0)))
#     for idx in range(hzn):
#         x_temp = upd_sim(x, u[idx,:], fi_flag, dt, nlplant)
#         x_seq[idx,:] = x_temp
#     return x_seq 

output_vars = [6,7,8,9,10,11]

# trim aircraft
h_t = 10000
v_t = 700

x, opt_res = trim(h_t, v_t, fi_flag, nlplant)

u = x[12:16]

A,B,C,D = linearise(x, u, output_vars, fi_flag, nlplant)

# vert stack
u_seq = np.array([u] * paras_mpc[0])
u_seq_flat = u_seq.reshape(u_seq.shape[0]*u_seq.shape[1],)

x_seq = lin_x_traj(paras_mpc[0], A, B, x, paras_mpc[1], u_seq)

MM, CC = calc_MC(paras_mpc[0], A, B, paras_mpc[1])

x_seq2 = np.matmul(MM, x) + np.matmul(CC, u_seq_flat)

######################TESTING##################



A = np.array([[1.1, 2],[0, 0.95]])
B = np.array([[0],[0.0787]])
C = np.array([-1,1])[np.newaxis]
hzn = 4

Q = np.matmul(C.T, C)
R = 0.01

H, F, G = calc_HFG(A, B, C, hzn, Q, R)

x0 = np.array([0,0])[np.newaxis]

L = -np.matmul(np.linalg.inv(H),F)

K = L[0,:][np.newaxis]

RHS = Q + np.matmul(K.T,K)

clp = A + np.matmul(B, K)



##############################################

# A = np.array([[-2, 1],[0, 1]])
# B = np.array([[1],[1]])
# C = np.array([1, 1])[np.newaxis]

# K = np.array([2, -1])[np.newaxis]

# RHS = np.matmul(C.T,C) + np.matmul(K.T,K)

# clp = A + np.matmul(B,K)

exit()

rng = np.linspace(time_start, time_end, int((time_end-time_start)/time_step))

# create storage
x_storage = np.zeros([len(rng),len(x)])
A = np.zeros([len(x),len(x),len(rng)])
B = np.zeros([len(x),len(u),len(rng)])
C = np.zeros([len(output_vars),len(x),len(rng)])
D = np.zeros([len(output_vars),len(u),len(rng)])

bar = progressbar.ProgressBar(maxval=len(rng)).start()

tic()

for idx, val in enumerate(rng):
    
    #----------------------------------------#
    #------------linearise model-------------#
    #----------------------------------------#
    
    [A[:,:,idx], B[:,:,idx], C[:,:,idx], D[:,:,idx]] = linearise(x, u, output_vars, fi_flag, nlplant)
    
    #----------------------------------------#
    #--------------Take Action---------------#
    #----------------------------------------#
    
    # MPC prediction using squiggly C and M matrices
    CC, MM = calc_MC(paras_mpc[0], A[:,:,idx], B[:,:,idx], time_step)
    
    
    #----------------------------------------#
    #--------------Integrator----------------#
    #----------------------------------------#    
    
    x = upd_sim(x, u, fi_flag, time_step, nlplant)
    
    #----------------------------------------#
    #------------Store History---------------#
    #----------------------------------------#
    
    x_storage[idx,:] = x
    
    bar.update(idx)

toc()

# In[]

#----------------------------------------------------------------------------#
#---------------------------------Visualise----------------------------------#
#----------------------------------------------------------------------------#

#%matplotlib qt

vis(x_storage, rng)

