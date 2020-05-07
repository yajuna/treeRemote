#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 27 20:03:51 2020

ref: https://scicomp.stackexchange.com/questions/20279/how-can-i-solve-wave-equation-for-circular-membrane-in-polar-coordinates

visualization ref with mayavi: https://docs.enthought.com/mayavi/mayavi/installation.html

This code solves equation (10) in Potter and Andresen

To Fix: 
    some issues with visualization
    add source term
    justify boundary and initial conditions
    Assume rho, c, k are all constants
    Check CFL condition
    the role of jn_zero and jv
    vectorize to increase efficiency (can wait till later)

@author: yajun
"""
import numpy as np
from scipy.special import jv, jn_zeros
import matplotlib.pyplot as plt
from matplotlib import cm

#from mayavi import mlab 

#%% Parameters; play with parameters, and see if warnings would go away
Nr = 50  
N_phi = 50  
N_steps = 100 
radius = 5.  
c = 50#4 # = k/(rho*c) in equation (2) of Potter and Anredsen paper
dphi = 2*np.pi/N_phi  
dphi2 = dphi * dphi
dr = radius/Nr 
dr2 = dr * dr 
dt = 0.001   
if (1 - 4*dt*c/dr2 - 4*dt*c/(dr2*dphi2))**2<1: # CFL condition. 
    # The maximum value for dt is dr*dphi/(2*c)
    print("CFL condition is satisfied")
else: 
    print("CFL condition not satisfied, shrink dt by 0.0001")
    dt = dt*0.0001
#%% Initial conditions; waiting for data
r = np.linspace(0, radius, Nr)  
phi = np.linspace(0, 2*np.pi, N_phi)  
R, Phi = np.meshgrid(r, phi) 
X = R*np.cos(Phi) 
Y = R*np.sin(Phi) 
Z = np.cos(Phi) * (R-R**3)
print("initial condition sin(phi)*(r-r^3)")
# T has to be float, can not make it int
T = np.zeros((N_steps, Nr, N_phi))
T[0, :, :] = Z.T 
#%% visualize initial condition in 3D



#%% Stepping
k1 = c*dt/dr2
for t in range(1, N_steps):
    for i in range(0, Nr-1):
        ri = max(r[i], 0.5*dr)  # To avoid the singularity at r=0
        k2 = c*dt/(2*ri*dr)
        k3 = c*dt/(dphi2*ri**2)
#        print(k2,k3,i)
        for j in range(0, N_phi-1):
            T[t, i, j] = T[t-1, i, j] \
            + k1*(T[t-1, i+1, j] - 2*T[t-1, i, j] + T[t-1, i-1, j])\
            + k2*(T[t-1, i+1, j] - T[t-1, i-1, j])\
            + k3*(T[t-1, i, j+1] - 2*T[t-1, i, j] + T[t-1, i, j-1])#\
#            + source terms add source terms here
        
        T[t, :, -1] = T[t, :, 0]  # Update the values for phi=2*pi, BC in phi
#        T[t,0,i]?? T[t,R,i]??
    print(np.max(abs(T[t,:,:])),t)
        
#%% visualization with matplolib BUGGY, NEED TO FIX
#fig = plt.figure()
#ax = fig.gca(projection='3d')
#
#u,v = np.meshgrid(X.T, Y.T)        
#surf = ax.plot_surface(u, v, 10*T[999], cmap = cm.coolwarm, linewidth = 0, antialiased = False)
#fig.colorbar(surf, shrink = 0.5, aspect = 5)
#plt.show()       
   
"""   
boundary condition for theta:    T[t, i, -1] = T[t, i, 0]
boundary condition for radius direction:    T[t, @R, phi] 
T[t, @center,] (assume) = 21.6 (assume it be some small varying function)
        
        
 """       
        
        
        
        
        
        
        
        
#%% visualization with mayavi
#surf = mlab.mesh(X.T, Y.T, 10*T[999], colormap='RdYlBu')
#mlab.show()

