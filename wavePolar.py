#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 27 14:54:53 2020

ref: https://scicomp.stackexchange.com/questions/20279/how-can-i-solve-wave-equation-for-circular-membrane-in-polar-coordinates

visualization ref with mayavi: https://docs.enthought.com/mayavi/mayavi/installation.html

@author: yajun
"""

import numpy as np
from scipy.special import jv, jn_zeros
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
#from mayavi import mlab


#%% Parameters
Nr = 50
N_phi = 50
N_steps = 200
radius = 5.
c = 0.5
dphi = 2*np.pi/N_phi
dr = 5./Nr
dt = 0.005  
if dt< dr*dphi/2/c:
    # The maximum value for dt is dr*dphi/(2*c)
    dt = dr*dphi/(4*c)

#%% Initial conditions
r = np.linspace(0, radius, Nr)
phi = np.linspace(0, 2*np.pi, N_phi)
R, phi = np.meshgrid(r, phi)
X = R*np.cos(phi)
Y = R*np.sin(phi)
kth_zero = jn_zeros(1, 1)
Z = np.cos(phi) * jv(1, kth_zero*R/radius)
T = np.zeros((N_steps, Nr, N_phi))
T[0, :, :] = Z.T
T[1, :, :] = Z.T

#%% Stepping
k1 = c*dt**2/dr**2
for t in range(2, N_steps):
    for i in range(0, Nr-1):
        for j in range(0, N_phi-1):
            ri = max(r[i], 0.5*dr)  # To avoid the singularity at r=0
            k2 = c*dt**2/(2*ri*dr)
            k3 = c*dt**2/(dphi*ri)**2
    
            T[t, i, j] = 2*T[t-1, i, j] - T[t-2, i, j] \
            + k1*(T[t-1, i+1, j] - 2*T[t-1, i, j] + T[t-1, i-1, j])\
            + k2*(T[t-1, i+1, j] - T[t-1, i-1, j])\
            + k3*(T[t-1, i, j+1] - 2*T[t-1, i, j] + T[t-1, i, j-1])

        T[t, i, -1] = T[t, i, 0]  # Update the values for phi=2*pi
#        print(k2,k3,i)
    print(np.max(abs(T[t,:,:])),t)
#%% visualization with matplolib BUGGY, NEED TO FIX
#fig = plt.figure()
#ax = fig.gca(projection='3d')
#
#u,v = np.meshgrid(X.T, Y.T)        
#surf = ax.plot_surface(u, v, 10*T[999], cmap = cm.coolwarm, linewidth = 0, antialiased = False)
#fig.colorbar(surf, shrink = 0.5, aspect = 5)
#plt.show()       
        
        
        
        
        
        
        
        
        
        
        
        
        
#%% visualization with mayavi
#surf = mlab.mesh(X.T, Y.T, 10*T[999], colormap='RdYlBu')
#mlab.show()