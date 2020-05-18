#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  6 16:43:14 2020

ref: https://scicomp.stackexchange.com/questions/20279/how-can-i-solve-wave-equation-for-circular-membrane-in-polar-coordinates

visualization ref with mayavi: https://docs.enthought.com/mayavi/mayavi/installation.html

This code solves equation (10) in Potter and Andresen


5/6 Update - Selina - Added in source terms and corrected some issues with counting.

QUESTION - It looks like the effect of direct solar radiation is very small. Is it worth
calculating it for every different j?

Here are the source terms when t = 0, i = 0, j = 0:Nphi-1
There is not much change...
    
107.59818210302105
107.59816991258005
107.59816219198879
107.59815266528275
107.59814413188482
107.59813347628265
107.59812683909415
107.59811826043644
107.59810688224756
107.59809410436733
107.59808548041946
107.59807626943598
107.5980676905669
107.59805658302504
107.59804633348497
107.5980345033526
107.59802335060242
107.5980137329442
107.5980017672836
107.59799431687024
107.59797991280155
107.5979684886933
107.59795697433499
107.5979438794254
107.59793742225385
107.59792739778064
107.59791854734705
107.59791159336403
107.59790346538597
107.59789475032362
107.5978876156279
107.59787871986605
107.59787113353202
107.59786137968607
107.59785311603359
107.59784272987743
107.5978345113087
107.59782222849056
107.59781053267484
107.59779928826643
107.59779143080482
107.59778416027036
107.59777558018493
107.59776618717089
107.59775887138804
107.59775218783159
107.59774274954931
107.59773385316213
107.59772500189787
    


@author: yajun
"""
import numpy as np
from scipy.special import jv, jn_zeros
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from sourceTerms import sourceTerms
#from mayavi import mlab 

#%% Parameters; play with parameters, and see if warnings would go away
Nr = 50  
N_phi = 50  
N_steps = 1000  
radius = 5.  
c = 0.5 # = k/(rho*c) in equation (2) of Potter and Anredsen paper
""" Changed line 38 from /(N_phi) to /(Nphi - 1) so that max phi = 2pi """
dphi = 2*np.pi/(N_phi - 1)
dr = 5./Nr  
dt = 0.005   
if dt< (dr*dphi)**2/2.0: # CFL condition. 
    # The maximum value for dt is dr*dphi/(2*c)
    dt = dr*dphi/(4*c)


## for Nr in range(2,500):
#%% Initial conditions; waiting for data
r = np.linspace(0, radius, Nr)  
phi = np.linspace(0, 2*np.pi, N_phi)  
R, phi = np.meshgrid(r, phi) 
X = R*np.cos(phi) 
Y = R*np.sin(phi)
kth_zero = jn_zeros(1, 1) 
Z = np.cos(phi) * jv(1, kth_zero*R/radius)
T = np.zeros((N_steps, Nr, N_phi),dtype = 'int64')
T[0, :, :] = Z.T 
#%% visualize initial condition in 3D

fig = plt.figure()
ax = fig.gca(projection='3d')

# Plot the surface.
surf = ax.plot_surface(R, phi, Z, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)

fig = plt.figure()
ax = fig.gca(projection='3d')

surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)

    
#%% Stepping
k1 = c*dt/dr**2
for t in range(1, N_steps):
    """ Increased the max limit for i and j by 1 
        In Python, the max is not included in the range, so before, it only 
        completed 49 iterations instead of 50 like we wanted
    """
    for i in range(0, Nr):
        #disp(i)
        for j in range(0, N_phi):
            #disp(j)
            # Calculate phi in radians (phi_0 = 0 rad)
            phi = j*dphi 
            disp(phi)
            #disp(sourceTerms(phi))
            ri = max(r[i], 0.5*dr)  # To avoid the singularity at r=0
            k2 = c*dt/(2*ri*dr)
            k3 = c*dt/(dphi*ri)**2
            """ k4 here """
            k4 = c*dt/dr
            T[t, i, j] = + T[t-1, i, j] \
            + k1*(T[t-1, i+1, j] - 2*T[t-1, i, j] + T[t-1, i-1, j])\
            + k2*(T[t-1, i+1, j] - T[t-1, i-1, j])\
            + k3*(T[t-1, i, j+1] - 2*T[t-1, i, j] + T[t-1, i, j-1])#\
#            + source terms 
            """ Added in source terms"""
            - k4*sourceTerms(phi)
        
        T[t, i, -1] = T[t, i, 0]  # Update the values for phi=2*pi, BC in phi
#        T[t,0,i]?? T[t,R,i]??
        
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

