# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

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
    vectorize to increase efficiency

@author: yajun
"""
"""
Questions about workflow:
    Is it correct to copy-paste the file from Git, edit it in my choice of IDE, then copy-paste my edited file back into the Github browser window? 
    (Is there a better option than copy-paste? I know I can push files from commandline, but I am confused about how.)
    Is Spyder a good choice of IDE? When do we use an IDE like Spyder as opposed to something like Jupyter Notebook?
    Where does Spyder save my file when I hit save? Does temporary mean it goes away when I close out or are they somewhere eating up storage space?
    When I clone the Git repository to my computer using the clone command from commandline, is that a live updating version?
    How to avoid having a lot of separate versions of the script saved locally on my computer that are hard to keep track of?
Questions about the script:
    Here is the error message I am getting when I attempt to run. I think this is the problem we were talking about before, where
    the first warning has to do with something probably being the wrong order of magnitude, and the second warning is because 
    something is a negative value that shouldn't be. Just to make sure, this is the expected result?
        /Users/selinateng/.spyder-py3/temp.py:73: RuntimeWarning: overflow encountered in double_scalars
          + k3*(T[t-1, i, j+1] - 2*T[t-1, i, j] + T[t-1, i, j-1])
        /Users/selinateng/.spyder-py3/temp.py:73: RuntimeWarning: invalid value encountered in double_scalars
          + k3*(T[t-1, i, j+1] - 2*T[t-1, i, j] + T[t-1, i, j-1])
      
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
N_steps = 1000  
radius = 5.  
c = 0.5 # k/(rho*c) in equation (2) of Potter and Anredsen paper
dphi = 2*np.pi/N_phi  
dr = 5./Nr  
dt = 0.005   
if dt< (dr*dphi)**2/2.0: # CFL condition. 
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

#%% Stepping
k1 = c*dt/dr**2
for t in range(1, N_steps):
    for i in range(0, Nr-1):
        for j in range(0, N_phi-1):
            ri = max(r[i], 0.5*dr)  # To avoid the singularity at r=0
            k2 = c*dt/(2*ri*dr)
            k3 = c*dt/(dphi*ri)**2
            T[t, i, j] = + T[t-1, i, j] \
            + k1*(T[t-1, i+1, j] - 2*T[t-1, i, j] + T[t-1, i-1, j])\
            + k2*(T[t-1, i+1, j] - T[t-1, i-1, j])\
            + k3*(T[t-1, i, j+1] - 2*T[t-1, i, j] + T[t-1, i, j-1])

        T[t, i, -1] = T[t, i, 0]  # Update the values for phi=2*pi

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

