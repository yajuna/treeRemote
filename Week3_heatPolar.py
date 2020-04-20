#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 27 20:03:51 2020

Week 3 - Selina Teng
Goal 1: figure out runtime warnings on line 66 (due to variable k3)
    Blind testing -> with Nr and N_phi held constant at 50, we encounter the runtime warning
    once N_steps exceeds 133.
    Resolved matrix by changing the type of matrix T to int64. Warning occurred because
    we were reaching integer values of magnitude larger than 2,147,483,648, which is 
    unacceptable in the int32 datatype. 
    Source: https://stackoverflow.com/questions/7559595/python-runtimewarning-overflow-encountered-in-long-scalars
Goal 2: add source terms to lines 63-66
    Copied over code from sourceVars.py (In the future, should we import that part from another file?)
    Line 129 - added the source terms to T 
Still need to do:
    Express changing variable i in the solar radiation function (i is the angle of incidence)
    Test for correctness
Miscellaneous changes/questions
    Commented out imported libraries that are unused for now
    Is it stylistically ok that I added functions in this script?


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

#%% Clear existing variables and console
try:
    from IPython import get_ipython
    get_ipython().magic('clear')
    get_ipython().magic('reset -f')
except:
    pass

#%% Imports
import numpy as np
import math
from scipy.special import jv, jn_zeros
#from scipy.integrate import odeint
#import matplotlib.pyplot as plt
#from matplotlib import cm
#from mpl_toolkits.mplot3d import Axes3D
#from mayavi import mlab

#%% Parameters
Nr = 50
N_phi = 50
N_steps = 1000
print("Number of steps: " + str(N_steps))
radius = 5.  
c = 0.5 # k/(rho*c) in equation (2) of Potter and Anredsen paper
dphi = 2*np.pi/N_phi  
dr = 5./Nr  
dt = 0.005   
if dt< (dr*dphi)**2/2.0: # CFL condition. 
    # The maximum value for dt is dr*dphi/(2*c)
    dt = dr*dphi/(4*c)

#%% Initial conditions

# Heat transfer
r = np.linspace(0, radius, Nr)  
phi = np.linspace(0, 2*np.pi, N_phi)  
R, phi = np.meshgrid(r, phi) 
X = R*np.cos(phi) 
Y = R*np.sin(phi)
kth_zero = jn_zeros(1, 1) 
Z = np.cos(phi) * jv(1, kth_zero*R/radius)
T = np.zeros((N_steps, Nr, N_phi),dtype='int64') #Changing type to int64 resolves warning
T[0, :, :] = Z.T 
#Conduction
r_max = 1           #Max. radius of tree
ro = 1              #Wood density (kg/m^3)
c = 1               #Specific heat (J/kg*K)
T_in = 1            #Temperature of tree interior (K)
k = 1               #Thermal conductivity (W/(m*K))
#Convective heat loss 
T_sfc = 1           #Temperature of tree surface (K)
T_air = 1           #Temperature of surrounding air (K)
L = 1               #Vertical height (m)
R = 1               #Tree radius (m)
theta1 = np.pi      #Difference between wind direction and the aspect of the surface point (radians)
u = 1               #Windspeed (m/s)   
#Solar radiative heating
tau = 0.7           #Atmospheric transmissivity (unitless) - at midlatitude 
                    #& low elevation, approx 0.76-0.81                   
eta = 0.8           #Atmospheric absorption parameter (unitless) - at 
                    #midlatitude & low elevation, approx 0.80-0.84                    
Z = 2 * np.pi       #Solar zenith angle (radians)
i = 2 * np.pi       #Angle of incidence, the angle between the direction of 
                    #sunlight and the local normal to the tree's surface                   
alpha = 0           #Albedo of surface
#Long wave radiation
    #This function also uses T_sfc and T_air (see #2)

#%% Main
def main():
    result = stepping()

#%% Stepping
def stepping():
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
                + dt/dr*sourceVars()    
            T[t, i, -1] = T[t, i, 0]  # Update the values for phi=2*pi, BC in phi
    return T;

#%% Source variables 
def sourceVars():   
    H = convectiveHeatLoss(T_sfc, T_air, L, R, theta1, u)
    S = solarRadiationHeating(tau, eta, Z, i, alpha)
    IR = longWaveRadiation(T_sfc, T_air)
    tot = H + S + IR
    return tot

def convectiveHeatLoss(T_sfc, T_air, L, R, theta1, u):
    if (theta1 < (2 * np.pi)):   #Theta is the lesser of 2pi radians or theta1
        theta = theta1 
    else: 
        theta = 2 * np.pi
    h_free = 18.293 * abs(T_sfc - T_air) / (L * T_air**3)**0.25 * (T_air + 97.77) / math.sqrt(179.02 + T_air)
    h_forced = 3.458 * (T_air + T_sfc - 0.74)**0.49 * (u / (R * T_air))**0.5 * (1 - (theta / 90)**3)
    h = h_free + h_forced
    H = h * (T_sfc - T_air) 
    return H

def solarRadiationHeating(tau, eta, Z, i, alpha):
    #Total solar radiation incident at a specific surface point is the sum of 
    #direct and diffuse solar radiation
    S_0 = 1368          #Solar constant (W/m^2)
    #Should do loop for multiple i's
    if (i > 90):
        S_dir = 0
    else:
        S_dir = S_0 * tau**(1 / math.cos(Z)) * math.cos(i)    
    S_dif = S_0 * math.cos(Z) / 3 * (1 + math.cos(Z)) * (eta - (tau)**(1 / math.cos(Z)))
    S = (S_dir + S_dif)*(1 - alpha)
    return S

def longWaveRadiation(T_sfc, T_air):
    sigma = 5.67 * 10**(-8)     #Stefan Boltzmann constant (W/(m^2*K^4))
    IR_out = sigma * T_sfc**4
    IR_in = sigma * T_air**4
    IR = IR_in - IR_out
    return IR

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

if __name__ == '__main__':
    main()