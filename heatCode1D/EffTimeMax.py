#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  3 20:07:33 2020

An efficient code to spit out time that gives the max temp differences

combines heat1dK.py and tempMax.py

Find the index of N max elements. For Numpy version higher than 1.8 (currently 1.18.4)
https://stackoverflow.com/questions/6910641/how-do-i-get-indices-of-n-maximum-values-in-a-numpy-array

Two peaks in temp curve might not be captured by two max values (local max and min)
@author: yajun
"""
import numpy as np
from scipy import sparse
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D # matplotlib version 3.1.0
# import source terms and bdry conditions
from sourceS import * # source term at bdry
from Temp_dataVec import * # boundary at tree bark

config = dict()

config['timeSteps'] = 1000 # n
mu,sigma = 1, 0.001 # mean and standard deviation
config['thermalConductivity'] = 0.12*np.ones(12) # k
config['heatCapacity_rhoc'] = 1.7  # rhoc
config['time'] = np.linspace(0, 1000, 50, endpoint = False)

m = 75 # number of grid points

# %% define bdry and initial conditions. 

# initial condition being const. temp at 7 am
def eta(m):
    return ((56 - 32) * 5/ 9 + 273.15) * np.ones(m)

## tree center bdry condition is homogeneous Neumann condition. In matrix

## tree bark with Dirichlet condition for temperature
def g1(t):
    return tTemp[t]

# source term at tree bark
def gs(t):
    
    return sourceTermsSvalue(t)

#%% main function

def tempTime(m):
    
    n = 1000
    k = 0.12*np.ones(12)
    k = np.asarray(k)
    rhoc = 1.7
    
    r = np.linspace(0, 1, m, endpoint=False)
    t = np.linspace(0, 12, n, endpoint=False)
    
    dr = r[1] - r[0]
    dt = t[1] - t[0]
    
    # interpolate thermal conductivity. Assuming data measured uniformly
    kinterp = np.interp(r, np.linspace(0,k.size - 1, k.size), k)
    a = rhoc / kinterp
    
    # parameters defined as in paper
    beta = np.ones(m)
    for j in range(beta.size):
        beta[j] = dt / (a[j] * dr ** 2)

    alpha = np.ones(m - 1)
    for j in range(alpha.size - 1):
        alpha[j] = dt / (4 * a[j + 1] * r[j + 1] * dr) - 0.5 * beta[j + 1] + dt / (4 * dr ** 2) * (1. / a[j + 2] - 1. / a[j])
    # last alpha value separately defined, with different approximation for dk/dr, reflected on index of a
    alpha[-1] = dt / (4 * a[-1] * r[-1] * dr) - 0.5 * beta[-1] + dt / (2 * dr ** 2) * (1. / a[-1] - 1. / a[-2])  
# use following for neumann condition at outer boundary
# alpha[-1] = -beta

    gamma = np.ones(m - 1)
    for j in range(1, gamma.size):
        gamma[j] = dt / (4 * a[j] * r[j] * dr) + 0.5 * beta[j] + dt / (4 * dr ** 2) * (1. / a[j + 1] - 1. / a[j - 1])
    gamma[0] = beta[0]  # define inner neumann bdry condition 1st row

    tridiag = sparse.diags([alpha, np.ones(m) + beta, -gamma], [-1, 0, 1], shape=(m, m)).toarray()
#%%
    soln = []
    U0 = eta(m) 
#  Solving Ax1 = Bx0 + b, this is B. 
    B = sparse.diags([-alpha, np.ones(m) - beta, gamma], [-1, 0, 1], shape=(m, m)).toarray()
    soln.append(U0)

# %% main time stepping: compute rhs = Bx0 + b, then solve for x1 with Ax1 = rhs
    for i in range(n - 1):
        rhs = B.dot(U0)
    # dirichlet bdry condition with g1, source at bdry with gs; average gs at time i and i + 1
        rhs[-1] = rhs[-1] + gamma[-1] * (g1(i) + g1(i + 1)) + dt/(2 * a[-1] * dr * k[-1]) * (gs(i) + gs(i + 1))
 
    # if neumann condition for g1(t)
    # rhs[-1] = rhs[-1] + 2* dr *(dt / (4 * a * r[-1] * dr) + 0.5 * beta) * (g1(i) + g1(i + 1)) + dt/(2*a*dr*k)*(gs(i) + gs(i + 1))
        U1 = np.linalg.solve(tridiag, rhs)  # sparse.linalg.lsqr(tridiag, rhs)
        U0 = U1
        soln.append(U0)

# %% print solutions
    soln_plot = np.asarray(soln)
    
    tempDiff = np.zeros(n)
    for j in range(n):
        tempDiff[j] = np.abs(soln_plot[j,0] - soln_plot[j,-1])

# to find the time for max temp    
#    MaxTempIndex = np.argmax(tempDiff) 
#    print("Max temperature difference occurs at", MaxTempIndex, "time step", "with grid point number ", m, ", the difference is", tempDiff[MaxTempIndex])   

# to find the time for two max temp replace 2 by N to find N max values
    ind = np.argpartition(tempDiff, -2)[-2:]
    print("Max temperature difference occurs at", ind, "time step", "with grid point number ", m, ", the difference is", tempDiff[ind])   
         
    
    return 

#%% for loop to get time for max temp difference
    
for m in range(50, 500, 50):
    tempTime(m)
    

















