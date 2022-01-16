#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 12 09:43:16 2021

This code is a temporary version to test parameters and coding variations. See heatMain.py
for parameter explanations. 
****DO NOT ADD, PUSH, COMMIT****

to run:
run testHeat
c = temp(config)  

to test parameter:
run testHeat
for j in range(1,500,10):
    config['heatCapacity_rhoc'] = j
    print("heatCapacity_rhoc is", j)
    c = temp(config)

10/4: omit TNTvec results; use experimental data to generate temperatures; changes at
      boundary condition g1(t)
@author: yajun
"""

import numpy as np
from scipy import sparse
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D # matplotlib version 3.1.0
# import source terms and bdry conditions
import source as stree # source term at bdry
import TNTvec as tntv # boundary at tree bark

#%%
config = dict()
config['gridPoints'] = 160
config['timeSteps'] = 1000
#mu,sigma = 1, 0.001 # mean and standard deviation
config['thermalConductivity'] = 0.12*np.ones(6) # 0.12*np.random.normal(mu, sigma, 6)
config['heatCapacity_rhoc'] = 510*1380
config['at_point'] = 150
config['point_pair'] = [1,-2]
config['time'] = np.linspace(0, 1000, 50, endpoint = False)
config['output'] = 'temp'
config['visualization dimension'] = 1 # 2, 3

#%%
def temp(config):
    m = config['gridPoints']
    n = config['timeSteps']
    k = config['thermalConductivity']
    k = np.asarray(k)
    rhoc = config['heatCapacity_rhoc']
    at_point = config['at_point'] 
    point_pair = config['point_pair']
    vis = config['visualization dimension']
    output = config['output']
    time = config['time']
    
    timeInt = time.astype(int)
    timelist = timeInt.tolist()
    
    r = np.linspace(0, 0.36 * 0.5, m, endpoint=False)    
    t = np.linspace(0, 24, n, endpoint=False)
    
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
    print("number of time steps, or length of solution list", len(soln))
    soln_plot = np.asarray(soln)
    print("max and min of soln_plot", soln_plot.max(), soln_plot.min())
    soln_diff = soln_plot[:,point_pair[0]]-soln_plot[:,point_pair[1]]
    print("shape of soln_diff, marks temp diff between two points all day", soln_diff.shape)
    if np.max(np.abs(soln_diff))>1e-3:
        print("max difference throughout the day: core and bark at", point_pair, "=", np.max(np.abs(soln_diff)))
##%% visualize
#    plt.plot(t, soln_plot[:, at_point], '.r-')
#    plt.title('Temperature with combined source term at grid point r=%i ' % at_point)
#    plt.axis([0,12,280,305])
#    plt.xlabel('Time since 7:00am (hrs)')
#    plt.ylabel('Temperature Distribution (K)')
#    plt.savefig('/home/yajun/Documents/treePower/figs/' + 'StempK' + str(at_point) + '.eps', format='eps', dpi=300,bbox_inches='tight')
#    plt.show()

#    return np.max(soln_plot[-1, :]), np.min(soln_plot[-1, :])
    return np.max(np.abs(soln_plot[:,0]-soln_plot[:,-1]))
#%% 3D visualization, plot a wireframe.
#    fig = plt.figure()
#    ax = fig.add_subplot(111, projection='3d')
#
#    X, Y = np.meshgrid(r, timelist)
#    solnFig = soln_plot[timelist, :]
##    zs = np.array([solnFig for r, timelist in zip(np.ravel(X), np.ravel(Y))])
#    Z = solnFig.reshape(X.shape)
#
#    ax.plot_wireframe(X, Y, Z, rstride=1, cstride=1)
#
#    ax.set_xlabel('Radius grid')
#    ax.set_ylabel('Time step')
#    ax.set_zlabel('Temperature (K)')
#
#    ax.set_title('Temperature distribution')
##    plt.savefig('/home/yajun/Documents/treePower/figs/' + 'Stemp' + str(at_point) + '.eps', format='eps', dpi=300,bbox_inches='tight')
#    plt.show()
#
#    return soln_plot
######################################################################
# %% define bdry and initial conditions. 

# initial condition being const. temp at 7 am
def eta(m):
    return 300.40 * np.ones(m)# ((56 - 32) * 5/ 9 + 273.15) * np.ones(m)

## tree center bdry condition is homogeneous Neumann condition. In matrix
    
# tree bark with Dirichlet condition for temperature
def g1(t):
    return tntv.bdryArray[t]
    return stree.barkTemp[t]
# source term at tree bark
def gs(t):
    
    return stree.sourceTerm[t]