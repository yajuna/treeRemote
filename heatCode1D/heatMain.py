#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 27 13:15:36 2021

this code incorporates the main code heat1dK.py, and replaces current conditions 
by data generated by code TNTVec.py (h for heat convective coefficient, 
and Tbfinal as the temperature at the tree bark g1 at each time generated by measured
tree center data), weather data in TNTweather.csv. Experimental data on Feb 16

To validate computation, use data from Feb x, x >=17. 

Parameters replaced by Protasio's data on mango trees.

correct value: for soft wood, we should use 
1: density 510 kg/m^3
2: specific heat 1380 J/(kg * K)

diameter: 36 cm = 0.36 m; grid size = .36 * 0.5 / numberGridPoints

function can output temperature of the tree trunk at chosen grid, or temperature
differences between two grid points.

to run:
run heatMain
c = temp(config)  

to test parameter for bark grid (centers are roughly the same)
run heatMain
for j in range(0, -30, -1):
    config['point_pair'] = [5,j]
    c = temp(config)

@author: yajuna
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
config['at_point'] = 38
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

#%% print solutions
#    print(soln)
    soln_plot = np.asarray(soln)
#    print("max and min of soln at ", at_point, " = ", np.max(soln_plot[:, at_point]), np.min(soln_plot[:, at_point]))
    
    soln_diff = soln_plot[:,point_pair[0]]-soln_plot[:,point_pair[1]]
    if np.max(np.abs(soln_diff))>1e-3:
        print("max difference throughout the day: core and bark at", point_pair, "=", np.max(np.abs(soln_diff)))

#%% 2D visualize
    if vis == 2:
        plt.plot(t, soln_diff, '.r-')
        message = f"Temperature Difference at grid points {point_pair}"
        plt.title(message)
        plt.axis([0,24,-1,1])
        plt.xlabel('Time (hrs)')
        plt.ylabel('Temperature Difference (K)')
#        plt.savefig('/home/yajun/Documents/treePower/figs/' + 'StempK' + str(at_point) + '.eps', format='eps', dpi=300,bbox_inches='tight')
        plt.show()

#%% 3D visualization, plot a wireframe.
    elif vis == 3:    
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        X, Y = np.meshgrid(r, timelist)
        solnFig = soln_plot[timelist, :]
#    zs = np.array([solnFig for r, timelist in zip(np.ravel(X), np.ravel(Y))])
        Z = solnFig.reshape(X.shape)
        ax.plot_wireframe(X, Y, Z, rstride=1, cstride=1)
        ax.set_xlabel('Radius grid')
        ax.set_ylabel('Time step')
        ax.set_zlabel('$\Delta T$ (K)')
        message = f"Temperature Distribution (K)" # if include variable, message = f" variable is {variable}"
        ax.set_title(message)
#    plt.savefig('/home/yajun/Documents/treePower/figs/' + 'Stemp' + str(at_point) + '.eps', format='eps', dpi=300,bbox_inches='tight')
        plt.show()
        
#%% no visualization for quick testing        
    else: 
        None  

#%% output all solution with 'tmep' or only temp difference between two grid points
    if output == 'temp':
        return soln_plot
    else:
        return soln_diff
######################################################################
# %% define bdry and initial conditions. 
# initial condition being const. approx. temp at 12 am 
def eta(m):
    return 300.40 * np.ones(m)
## tree center bdry condition is homogeneous Neumann condition. In matrix

## tree bark with Dirichlet condition for temperature
def g1(t):
    return tntv.bdryArray[t]
# source term at tree bark
def gs(t):
    return stree.sourceTerm[t]
