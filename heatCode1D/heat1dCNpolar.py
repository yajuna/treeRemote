#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 10:04:31 2020
Rewrite heat1D_CNpolar_sourceUpdate3.py as a function. 

1D code for heat equation in polar, with Crank-Nicolson
Ref: LeVeque Chapter 9, pg 183
Handling of source term: a primer on PDEs, Salsa et al. pg. 99
u_t = a*(1/r*u_r + u_rr)
a = rho*c/k
u(x,0)=eta(x) initial condition
d u(0,t)/d r=0 Neumann condition in the center
u(R,t)=g1(t) Dirichlet condition at the bark of tree
ri = i*dr, tn = n*dt, dr = Delta r, dt = Delta t
gs(t) is a source term only occuring at the boundary
CN:
stability to test; standard heat equation CN stable for any dt>0. Take dt = O(dr)

Update issues: 
    1. plot temperature at center: increase number of grid points
    observation: when chaning gridPoints from 50 to 500, figures do not coincide (eg. 380/500
    and 38/50 are not the same).
    2. variable conductivity
    
to test parameters:
run heat1dCNpolar

for j in range(1,50,5):
    config['at_point'] = j
    c = temp(config)    
@author: yajun
"""
import numpy as np # version 1.16.4
from scipy import sparse # version 1.3.0
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D # matplotlib version 3.1.0

# import source terms and bdry conditions
from sourceS import * # source term at bdry
from Temp_dataVec import * # boundary at tree bark
#%%
config = dict()
config['gridPoints'] = 50
config['timeSteps'] = 1000
config['thermalConductivity'] = 0.12
config['heatCapacity_rhoc'] = 1.7
config['at_point'] = 38
config['time'] = np.linspace(0, 1000, 50, endpoint = False)
#%%
def temp(config):
    m = config['gridPoints']
    n = config['timeSteps']
    k = config['thermalConductivity']
    rhoc = config['heatCapacity_rhoc']
    at_point = config['at_point'] 
    
    time = config['time']
    timeInt = time.astype(int)
    timelist = timeInt.tolist()
    
    r = np.linspace(0, 1, m, endpoint=False)
    t = np.linspace(0, 12, n, endpoint=False)
    
    dr = r[1] - r[0]
    dt = t[1] - t[0]
    a = rhoc / k
    
    # parameters defined as in paper
    beta = dt / (a * dr ** 2)

    alpha = np.ones(m - 1)
    for j in range(alpha.size):
        alpha[j] = dt / (4 * a * r[j + 1] * dr) - 0.5 * beta
# use following for neumann condition at outer boundary
# alpha[-1] = -beta

    gamma = np.ones(m - 1)
    for j in range(1, gamma.size):
        gamma[j] = dt / (4 * a * r[j] * dr) + 0.5 * beta
    gamma[0] = beta  # define inner neumann bdry condition 1st row

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
        rhs[-1] = rhs[-1] + (dt / (4 * a * r[-1] * dr) + 0.5 * beta) * (g1(i) + g1(i + 1)) + dt/(2*a*dr*k)*(gs(i) + gs(i + 1))

    # if neumann condition for g1(t)
    # rhs[-1] = rhs[-1] + 2* dr *(dt / (4 * a * r[-1] * dr) + 0.5 * beta) * (g1(i) + g1(i + 1)) + dt/(2*a*dr*k)*(gs(i) + gs(i + 1))
        U1 = np.linalg.solve(tridiag, rhs)  # sparse.linalg.lsqr(tridiag, rhs)
        U0 = U1
        soln.append(U0)

# %% print solutions
#    print(soln)
    soln_plot = np.asarray(soln)
#    print(soln_plot.shape)
    print("max and min of soln at ", at_point, " = ", np.max(soln_plot[:, at_point]), np.min(soln_plot[:, at_point]))
    
#%% visualize
#    plt.plot(t, soln_plot[:, at_point], '.r-')
#    plt.title('Temperature with combined source term at grid point r=%i ' % at_point)
#    plt.axis([0,12,280,305])
#    plt.xlabel('Time since 7:00am (hrs)')
#    plt.ylabel('Temperature Distribution (K)')
#    plt.savefig('/home/yajun/Documents/treePower/figs/' + 'Stemp' + str(at_point) + '.eps', format='eps', dpi=300,bbox_inches='tight')
#    plt.show()
#%% 3D visualization, plot a wireframe.
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    X, Y = np.meshgrid(r, timelist)
    solnFig = soln_plot[timelist, :]
#    zs = np.array([solnFig for r, timelist in zip(np.ravel(X), np.ravel(Y))])
    Z = solnFig.reshape(X.shape)

    ax.plot_wireframe(X, Y, Z, rstride=1, cstride=1)

    ax.set_xlabel('Radius grid')
    ax.set_ylabel('Time step')
    ax.set_zlabel('Temperature (K)')

    ax.set_title('Temperature distribution')
#    plt.savefig('/home/yajun/Documents/treePower/figs/' + 'Stemp' + str(at_point) + '.eps', format='eps', dpi=300,bbox_inches='tight')
    plt.show()

    return soln_plot
######################################################################
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















