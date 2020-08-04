#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tuesday 07 14 14:55:12 2020
1D code for heat equation in polar, with Crank-Nicolson
Ref: LeVeque Chapter 9, pg 183
Handling of source term: a primer on PDEs, Salsa et al.

u_t = a*(1/r*u_r + u_rr)
a = rho*c/k
u(x,0)=eta(x) initial condition
d u(0,t)/d r=0 Neumann condition in the center
u(R,t)=g1(t) Dirichlet condition at the bark of tree
ri = i*dr, tn = n*dt, dr = Delta r, dt = Delta t
gs(t) is a source term only occuring at the boundary

CN:
stability to test; standard heat equation CN stable for any dt>0. Take dt = O(dr)
@author: yajun
"""

import numpy as np
from scipy import sparse
from matplotlib import pyplot as plt
from combinedHeatsource_updated import *
from Temp_data import *

# number of grid points, space grid points r.
m = 50
r = np.linspace(0, 1, m, endpoint=False)

# number of time steps
n = 1000
t0 = 0
t = np.linspace(t0, 12, n, endpoint=False)

# define dr and dt

dr = r[1] - 0  # grid size
dt = t[1] - t0  # time step

###
k = 0.12
a = 1.7/k

# parameters defined as in paper
beta = dt / (a * dr ** 2)

alpha = np.ones(m - 1)
for j in range(alpha.size):
    alpha[j] = dt / (4 * a * r[j + 1] * dt) - 0.5 * beta


gamma = np.ones(m - 1)
for j in range(1, gamma.size):
    gamma[j] = dt / (4 * a * r[j + 1] * dt) + 0.5 * beta
#    print(r[j+1])
gamma[0] = beta  # define neumann bdry condition 1st row


######################################################################
# %% define bdry and initial conditions. Just examples, to modify

# initial condition being a small cos func
def eta(r):
    return tempArray(r)


## tree center bdry condition is homogeneous Neumann condition. In matrix

## tree bark with Dirichlet condition for temperature

def g1(t):
    # gaussian(x, mu, sig)
    mu = 7 # start meansuring at 7 am, peak at 14, unit kelvin
    sig = 1
    return np.exp(-np.power(t - mu, 2.) / (2 * np.power(sig, 2.))) + 293

# source term at tree bark
def gs(t):
    # gaussian(x, mu, sig)

    return sourceTermsSvalue(t)


#######################################################################

# %% define matrices in time stepping. Diagonal/super diagonal depend on ri

# diag = np.ones(m) + beta

tridiag = sparse.diags([alpha, np.ones(m) + beta, -gamma], [-1, 0, 1], shape=(m, m)).toarray()

# %%

soln = []

# IC
U0 = eta(r)

# %%
#  Solving Ax1 = Bx0 + b, this is B. No flux results in different last row from before

B = sparse.diags([-alpha, np.ones(m) - beta, gamma], [-1, 0, 1], shape=(m, m)).toarray()



soln.append(U0)

# %% main time stepping: compute rhs = Bx0 + b, then solve for x1 with Ax1 = rhs

for i in range(n - 1):
    rhs = B.dot(U0)
    # dirichlet bdry condition with g1, source at bdry with gs
    # rhs[-1] = rhs[-1] + (dt / (4 * a * r[-1] * dr) + 0.5 * beta) * (g1(i) + g1(i + 1))
    rhs[-1] = rhs[-1] + (dt / (4 * a * r[-1] * dr) + 0.5 * beta) * (g1(i) + g1(i + 1)) + dt/(2*a*dr*k)*(gs(i) + gs(i + 1))
    U1 = np.linalg.solve(tridiag, rhs)  # sparse.linalg.lsqr(tridiag, rhs)
    U0 = U1
    soln.append(U0)
#    print(np.max(soln[i, :]))
# %%
print(soln)
soln_plot = np.asarray(soln)

print(soln_plot.shape)

print("max and min of soln at final step = ", np.max(soln_plot[-1, :]), np.min(soln_plot[-1, :]))

grid_point = 29
plt.plot(t, soln_plot[:, grid_point], '.r-')
plt.title('Temperature with combined source term at grid point r=%i ' % grid_point)
plt.xlabel('Time since 7:00am (hrs)')
plt.ylabel('Temperature Distribution')
plt.show()