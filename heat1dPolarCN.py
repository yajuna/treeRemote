#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tuesday 07 14 14:55:12 2020

1D code for heat equation, with Crank-Nicolson

Ref: LeVeque Chapter 9, pg 183

u_t = kappa*(1/r*u_r + u_rr) (with kappa = 1 WLOG)

u(x,0)=eta(x)
u(0,t)=g0(t)
u(1,t)=g1(t)

xi = ih, tn = nk, h = Delta x, k = Delta t

CN:
(U^n+1_i - U^n_i)/k =k/2h*1/r_i(U^n_i+1 - U^n_i + U^n+1_i+1 - U^n+1_i) + 1/2h^2(U^n_i-1 - 2U^n_i + U^n_i+1 + U^n+1_i-1 - 2U^n+1_i + U^n+1_i+1) + source term  

rewrite for main time stepping:
-k/2h^2 * T^n+1_i-1 + (1 + k/2h * 1/ri + k/h^2) * T^n+1_i + (-k/h^2 * 1/ri - k/2h^2) * T^n+1_i+1
=
k/2h^2 * T^n_i-1 + (1 - k/2h * 1/ri - k/h^2) * T^n_i + (k/h^2 * 1/ri + k/2h^2) * T^n_i+1

bdry condition: 

rhs[0] = k/2h^2 (g0(tn) + g0(tn + k)) + (1-k/2h * 1/r1 - k/h^2)* T^n_1 + k/2h * 1/r1 + k/2h^2) * T^n_2

rhs[-1] = k/2h^2* T^n_m-1 + (1-k/2h * 1/rm - k/h^2) * T^n_m + (k/2h * 1/rm + k/2h^2) * (g1(tn) + g1(tn + k))
 

stability to test; standard heat equation CN stable for any k>0. Take k = O(h)
@author: yajun
"""

import numpy as np
from scipy import sparse
from matplotlib import pyplot as plt

# number of grid points, space grid points x
m = 50
x = np.linspace(0, 1, m)

# number of time steps
n = 1000
t0 = 0
t = np.linspace(t0, 8, n)

# define h and k

h = 1/m # h is grid size
k = t[1] - t0 # k is time step

# please note that r has changed, r1 and r2 are also new
r = k / (h ** 2)
r1 = 0.5 * r # r1 = k / (2 * h ** 2)
r2 = k / (2 * h) 

######################################################################
# %% define bdry and initial conditions. Just examples, to modify

# initial condition being a small cos func
def eta(x):
    return 0.005 * np.cos(x)


# tree center being a small cos func; outside with source term being a Gaussian func
def g0(t):
    return 0.005 * np.cos(t)


# bdry condition at tree bark ~ source term?
def g1(t):
    # gaussian(x, mu, sig)
    mu = 4 # start meansuring at 8 am, peak at 12
    sig = 1
    return np.exp(-np.power(t - mu, 2.) / (2 * np.power(sig, 2.))) + 3

#######################################################################

# %% define matrices in time stepping. Diagonal/super diagonal depend on ri
a = -r1 * np.ones(m-1) 
b = np.ones(m)
for j in range(m):
    b[j] = 1 + r2 / (j * h + h) + r1
c = np.ones(m-1)
for j in range(m-1):
    c[j] = -r2 / ( j * h + h) - r1

tridiag = sparse.diags([a, b, c], [-1, 0, 1], shape=(m, m)).toarray()

#%%

soln = []

# IC
U0 = eta(x)

#%%
#  Solving Ax1 = Bx0, this is B. Need to modify bdry condition each step
aa = r1 * np.ones(m-1) 
bb = np.ones(m)
for j in range(m):
    bb[j] = 1 - r2 / (j * h + h) - r1
cc = np.ones(m-1)
for j in range(m-1):
    cc[j] = r2 / ( j * h + h) + r1

B = sparse.diags([aa, bb, cc], [-1, 0, 1], shape=(m, m)).toarray()


soln.append(U0)

# %% main time stepping

for i in range(n-1): 
    rhs = B.dot(U0)
    rhs[0] = rhs[0] + r1 * (g0(t0 + i*k) + g0(t0 + (i+1)*k))
    rhs[-1] = rhs[-1] + (r2 / (m * h) + r1) * (g1(t0 + i*k) + g1(t0 + (i+1)*k))
    U1 = np.linalg.solve(tridiag, rhs)  # sparse.linalg.lsqr(tridiag, rhs)
    U0 = U1
    soln.append(U0)
#%%

soln_plot = np.asarray(soln)

print(soln_plot.shape)

print("max and min of soln at final step = ",np.max(soln_plot[-1,:]),np.min(soln_plot[-1,:]))

grid_point = 23
plt.plot(t, soln_plot[:,grid_point], '.r-')
plt.title('Temperature with Gaussian source term with time step n=%i, at grid point x=%i ' %(n,grid_point) )
plt.xlabel('time')
plt.ylabel('temperature distribution')
plt.show()

"""
Test Gaussian: https://stackoverflow.com/questions/14873203/plotting-of-1-dimensional-gaussian-distribution-function

from matplotlib import pyplot as plt
import numpy as np

def gaussian(x, mu, sig):
    return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

x_values = np.linspace(-3, 3, 120)
for mu, sig in [(-1, 1), (0, 2), (2, 3)]:
    plt.plot(x_values, gaussian(x_values, mu, sig))

plt.show()
"""
