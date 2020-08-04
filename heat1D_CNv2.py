#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 27 12:47:12 2020

1D code for heat equation, with Crank-Nicolson

Ref: LeVeque Chapter 9

u_t=ku_xx

u(x,0)=eta(x)
u(0,t)=g0(t)
u(1,t)=g1(t)

xi = ih, tn = nk, h = Delta x, k = Delta t

centered difference:
(U^n+1_i - U^n_i)/k = 1/h^2 (U^n_i-1 - 2U^n_i + U^n_i+1)

CN:
(U^n+1_i - U^n_i)/k = 1/2 (D^2U^n_i + D^2U^n+1_i) 
=1/2h^2(U^n_i-1 - 2U^n_i + U^n_i+1 + U^n+1_i-1 - 2U^n+1_i + U^n+1_i+1)   

rewrite to be (r=k/2h^2)
-r U^n+1_i-1 + (1+2r) U^n+1_i - r U^n+1_i+1 = 
r U^n_i-1 + (1-2r) U^n_i + r U^n_i+1

stable for any k>0. Take k = O(h)
@author: yajun
"""

import numpy as np
from scipy import sparse

# number of grid points, space grid points x
m = 5
x = np.linspace(0, 1, m)

# number of time steps
n = 100
t0 = 0

# define h and k

h = 1/m #changed this so x always starts and ends at 0 and 1 no matter what m is
k = 1/2*h ** 2 #changed this to ensure k <= 1/2*h^2


r = k / (2 * h ** 2)


# %% define bdry and initial conditions. Just examples, to modify

# initial condition being a small cos func
def eta(x):
    return 0.005 * np.cos(x)


# tree center being a small cos func; outside being a big cos func
def g0(t):
    return 0.005 * np.cos(t)


def g1(t):
    return 0.05 * np.cos(t)


# %% define matrices in time stepping. eqn (9.9)
tridiag = sparse.diags([-r, 1 + 2 * r, -r], [-1, 0, 1], shape=(m, m)).toarray()

# %% main time stepping


soln = []
U0 = eta(x)
rhs = r * np.roll(U0, -1) + (1 - 2 * r) * U0 + r * np.roll(U0, 1)
rhs[0] = r * (g0(t0) + g0(t0 + k)) + (1 - 2 * r) * U0[0] + r * U0[1]
rhs[-1] = r * (g1(t0) + g1(t0 + k)) + (1 - 2 * r) * U0[-1] + r * U0[-2]
soln.append(U0)
print("soln for t = 0: " + str(soln));
for i in range(n):
    print()
    U1 = np.linalg.solve(tridiag, rhs)  # sparse.linalg.lsqr(tridiag, rhs)
    U0 = U1
    soln.append(U0)
    print("soln at " + str(i*k) + "< t < " + str((i+1)*k) + " = " + str(soln[i+1])) #changed to display solution for each time interval
    rhs = r * np.roll(U0, -1) + (1 - 2 * r) * U0 + r * np.roll(U0, 1)
    rhs[0] = r * (g0(t0+i*k) + g0(t0 + (i+1)*k)) + (1 - 2 * r) * U0[0] + r * U0[1] #changed what g0 is evaluated at
    rhs[-1] = r * (g1(t0+i*k) + g1(t0 + (i+1)*k)) + (1 - 2 * r) * U0[-1] + r * U0[-2] #changed what g1 is evaluated at

#print(soln)