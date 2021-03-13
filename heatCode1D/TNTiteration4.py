#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 28 21:41:54 2021

Problem with *arg variables combined with partial. 
Does not feed correct value.

Investigate when I have time. TNTiteration2.py works with *kwargs

@author: yajuna
"""
import numpy as np
from functools import partial
from scipy.optimize import broyden1
"""
param = {"Ta": 29 + 273.2, "Va": 10, "qrads": 650, "Pr": 0.707, "Ka": 26.3 * 10 ** -3, "Kt": 0.11,
         "nu": 15.89 * 10 ** -6, "epsilon": 0.8, "sigma": 5.67 * 10 ** -8, "C": 0.193, "m": 0.618, "rb": 0.2,
         "L": 10, "DeltaT": 2, "DeltaR": 100 / 1000}
"""
param = [29 + 273.2, 10, 650, 0.707, 26.3e-3, 0.11,
          15.89e-6, 0.8, 5.67e-8, 0.193, 0.618, 0.2,
          10, 2, 100 / 1000]
def F(Tb, *param):
    Ta, Va, qrads, Pr, Ka, Kt, nu, epsilon, sigma, C, m, rb, L, DeltaT, DeltaR = param
   
    r1 = rb - DeltaR
    Re = Va * (2 * rb) / nu
    Nu = C * (Re ** m) * Pr ** (1 / 3)  # Nusselt number
    h = Nu * Ka / (rb * 2)  # W/m2-k heat transfer coefficient
    Rcond = np.log(rb / r1) / (2 * np.pi * L * Kt)  # K/W
    qcond = DeltaT / Rcond  # W
    qeq = qcond-qrads
    return Ta + qeq / (2 * np.pi * L * h * rb + 
                       2 * L * epsilon * rb * sigma * np.pi * (Ta ** 3 + Ta ** 2 * Tb + Ta * Tb ** 2 + Tb ** 3)) - Tb
    

Tbinit = 22 + 273.2
F_partial = partial(F, *param)
Tb = broyden1(F_partial, Tbinit)
Tbfinal = Tb - 273.2