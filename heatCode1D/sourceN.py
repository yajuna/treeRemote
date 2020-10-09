#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  8 17:10:41 2020

from combinedHeatsource_update1.py

Only give source in the north

@author: yajun
"""

import matplotlib.pyplot as plt
import numpy as np

albedo = .3  # "For albedo, we assumed a value of 0.3" - Potter Anderson page 3

#  North variables
insolationN = [28, 30, 29, 27, 26, 25, 26, 27, 28, 30, 28, 14, 0]  # S_dir+S_dif
netInfN = [13, -9, -7, -2, -5, -7, -8, -6, -7, -14, -16, -28, -8]  # IR_in-IR_out
convectionN = [13, -9, -7, -2, -5, -7, -8, -9, -16, -18, -33, -30, -22]  # H
#conductionN = [15, -10, -8, -11, -8, -3, -1, -2, 1, 7, 20, 22, 20]  # measured quantity??
totalN = []

#  Make sure t = [0, 1, 2, ... , 12]
#  because we need 13 on the t-axis for the 13 data points the Potter Anderson paper provides.
t = np.linspace(0, 12, 13)

# x is the number of points the linear interpolation will have
# currently set to be 10x the size of t
x = np.linspace(0, 12, 1000)

#  interpolate the north data
#  Also turn into an array so it sums how I want it to for totalN
insolationNI = (1-albedo)*np.array(np.interp(x, t, insolationN))
# Scaled by (1-alpha)
netInfNI = np.array(np.interp(x, t, netInfN))  
convectionNI = np.array(np.interp(x, t, convectionN))

# Sum all the lines up
totalN = insolationNI + netInfNI + convectionNI

def sourceTermsNvalue(i):
    return totalN[i]