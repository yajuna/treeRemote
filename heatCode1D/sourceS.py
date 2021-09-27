#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  8 17:08:32 2020

from combinedHeatsource_update1.py

Only give source in the south

@author: yajun
"""
import matplotlib.pyplot as plt
import numpy as np

albedo = .3  # "For albedo, we assumed a value of 0.3" - Potter Anderson page 3

#  South variables
insolationS = [50, 150, 250, 325, 350, 370, 350, 325, 250, 150, 75, 25, 10] # S_dir+S_dif 
netInfS = [0, -25, -35, -50, -75, -80, -80, -75, -50, -40, -25, -10, -5]  # IR_in-IR_out
convectionS = [0, -30, -60, -75, -100, -125, -125, -100, -90, -65, -50, -30, -10]  # H
#conductionS = [0, -100, -150, -190, -195, -200, -180, -130, -80, -30, -28, -10, -6]  # measured quantity??
totalS = []

#  Make sure t = [0, 1, 2, ... , 12]
#  because we need 13 on the t-axis for the 13 data points the Potter Anderson paper provides.
t = np.linspace(0, 12, 13)

# x is the number of points the linear interpolation will have
# currently set to be 10x the size of t
x = np.linspace(0, 12, 1000)

#  interpolate the south data
#  Also turn into an array so it sums how I want it to for totalS
insolationSI = (1-albedo)*np.array(np.interp(x, t, insolationS))
# Scaled by (1-alpha)
netInfSI = np.array(np.interp(x, t, netInfS))  # 
convectionSI = np.array(np.interp(x, t, convectionS))

# Sum all the lines up
totalS = insolationSI + netInfSI + convectionSI

## to grab values
def sourceTermsSvalue(i):
    return totalS[i]