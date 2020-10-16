#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 16 12:06:05 2020

This code computes sensitivity of dependence on thermal conductivity. Code to study: heat1dK.py

@author: yajun
"""

import numpy as np
from scipy import sparse
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D # matplotlib version 3.1.0
# import source terms and bdry conditions
from sourceS import * # source term at bdry
from Temp_dataVec import * # boundary at tree bark
#%%
config = dict()
config['gridPoints'] = 50
config['timeSteps'] = 1000
mu,sigma = 1, 0.001 # mean and standard deviation
config['thermalConductivity'] = 0.12*np.random.normal(mu, sigma, 6)
config['heatCapacity_rhoc'] = 1.7
config['at_point'] = 38
config['time'] = np.linspace(0, 1000, 50, endpoint = False)

import heat1dK as h1

config['thermalConductivity'] = 0.12 * np.ones(12)

# constant thermal conductivity
c0 = h1.temp(config)

sigma = np.linspace(0.01, 0.09, 9)

c = []

for j in range(sigma.size):
    config['thermalConductivity'] = 0.12*np.random.normal(mu, sigma[j], 12)
    c.append(h1.temp(config))
    
sensitivity = []    

for j in range(sigma.size):
    sensitivity.append(np.max(np.abs((c0 - c[j])/c0)))
    
    
    



