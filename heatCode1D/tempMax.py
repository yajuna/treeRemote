#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 23 12:04:09 2020

This code searches for the time within 12 hours where the max of temp difference
between center and bark of tree occurs

Study time that gives max temp difference. Variable gridPoints.

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
config['gridPoints'] = 750
config['timeSteps'] = 1000
mu,sigma = 1, 0.001 # mean and standard deviation
config['thermalConductivity'] = 0.12*np.ones(12) #0.12*np.random.normal(mu, sigma, 6)
config['heatCapacity_rhoc'] = 1.7
config['at_point'] = 38
config['time'] = np.linspace(0, 1000, 50, endpoint = False)

import heat1dK as h1


c = h1.temp(config)

print("shape of array ", c.shape) # (1000,gridPoint)

tempDiff = np.zeros(c[0].shape[0]) # number of rows, number of time steps

for j in range(tempDiff.size):
    tempDiff[j] = np.abs(c[j,0] - c[j,-1])
    
MaxTempIndex = np.argmax(tempDiff) 

print("Max temperature difference occurs at", MaxTempIndex, "time step", ", the difference is", tempDiff[MaxTempIndex])   
    

# loop for different gridPoints. Ran out of memory
"""
c = []
for j in range(100):
    config['gridPoints'] = config['gridPoints'] + 50
    c.append(h1.temp(config))
    
print("shape of each array, and length of list", c[0].shape, len(c))    

"""











