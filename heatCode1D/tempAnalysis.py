#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 20 20:27:24 2021

data analysis of heatMain.py

@author: yajun
"""

from heatMain import *
c = temp(config)

# plot temp at grid point 5 throughout 24 hours, length 1000
c0 = c[:,5]
# plot temp at time step 5 for all grid points, length 160

i = 5
for j in range(3):
    tempInd = i * 10**j
    tempC = c[tempInd, :]
    print(tempC)





