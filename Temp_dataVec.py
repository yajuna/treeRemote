#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tuesday 04 08 09:32:45 2020

takes temp measurements from Tacoma, 7 am to 7 pm. Interpolate according to combinedHeatsource_updated.py by Michael Hockman. Use temp at 7 am as initial condition, use interpolation as boundary condition

@author: yajun

"""

import matplotlib.pyplot as plt
import numpy as np

# hourly temperature from 7 o'clock to 7 o'clock, converted to kelvin
T = np.array([56, 57, 59, 62, 64, 67, 69, 71, 73, 74, 75, 75, 74])
T = (T-32)*5/9+273.15

# same variables as combinedHeatsource_updated.py

#  Make sure t = [0, 1, 2, ... , 12]
#  because we need 13 on the t-axis for the 13 data points the Potter Anderson paper provides.
t = np.linspace(0, 12, 13)

# x is the number of points the linear interpolation will have
# currently set to be 10x the size of t
x = np.linspace(0, 12, 1000)

# interpolate temp data, length 1000

tTemp = np.array(np.interp(x, t, T))

"""
plt.plot(x[:], tTemp[:], "-c", label="air temp")
plt.axis([0,12,0,350])
plt.legend()
plt.xlabel("Time since 7:00am (hrs)")
plt.ylabel("Air temp (kelvin)")
#plt.grid()
plt.show()
"""


def TacoTemp(t0):
    t0 = int(t0 / 0.012)
    return tTemp[t0]

## temperature of air at outer boundary                
def TacomaTemp(i):
    return tTemp[i]


