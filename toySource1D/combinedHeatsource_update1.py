#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thursday 08 20 10:28:12 2020
re interpretation of Fig.4 data from Potter and ANdresen. Based on Michael Hockman's code
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


#  South variables
insolationS = [50, 150, 250, 325, 350, 370, 350, 325, 250, 150, 75, 25, 10] # S_dir+S_dif 
netInfS = [0, -25, -35, -50, -75, -80, -80, -75, -50, -40, -25, -10, -5]  # IR_in-IR_out
convectionS = [0, -30, -60, -75, -100, -125, -125, -100, -90, -65, -50, -30, -10]  # H
#conductionS = [0, -100, -150, -190, -195, -200, -180, -130, -80, -30, -28, -10, -6]  # measured quantity??
totalS = []

#  interpolate the south data
#  Also turn into an array so it sums how I want it to for totalS
insolationSI = (1-albedo)*np.array(np.interp(x, t, insolationS))
# Scaled by (1-alpha)
netInfSI = np.array(np.interp(x, t, netInfS))  # 
convectionSI = np.array(np.interp(x, t, convectionS))

# Sum all the lines up
totalS = insolationSI + netInfSI + convectionSI



print(totalS)
print(type(totalS))
print(np.size(totalS))



#  Uncomment if you want to see the plot
plt.plot(x[:], totalS[:], "-g", label="south")
plt.plot(x[:], totalN[:], "-b", label="north")
plt.axis([0,12,-60,150])
plt.legend()
plt.xlabel("Time since 7:00am (hrs)")
plt.ylabel("Energy flux $(W/m^2)$")
plt.grid()
plt.savefig('/home/yajun/Documents/treePower/combinedHeatSource.eps', format='eps', dpi=300)
plt.show()



#  input: t0, the time at which the source term will be evaluated at
#  return: North, the total source terms at time t0
#  if it can't be evaluated at t0 exactly, it will find the closest available point to evaluate
def sourceTermsN(t0):
    t0 = int(t0 / 0.012)
    return totalN[t0]

## to grab values
def sourceTermsNvalue(i):
    return totalN[i]

# 0.008 is dt from length of t vector

#  input: t0, the time at which the source term will be evaluated at
#  return: S evaluated at t0, the total source terms at time t0
#  if it can't be evaluated at t0 exactly, it will find the closest available point to evaluate
def sourceTermsS(t0):
    t0 = int(t0/0.012)
    print(t0)
    return totalS[t0]

## to grab values
def sourceTermsSvalue(i):
    return totalS[i]

# Simply returns an array of size t instead of a single element at point t0
# t is an input to make it match with the Heat1D file
def sourceTermsArrayS(t):
    return totalS[np.arange(1, np.size(t))]


# Simply returns an array of size t instead of a single element at point t0
# t is an input to make it match with the Heat1D file
def sourceTermsArrayN(t):
    return totalN[np.arange(1, np.size(t))]