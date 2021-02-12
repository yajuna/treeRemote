# -*- coding: utf-8 -*-
"""
Created on Wed Feb 10 15:07:51 2021

@author: Mhock
"""


from TNTiteration3 import *
import scipy.optimize
from functools import partial
import time
from random import *

#Setting up variables
param = {"Ta": 29 + 273.2, "Va": 10, "qrads": 650, "Pr": 0.707, "Ka": 26.3 * 10 ** -3, "Kt": 0.11,
         "nu": 15.89 * 10 ** -6, "epsilon": 0.8, "sigma": 5.67 * 10 ** -8, "C": 0.193, "m": 0.618, "rb": 0.2,
         "L": 10, "DeltaT": 2, "DeltaR": 100 / 1000}
TbB = []
tB = []
TbR = []
tR = []
n = 10**5

#Testing Broyden1 method
for i in range(n):
    #Random temperature guess between 100 and 400 degrees kelvin
    TempGuess = randrange(100, 400)
    #Randomly change parameter Ta by +-20%
    param["Ta"] = param["Ta"]*(1+0.2*randrange(-1, 1))
    tic = time.perf_counter()
    TbB.append(broyden1(partial(F, **param), TempGuess))
    toc = time.perf_counter()
    tB.append(toc-tic)

#Calculate average time it takes for broyden1 method to solve
tAvgB = sum(tB)/len(tB)


#Resetting Ta to it's original value
param["Ta"] = 29+293.2

#Testing root method, basically the same as the first loop
for i in range(n):
    TempGuess = randrange(100, 400)
    param["Ta"] = param["Ta"]*(1+0.2*randrange(-1, 1))
    tic = time.perf_counter()
    TbR.append(scipy.optimize.root(partial(F, **param), TempGuess))
    toc = time.perf_counter()
    tR.append(toc-tic)

#Calculate average time for root method
tAvgR = sum(tR)/len(tR)

#Print values
print("Average time for Bryoden1 =", str(tAvgB))
print("Average time for optimize.root =", str(tAvgR))