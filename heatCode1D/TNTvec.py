#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 18 15:24:43 2021

Question: "DeltaR": 100 / 1000 why is it?? change to 0.18/160, gridPoint = 160??

The following are modifications from TNTiteration2.py

input Ta as a vector: read center temperature as cR.coreTemp16np, 
then use Ta[j] to compute Tb

input Va as a vector: read windspeed cR.windspeed16np
then use Va[j] to compute Tb

gridPoint = 50, rb = radius of tree at bark = 0.2

##########
Main function to call is chtbt (convective heat transfer, back temp) 

EX. 
import TNTvec as tntV
h, Tb = tntV.chtbt()

## to generate h and bdry, directly - run TNTvec -

line 79: r1 = deltaX 0.18/int(gridPoint)
then DeltaT is temp diff corresponding to r1

diameter of mango tree about 36 cm, rb = 0.18

Conductivity computed by Hee-Seok

@author: yajuna
"""
import numpy as np
from functools import partial
from scipy.optimize import broyden1
from matplotlib import pyplot as plt

import csvReader as cR

"""
interpretate coreTemp16npinte and windspeed16npinte data, then

h = []
bdry = []
for j in range(interpretated.shape):
    param = {"Ta": coreTemp16npinte[j], "Va": windspeed16npinte[j], "qrads": 650, "Pr": 0.707, "Ka": 26.3e-3, "Kt": 0.11,
         "nu": 15.89e-6, "epsilon": 0.8, "sigma": 5.67e-8, "C": 0.193, "m": 0.618, "rb": 0.18,
         "L": 10, "DeltaT": 2, "DeltaR": 0.2 / 50, "timeSteps": 1000}
    h, Tb = tnt2.chtbt()
    h.append(h)
    bdry.append(Tb)
"""
#param = {"Ta": cR.coreTemp16np, "Va": cR.windspeed16np, "qrads": 650, "Pr": 0.707, "Ka": 26.3e-3, "Kt": 0.11,
#         "nu": 15.89e-6, "epsilon": 0.8, "sigma": 5.67e-8, "C": 0.193, "m": 0.618, "rb": 0.2,
#         "L": 10, "DeltaT": 2, "DeltaR": 100 / 1000, "timeSteps": 1000}

def F(Tb, **param):
    Ta = param["Ta"]
    qrads = param["qrads"]
    Kt = param["Kt"]
    epsilon = param["epsilon"]
    sigma = param["sigma"]
    L = param["L"]
    DeltaT = param["DeltaT"]
    DeltaR = param["DeltaR"]
    
    rb = param["rb"]
    
    Pr = param["Pr"]
    Ka = param["Ka"]
    C = param["C"]
    m = param["m"]
    Va = param["Va"]
    nu = param["nu"]
    
    Re = Va * (2 * rb) / nu
    Nu = C * (Re ** m) * Pr ** (1 / 3)  # Nusselt number
    h = Nu * Ka / (rb * 2)  # W/m2-k heat transfer coefficient
    
    r1 = rb - DeltaR
    Rcond = np.log(rb / r1) / (2 * np.pi * L * Kt)  # K/W
    qcond = DeltaT / Rcond  # W
    qeq = qcond - qrads
    return Ta + qeq / (2 * np.pi * L * h * rb + 
                       2 * L * epsilon * rb * sigma * np.pi * (Ta ** 3 + Ta ** 2 * Tb + Ta * Tb ** 2 + Tb ** 3)) - Tb

def heatTransferCoeff(**param):
    rb = param["rb"]
    
    Pr = param["Pr"]
    Ka = param["Ka"]
    C = param["C"]
    m = param["m"]
    Va = param["Va"]
    nu = param["nu"]
    
    Re = Va * (2 * rb) / nu
    Nu = C * (Re ** m) * Pr ** (1 / 3)  # Nusselt number
    return Nu * Ka / (rb * 2)  # W/m2-k heat transfer coefficient

Tbinit = 300.40 #22 + 273.2

def chtbt():
    h = heatTransferCoeff(**param)
    Tb = broyden1(partial(F, **param), Tbinit)
    Tbfinal = Tb # Tb - 273.2 if celsius
    return h, Tbfinal

#%% weather data collected
windSpeed = np.interp(np.linspace(0,24,1000), np.linspace(0,24,cR.windspeed16np.size),cR.windspeed16np)
albedo = .3  # "For albedo, we assumed a value of 0.3" - Potter Anderson page 3
radiation = (1 - albedo) * np.interp(np.linspace(0,24,1000), np.linspace(0,24,cR.solar16np.size), cR.solar16np)

hVec = []
bdry = []

#%% Tb generated with experimental data
barkTemp = np.interp(np.linspace(0,24,1000), np.linspace(0,24,cR.barkTemp16np.size),cR.barkTemp16np)
coreTemp = np.interp(np.linspace(0,24,1000), np.linspace(0,24,cR.coreTemp16np.size),cR.coreTemp16np)
diffTemp = barkTemp - coreTemp

for j in range(1000):
    param = {"Ta": coreTemp[j], "Va": windSpeed[j], "qrads": radiation[j], "Pr": 0.707, "Ka": 26.3e-3, "Kt": 0.12,
         "nu": 15.89e-6, "epsilon": 0.8, "sigma": 5.67e-8, "C": 0.193, "m": 0.618, "rb": 0.18,
         "L": 10, "DeltaT": diffTemp[j], "DeltaR": 100/1000}#, "timeSteps": 1000}
    h, Tb = chtbt()
    hVec.append(h)
    bdry.append(Tb)
    
hArray = np.asarray(hVec)
bdryArray = np.asarray(bdry)
hArrayEx = hArray
bdryArrayEx = bdryArray

print("h and Tb from experimental data")

#%% visualize

t = np.linspace(0,24,1000,endpoint = False)

plt.plot(t, bdryArrayEx, 'r', label = "Generated")

plt.plot(t,barkTemp,'b',label="Measured") #

plt.legend(loc = 'lower left')
plt.title('Generated Bark Temperature vs Measured Bark Temperature')
plt.axis([0,24,290,315])
plt.xlabel('Time (hrs)')
plt.ylabel('Temperature (K)')
plt.savefig('/home/yajun/Documents/treePower/TNT/TNTfig/' + 'PostempDiff' + 'ExBark' + '.eps', format='eps', dpi=300,bbox_inches='tight')
plt.show()    
    

#%% Run either lines 121-154, or lines 162-195


#%% Tb generated with estimate data (no measurements of tree temp)

#tempAir = np.interp(np.linspace(0,24,1000), np.linspace(0,24,cR.temp16np.size),cR.temp16np)
#deltaT = np.average(diffTemp)
#
#for j in range(1000):
#    param = {"Ta": tempAir[j], "Va": windSpeed[j], "qrads": radiation[j], "Pr": 0.707, "Ka": 26.3e-3, "Kt": 0.12,
#         "nu": 15.89e-6, "epsilon": 0.8, "sigma": 5.67e-8, "C": 0.193, "m": 0.618, "rb": 0.18,
#         "L": 10, "DeltaT": deltaT, "DeltaR": 100/1000, "timeSteps": 1000}    
#    h, Tb = chtbt()
#    hVec.append(h)
#    bdry.append(Tb)
#    
#hArray = np.asarray(hVec)
#bdryArray = np.asarray(bdry)
#
#hArrayEs = hArray
#bdryArrayEs = bdryArray
#
#print("h and Tb from estimate data")
#
##%% visualize
#
#t = np.linspace(0,24,1000,endpoint = False)
#
#plt.plot(t, bdryArrayEs, 'r', label = "Estimated")
#
#plt.plot(t,barkTemp,'b',label="Measured") #
#
#plt.legend(loc = 'lower left')
#plt.title('Estimated Bark Temperature vs Measured Bark Temperature')
#plt.axis([0,24,290,305])
#plt.xlabel('Time (hrs)')
#plt.ylabel('Temperature ()')
#plt.savefig('/home/yajun/Documents/treePower/TNT/TNTfig/' + 'tempDiff' + 'EsBark' + '.eps', format='eps', dpi=300,bbox_inches='tight')
#plt.show()    
    