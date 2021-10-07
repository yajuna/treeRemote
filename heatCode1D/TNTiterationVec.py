#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 12 13:17:17 2021

@author: mhockman
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
from functools import partial
from scipy.optimize import broyden1
from matplotlib import pyplot as plt

# Ta, Va, qrads-TNTweather
#Ta = 29+273.2

param = {"Ta": 300.4, "Va": 5, "qrads": 950, "Pr": 0.707, 
         "Ka": 26.3 * 10 ** -3, "Kt": 0.12, "nu": 15.89 * 10 ** -6,
         "epsilon": 0.8, "sigma": 5.67 * 10 ** -8, "C": 0.193,
         "m": 0.618, "rb": 0.18, "L": 10, "DeltaT": 1,
         "DeltaR": 100 / 1000}

def F(Tb, **param):
    Ta = param["Ta"]
    Va = param["Va"]
    qrads = param["qrads"]
    Pr = param["Pr"]
    Ka = param["Ka"]
    Kt = param["Kt"]
    nu = param["nu"]
    epsilon = param["epsilon"]
    sigma = param["sigma"]
    C = param["C"]
    m = param["m"]
    rb = param["rb"]
    L = param["L"]
    DeltaT = param["DeltaT"]
    DeltaR = param["DeltaR"]
    r1 = rb - DeltaR
    Re = Va * (2 * rb) / nu
    Nu = C * (Re ** m) * Pr ** (1 / 3)  # Nusselt number
    h = Nu * Ka / (rb * 2)  # W/m2-k heat transfer coefficient
    Rcond = np.log(rb / r1) / (2 * np.pi * L * Kt)  # K/W
    qcond = DeltaT / Rcond  # W
    qeq = qcond-qrads
    return Ta + qeq / (2 * np.pi * L * h * rb + 
                       2 * L * epsilon * rb * sigma * np.pi * (Ta ** 3 
                        + Ta ** 2 * Tb + Ta * Tb ** 2 + Tb ** 3)) - Tb


#Going to go through Salina's file and grab all the core temperatures
#and use them to find the bark temperatures, and store them all in 
#a vector. 
#Tb will hold the answer in Kelvin
Tbinit = 300.40
file_name = 'TNTtemp.csv'
# computed bark temp
Tb = []
# measured bark temp
TbActual = []
with open(file_name) as file:
    for line in file:
        tempLine = line.split(',')
        if(tempLine[1] != 'date'):
            if(tempLine[-1][-1] == '\n'):
                tempLine[-1] = tempLine[-1][0:-1]
            tempTb = float(tempLine[-1])
            TbActual.append(tempTb)
            tempTa = float(tempLine[-2])
            print('Ta is',tempTa,'Tb is',tempTb, 'line is', line)
            param["Ta"] = tempTa
            tempSltn = broyden1(partial(F, **param), Tbinit)
            Tb.append(float(tempSltn))
            

#Calculate difference between actual bark temperature values
error1 = []

for i in range(len(Tb)):
#    tempError = (Tb[i]-TbActual[i])/TbActual[i]
    tempError = Tb[i]-TbActual[i]
    error1.append(tempError)
    
Tbnp = np.asarray(Tb)
TbAnp = np.asarray(TbActual)
error = np.asarray(error1)
#%% visualize

t = np.linspace(0,24,len(Tb),endpoint = False)

plt.plot(t, Tbnp, 'r', label = "Estimated")

plt.plot(t, TbAnp,'b', label = "Measured") #

#plt.plot(t, error, 'g', label = "Error")

plt.legend(loc = 'lower left')
plt.title('Estimated Bark Temperature vs Measured Bark Temperature')
plt.axis([0,24,295,305])
plt.xlabel('Time (hrs)')
plt.ylabel('Temperature (K)')
plt.savefig('/home/yajun/Documents/treePower/TNT/TNTfig/' + 'TNTiterationVec' + '.eps', format='eps', dpi=300,bbox_inches='tight')
plt.show()  
