#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
from functools import partial
from scipy.optimize import broyden1


param = {"Ta": 29 + 273.2, "Va": 10, "qrads": 650, "Pr": 0.707, 
         "Ka": 26.3 * 10 ** -3, "Kt": 0.11, "nu": 15.89 * 10 ** -6,
         "epsilon": 0.8, "sigma": 5.67 * 10 ** -8, "C": 0.193,
         "m": 0.618, "rb": 0.2, "L": 10, "DeltaT": 2,
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
#Tb will hold the answer in Kelvin, TbFinal in Celcius
Tbinit = 22 + 273.2
file_name = 'protasio_tree_data.csv'
Tb = []
TbFinal = []
TbActual = []
with open(file_name) as file:
    for line in file:
        tempLine = line.split(',')
        if(tempLine[0] != 'datetime'):
            if(tempLine[-1][-1] == '\n'):
                tempLine[-1] = tempLine[-1][0:-1]
            tempTb = float(tempLine[-1])
            TbActual.append(tempTb)
            tempTa = float(tempLine[-2])
            param["Ta"] = tempTa
            tempSltn = broyden1(partial(F, **param), 100)
            Tb.append(float(tempSltn))
            TbFinal.append(tempSltn-273.2)
            

#print(Tb)
#print(TbFinal)

#Calculate difference between actual bark temperature values
error = []

for i in range(len(Tb)):
    tempError = (Tb[i]-TbActual[i])/TbActual[i]
    error.append(tempError)

print(error)