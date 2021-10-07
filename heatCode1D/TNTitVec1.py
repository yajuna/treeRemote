# -*- coding: utf-8 -*-
"""
Created on Sun Aug 15 14:19:55 2021
Modification of TNTiterationVec.py, with correct sampling consistent with heatMain.py

To replace TNTvec.py
modifications done:
    1. Replace param with updated values
    2. Replace Tbinit
    
    4. Added conversion of list to numpy array
    5. Added visualization

@author: mhockman
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
from functools import partial
from scipy.optimize import broyden1
import datetime as dt
from matplotlib import pyplot as plt

param = {"Ta": 300.4, "Va": 10, "qrads": 950, "Pr": 0.707, 
         "Ka": 26.3 * 10 ** -3, "Kt": 0.12, "nu": 15.89 * 10 ** -6,
         "epsilon": 0.8, "sigma": 5.67 * 10 ** -8, "C": 0.193,
         "m": 0.618, "rb": 0.2, "L": 10, "DeltaT": 1,
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
Tbinit = 300.40
file_name = 'TNTtemp.csv'
Tb = []
TbActual = []
t = []
with open(file_name) as file:
    # First line is header, [1:] takes all lines after header
    # lines = list of lists. For each list in lines,
        # list[0] = index of data entry, starting at 0
        # list[1] = date in YYYY-MM-DD format
        # list[2] = time in 12H format
        # list[3] = core temp in Kelvin = list[-2]
        # list[4] = bark temp in Kelvin = list[-1]
    lines = [line.split(',') for line in file][1:]
    TbActual = [float(line[-1]) for line in lines]
    TbActual = TbActual[367:782]
    TaActual = [float(line[-2]) for line in lines]
    TaActual = TaActual[367:782]
    t = [line[2] for line in lines]
    # Since time doesn't have AM and PM, have to break up the data into
    # 3 chunks and put them all together to get a 24hr+ sample period
    # with initial value t0 given by TNTtemp.csv first entry
    t1 = t[367:575]
    t0a = dt.datetime.strptime(t1[0], '%I:%M:%S')
    t0 = 1/60*t0a.minute+(1/60)**2*t0a.second
    t1a = [dt.datetime.strptime(i, '%I:%M:%S') for i in t1]
    t1 = [i.hour+1/60*i.minute+(1/60)**2*i.second for i in t1a]
    # second time chunk
    t2 = t[575:782]
    t1fa= dt.datetime.strptime(t2[0], '%I:%M:%S')
    t1f = 12
    t2a = [dt.datetime.strptime(i, '%I:%M:%S') for i in t2]
    t2 = [t1f+i.hour+1/60*i.minute+(1/60)**2*i.second for i in t2a]
    # final time sample
    t = t1+t2
    # tt = sample time we use in heatMain
    tt = np.linspace(0, 24, 1000, endpoint=False)
    # TaActual2 = core temp during interpolated for times tt
    TaActual2 = np.interp(tt, t, TaActual[0:len(t)])
    # solve for bark temperatures for sample times tt
    for Ta in TaActual2:
        param["Ta"] = Ta
        tempSltn = broyden1(partial(F, **param), Tbinit)
        Tb.append(tempSltn)
        
            

#Calculate difference between actual bark temperature values
error = []

for i in range(len(Tb)):
#    tempError = (Tb[i]-TbActual[i])/TbActual[i]
    tempError = Tb[i]-TbActual[i]
    error.append(tempError)

#Calculate difference between actual bark temperature values
minCalculated = min(Tb)
minActual = min(TbActual)
maxCalculated = max(Tb)
maxActual = max(TbActual)

errorMin = (minCalculated-minActual)/minActual*100
errorMax = (maxCalculated-maxActual)/maxActual*100

print("min value error = {:0.4f}%".format(errorMin))
print("max value error = {:0.4f}%".format(errorMax))

Tbnp = np.asarray(Tb)
TbAnp = np.asarray(TbActual)
error = np.asarray(error)

print("max error in Tb and TbActual", error.max())

print("size of Tb and TbActual", Tbnp.size, TbAnp.size)

t = np.linspace(0,24,len(Tb),endpoint = False)

print(TbAnp.size)

plt.plot(t, Tbnp, 'r', label = "Estimated")

plt.plot(t, TbAnp,'b', label = "Measured") #

#plt.plot(t, error, 'g', label = "Error")

plt.legend(loc = 'lower left')
plt.title('Estimated Bark Temperature vs Measured Bark Temperature')
plt.axis([0,24,295,305])
plt.xlabel('Time (hrs)')
plt.ylabel('Temperature ()')
#plt.savefig('/home/yajun/Documents/treePower/TNT/TNTfig/' + 'TNTiterationVec' + '.eps', format='eps', dpi=300,bbox_inches='tight')
plt.show()