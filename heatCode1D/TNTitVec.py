# -*- coding: utf-8 -*-
"""
Created on Sun Aug 15 14:19:55 2021

Change legends to reflect plotted curves

@author: mhockman
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
from functools import partial
from scipy.optimize import broyden1
import datetime as dt
from matplotlib import pyplot as plt
import source as s


param = {"Ta": 300.4, "Va": 5, "qrads": 950, "Pr": 0.707,
         "Ka": 26.3 * 10 ** -3, "Kt": 0.11, "nu": 15.89 * 10 ** -6,
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
#Tb will hold the answer in Kelvin, TbFinal in Celcius
Tbinit = 300.4
file_name = 'TNTtemp.csv'
Tb = []
TbCelcius = []
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
    TaAAvg = sum(TaActual)/len(TaActual)
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
    #for Ta in TaActual2:  ####Uncomment to use variable core temp####
    for ws in s.windspeed: ####Uncomment to use variable windspeed####
#    for q in s.solar:   #### Uncomment to use variable solar rads ####
        #param["Ta"] = Ta   ####Uncomment to use variable core temp####
        param["Va"] = ws   ####Uncomment to use variable windspeed####
#        param["qrads"] = q  #### Uncomment to use variable solar rads ####
        tempSltn = broyden1(partial(F, **param), Tbinit)
        Tb.append(tempSltn)
        TbCelcius.append(tempSltn-273.2)
            

#print(Tb)
#print(TbCelcius)

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
TaAnp = np.asarray(TaActual2)
#error = np.asarray(error)

#print("max error in Tb and TbActual", error.max())


tActual = np.linspace(0,24,len(TbActual),endpoint = False)

### TESTING WITH VARIABLE CORE TEMPERATURE
# plt.plot(tt, Tbnp, 'r', label = "Estimated")
# plt.plot(tt, TaAnp, 'g', label = "Measured Core Temperature")
# plt.plot(tActual, TbAnp,'b', label = "Measured") #

# #plt.plot(t, error, 'g', label = "Error")

# plt.legend(loc = 'lower left')
# plt.title('Estimated Bark Temperature vs Measured Bark Temperature')
# plt.axis([0,24,295,305])
# plt.xlabel('Time (hrs)')
# plt.ylabel('Temperature ()')
# #plt.savefig('/home/yajun/Documents/treePower/TNT/TNTfig/' + 'TNTiterationVec' + '.eps', format='eps', dpi=300,bbox_inches='tight')
# plt.show()


##### TESTING WITH VARIABLE WINDSPEED
plt.plot(tt, Tbnp, 'r', label = "Estimated")
plt.plot(tActual, TbAnp, 'b', label = "Measured")

plt.legend(loc = 'lower left')
plt.title('Estimated Bark Temperature vs Measured Bark Temperature')
plt.axis([0,24,280,305])
plt.xlabel('Time (hrs)')
plt.ylabel('Temperature ()')
plt.show()

plt.plot(tt, s.windspeed, 'g')
plt.title("Windspeed")
plt.axis([0, 24, 0, 6])
plt.show()

#### TESTING WITH VARIABLE SOLAR RADIATION
#plt.plot(tt, Tbnp, 'r', label = "Estimated")
#plt.plot(tActual, TbAnp, 'b', label = "Measured")
#plt.axis([0,24,297,303])
#
#plt.legend(loc = 'lower left')
#plt.title('Estimated Bark Temperature vs Measured Bark Temperature')
#plt.axis([0,24,297,303])
#plt.xlabel('Time (hrs)')
#plt.ylabel('Temperature ()')
#plt.show()
#
#plt.plot(tt, s.solar, 'g')
#plt.title("Solar Radiation")
#plt.axis([0, 24, 0, max(s.solar)])
#plt.show()