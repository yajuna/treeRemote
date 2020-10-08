#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May  3 17:45:53 2020

@author: selinateng

This code calculates the source terms in Potter & Andresen.

Sign convention - the latitude and longitude are positive for North and East, 
respectively. The azimuth angle is 0 at North and positive towards the east. 
The zenith angle is 0 at vertical and positive towards the horizon.

Tacoma coordinates: 47.2529° N, 122.4443° W

Angles should be in radians.
"""
import numpy as np
from pylab import *
from sunPosition import sunpos
from datetime import datetime
import pytz
from pysolar.solar import *
from Gaussian import sourceTerms as curve # not sure what this does anymore


#%% Initial conditions for source vars
ro = 1              #Wood density (kg/m^3)
c = 1               #Specific heat (J/kg*K)
T_in = 1            #Temperature of tree interior (K)
k = 1               #Thermal conductivity (W/(m*K))
T_sfc = 1           #Temperature of tree surface (K)
T_air = 1           #Temperature of surrounding air (K)


T_air = curve(0)

L = 1               #Vertical height (m)
R = 0.09            #Tree radius (m)
theta1 = np.pi      #Difference between wind direction and the aspect of the surface point (radians)
u = 1               #Windspeed (m/s)   
tau = 0.7           #Atmospheric transmissivity (unitless) - at midlatitude 
                    #& low elevation, approx 0.76-0.81                   
eta = 0.8           #Atmospheric absorption parameter (unitless) - at 
                    #midlatitude & low elevation, approx 0.80-0.84                    
#i = 2 * np.pi       #Angle of incidence, the angle between the direction of 
                    #sunlight and the local normal to the tree's surface                   
alpha = 0.3         #Albedo of surface
LAT = 47.2529
LON = -122.4443     # Latitude and longitude of Tacoma
phi_0 = 0

#%% Init Conditions based on Potter/Andresen (pine log in Boulder, CO)
R = 0.09
L = 2.4
alpha = 0.3
latitude_deg = 40.015
longitude_deg = -105.27
tau = 0.78
nu = 0.80
ro = 686
c = 2470
k = 0.46

#%%
def convectiveHeatLoss():
    if (theta1 < (2 * np.pi)):   #Theta is the lesser of 2pi radians or theta1
        theta = theta1 
    else: 
        theta = 2 * np.pi
    h_free = 18.293 * abs(T_sfc - T_air) / (L * T_air**3)**0.25 * (T_air + 97.77) / np.math.sqrt(179.02 + T_air)
    h_forced = 3.458 * (T_air + T_sfc - 0.74)**0.49 * (u / (R * T_air))**0.5 * (1 - (theta / 90)**3)
    h = h_free + h_forced
    H = h * (T_sfc - T_air) 
    return H


#%%
def longWaveRadiation():
    sigma = 5.67 * 10**(-8)     #Stefan Boltzmann constant (W/(m^2*K^4))
    IR_out = sigma * T_sfc**4
    IR_in = sigma * T_air**4
    IR = IR_in - IR_out
    return IR

def solarRadiationHeating():
    date = pytz.timezone('America/Vancouver').localize(datetime.datetime.now())
    altitude_deg = get_altitude(latitude_deg, longitude_deg, date)
    radiation = pysolar.solar.radiation.get_radiation_direct(date, altitude_deg)
    return radiation

def sourceTerms(atSurface):
    if (atSurface):
        H = convectiveHeatLoss()
        S = solarRadiationHeating()
        IR = longWaveRadiation()
        tot = H + S + IR
        print(H, S, IR)
    else:
        tot = 0
    return tot

# Main function for testing output 
def main():
    sourceTerms(atSurface=True)

if __name__ == "__main__":
    main()
