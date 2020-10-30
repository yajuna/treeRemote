#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  8 17:08:32 2020
from combinedHeatsource_update1.py
Only give source in the south
@author: yajun
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pysolar as sun
from datetime import *
from dateutil import relativedelta
import pytz
from HeatSource import *


albedo = .3  # "For albedo, we assumed a value of 0.3" - Potter Anderson page 3


#  South variables

def insolation():
    data = today['solar_radiation'].tolist()
    insolationS = data[7:20] # 7am to 7pm
    print(insolationS)
    return insolationS


def netInfS():
    sb_const = 5.67e-8 # Stefan-Boltzmann constant (W/m^2/K^4)
    data_air_temp = today['temperature'].tolist()
    data_tree_temp = [0]*12
    air_4 = [n ** 4 for n in data_air_temp]
    tree_4 = [n ** 4 for n in data_tree_temp]
    netInfS = [(air_i - tree_i)*sb_const for air_i, tree_i in zip(air_4, tree_4)]
    print(data_air_temp)
    print(netInfS)
    netInfS = [0, -25, -35, -50, -75, -80, -80, -75, -50, -40, -25, -10, -5]  # IR_in-IR_out
    return netInfS


def convectionS():
    h=1 # convective heat transfer coefficient (placeholder)
    data_air_temp = today['temperature'].tolist()
    data_tree_temp = [0]*12
    convectionS = h*(data_tree_temp[7:20] - data_air_temp[7:20])
    convectionS = [0, -30, -60, -75, -100, -125, -125, -100, -90, -65, -50, -30, -10]
    return convectionS


## to grab values
def sourceTermsSvalue(i, source_type):
    totalS = []

    #  Make sure t = [0, 1, 2, ... , 12]
    #  because we need 13 on the t-axis for the 13 data points the Potter Anderson paper provides.
    t = np.linspace(0, 12, 13)
    
    # x is the number of points the linear interpolation will have
    # currently set to be 10x the size of t
    x = np.linspace(0, 12, 1000)
    
    #  interpolate the south data
    #  Also turn into an array so it sums how I want it to for totalS
    insolationSI = (1-albedo)*np.array(np.interp(x, t, insolationS))
    # Scaled by (1-alpha)
    netInfSI = np.array(np.interp(x, t, netInfS))  # 
    convectionSI = np.array(np.interp(x, t, convectionS))
    
    # Sum all the lines up
    totalS = insolation() + netInfS() + convectionS()
    
    return totalS[i]

netInfS()
print(totalS[0:13])

