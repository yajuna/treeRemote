#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  6 13:12:44 2020

data obtained from source_datadriven_inprogress.py 
(line 26, starts from 7am)

Only compute source from solar radiation.

to replace source term in heat1dK.py

@author: yajun
"""
import numpy as np

solar = [290.35, 289.84, 287.94, 286.03, 284.84, 284.05, 283.55, 283.24, 282.93, 282.23, 281.54, 280.64, 280.54, 280.43, 280.14, 279.84, 280.36, 282.04, 283.94, 285.03, 284.59, 284.71, 286.16, 286.45]

solarArray = np.asarray(solar)

solar7to7 = solarArray[:13]

#  Make sure t = [0, 1, 2, ... , 12]
#  because we need 13 on the t-axis for the 13 data points the Potter Anderson paper provides.
t = np.linspace(0, 12, 13)

# x is the number of points the linear interpolation will have
# currently set to be 10x the size of t
x = np.linspace(0, 12, 1000)

solarSource = np.array(np.interp(x, t, solar7to7))

def solar(i):
    return solarSource[i]
