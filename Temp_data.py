#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Michael Hockman
8/2/2020
Takes the temperature data and creates a function out of it
"""

import matplotlib.pyplot as plt
import numpy as np

T = np.array([59, 58, 57, 56, 55, 55, 55, 56, 57, 59, 62, 64, 67, 69, 71, 73, 74, 75, 75, 74, 70, 66, 63, 62, 60])
T = (T-32)*5/9+273.15
#  print(np.size(T))

#  Defines each point on the x-axis in which we have a corresponding data point from T
t = np.linspace(0, 24, 25)
#  print("t = " + str(t))


#  m is the number of grid points from heat1D Polar
#  t0 is a time to evaluate T at
#  Returns T (temp) evaluated at time t0
def temp(m, t0):
    #  Splits the domain 0 <= t <= 24 into m grid points
    x = np.linspace(0, 24, m)

    #  Creates the linear interpolation of the data
    y = np.interp(x, t, T)
    #  print("y = " + str(y))
    #  print("size y = " + str(np.size(y)))

    #  z is the length between data points on y. We are solving the equation z*m = 24
    #  which is saying (length between data points)*(total # of data points) = (last number in domain of x)
    z = 24 / m
    #  print("z = " + str(z))

    #  Solve the equation iz = t0 for i to find the index we need to evaluate y at
    #  rounded because you can't have a decimal index
    i = round(t0 / z)
    #  print("i = " + str(i))
    """
    #  Plot the full temperature data
    plt.plot(x[:], y[:])
    plt.title('1 day temperature, Tacoma WA')
    plt.xlabel('Hour')
    plt.ylabel('Temperature (F)')
    plt.xlim(0, 24)
    plt.ylim(280, 300)
    plt.show()
    """

    #  Return the temperature data from the correct index and convert to Kelvin
    return y[i]


#  Check to see if it gives good results
Z = temp(1000, 1)
print(Z)


def tempArray(t1):
    x = np.linspace(0, 24, np.size(t1)+1)

    #  Creates the linear interpolation of the data
    y = np.interp(x, t, T)
    print("size y = " + str(np.size(y)))

    print("size range = " + str(np.arange(1, np.size(t1))))
    return y[np.arange(1, np.size(t1)+1)]