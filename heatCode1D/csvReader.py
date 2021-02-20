#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 16 11:21:10 2021

csv reader for TNT_weather data. 
Reads specific column, then convert to np array

@author: yajuna
"""
import pandas
import numpy as np

colnames = ['time', 'temperature', 'windspeed']

data = pandas.read_csv('TNT_weather_data_JP_Jan21.csv', names=colnames)

time = data.time.tolist()

temp = data.temperature[1:].tolist()
temp1 = np.asarray([float(t) for t in temp])

windspeed = data.windspeed[1:].tolist()
windspeed1 = np.asarray([float(v) for v in windspeed])


# test value of center tree temperature
Ta = np.random.uniform(low=21.5 + 273.2, high=25 + 273.2, size=(50,))