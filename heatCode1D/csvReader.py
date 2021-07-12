#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 16 11:21:10 2021

csv reader for TNT data. Read Feb 16 from 0am to 24pm
Reads specific column, then convert to np array

TNTweather.csv has weather data from location, dating Feb 16-18: wind 
speed, air temperature, and solar radiation

TNTtemp.csv has temp date from location, dating Feb 16-18: core 
temperature and bark temperature

TNTweather: line 4-27, or 6-29?? Data is from Feb 16. 

@author: yajuna
"""
import pandas
import numpy as np

#%% Use data from Feb 16; Feb 15 data maybe incomplete

colnames = ['date', 'time', 'temperature', 'windspeed', 'solar']

dataWeather = pandas.read_csv('TNTweather.csv', names=colnames)

date16 = dataWeather.date[4:27].tolist()

time16 = dataWeather.time[4:27].tolist()

temp16 = dataWeather.temperature[4:27].tolist() # row 6 ro 29
temp16np = np.asarray([float(t) for t in temp16])

windspeed16 = dataWeather.windspeed[4:27].tolist()
windspeed16np = np.asarray([float(v) for v in windspeed16])

solar16 = dataWeather.solar[4:27].tolist()
solar16np = np.asarray([float(v) for v in solar16])

#%% Use data from Feb 16; Feb 15 and 19 data maybe incomplete

colnames = ['date1', 'time1', 'coreTemp', 'barkTemp']

dataTemp = pandas.read_csv('TNTtemp.csv', names=colnames)

datet16 = dataTemp.date1[367:781].tolist()

timet16 = dataTemp.time1[367:781].tolist()

coreTemp16 = dataTemp.coreTemp[367:781].tolist() # row 6 ro 29
coreTemp16np = np.asarray([float(t) for t in coreTemp16])

barkTemp16 = dataTemp.barkTemp[367:781].tolist()
barkTemp16np = np.asarray([float(v) for v in barkTemp16])
