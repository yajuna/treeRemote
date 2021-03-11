#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 16 11:21:10 2021

csv reader for TNT data. 
Reads specific column, then convert to np array

first portion is reference for format

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

#%% Use data from Feb 16; Feb 15 data maybe incomplete

colnames = ['time', 'temperature', 'windspeed', 'solar']

dataWeather = pandas.read_csv('weatherJoaoPessoaFeb15-19.csv', names=colnames)

time16 = dataWeather.time.tolist()

temp16 = dataWeather.temperature[5:28].tolist() # row 6 ro 29
temp16np = np.asarray([float(t) for t in temp16])

windspeed16 = dataWeather.windspeed[5:28].tolist()
windspeed16np = np.asarray([float(v) for v in windspeed16])

solar16 = dataWeather.solar[5:28].tolist()
solar16np = np.asarray([float(v) for v in solar16])