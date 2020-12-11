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


albedo = .3  # "For albedo, we assumed a value of 0.3" - Potter Andresen page 3


# =============================================================================
# Finished
# Code up source terms other than convection
# Put everything inside the heat source class
# Iron out functionality for everything to work from one call only specifying
# the range of dates and (optionally) the latitude and longitude

# To do 
# Put variable windspeed in the convection
# Figure out a good place to store misc. input variables 
# Improve documentation!!

# Questions for meeting
# What are lines 162-167 for? @Yajun
# Recommendations for how to validate results (e.g. graph?)
# Any update on getting data from Protasio for surface temp? (currently using
# an array of 1's as a placeholder)
# =============================================================================
class HeatSource:
    """
    Takes in a range of dates to study + latitude and longitude of the study
    location (uses coordinates for Tacoma, WA by default)
    Returns a dictionary where each key is a date and each value is an array
    of hourly, weather-dependent, total source term values for that day from
    7am to 7pm
    Source terms include solar radiation (aka insolation), convection (free and
    forced), and infrared (aka black body radiation)
    """

    def __init__(self, daterange, latitude=47.2529, longitude=-122.4443):
        """
        """
        # Coordinates for Tacoma, WA are default
        self.latitude = latitude
        self.longitude = longitude
 
        # Make a dict where each key is the name of a date, and each value is
        # the array of source term values from 7am to 7pm on that date
        result = {}
        weather_data = self.get_weather_data(daterange)
        for date in daterange:
            # Get weather data for the day
            day_name = date.strftime('%b %d, %Y')
            weather_on_day = weather_data[day_name]
            source_terms = self.sourceTerms(weather_on_day)
            result[day_name] = source_terms
    

    def get_weather_data(self, daterange):
         """
         Retrieves weather data: clear sky radiation (diffuse+direct), temperature and windspeed
             Weather datasets - https://www.kaggle.com/selfishgene/historical-hourly-weather-data
             Solar radiation method - https://pysolar.readthedocs.io/en/latest/#
         """
         # Read in weather data from historical hourly weather dataset
         def dateparse(x): return pd.datetime.strptime(x, '%Y-%m-%d %H:%M:%S')
         data_path_temp = '/Users/selinateng/Google Drive/UW 2020-21/Tree/tnt_weather_data/temperature.csv'
         data_handle_temp = pd.read_csv(data_path_temp, usecols=['datetime','Seattle'], parse_dates=[
                                    'datetime'], date_parser=dateparse, engine='python').fillna(0)
         data_path_wind = '/Users/selinateng/Google Drive/UW 2020-21/Tree/tnt_weather_data/wind_speed.csv'
         data_handle_wind = pd.read_csv(data_path_wind, usecols=['datetime','Seattle'], parse_dates=[
                                    'datetime'], date_parser=dateparse, engine='python').fillna(0)
         data_handle = data_handle_temp.merge(data_handle_wind, left_on='datetime', right_on='datetime', how='outer')
         # Make a df of weather data for each date and add it to a dict of dfs
         weather_data_dict = {}
         end_delta = relativedelta.relativedelta(days=1)
         for date in daterange:
             end_date = date + end_delta
             # Clip off the minutes and seconds since data is hourly
             end_date = end_date.replace(minute=0, second=0, microsecond=0)
             mask = (data_handle['datetime'] >= date) & (data_handle['datetime'] < end_date)
             df = data_handle.loc[mask]
             df.columns=['measure_date','temperature','windspeed']
             # Add in solar radiation data from pysolar
             solar_radiation = []
             for hour in range(24):
                 date_timezone = date.replace(tzinfo=pytz.timezone('America/Los_Angeles')) + timedelta(hours=hour)
                 altitude_deg = sun.solar.get_altitude(self.latitude, self.longitude, date_timezone)
                 radiation = sun.radiation.get_radiation_direct(date_timezone, altitude_deg)
                 solar_radiation.append(radiation)
             df['solar_radiation'] = solar_radiation
             # Save the 24hr data separately as a csv file
             df.to_csv(r'Source Term Weather ' + date.strftime('%b-%d-%Y') + '.csv')
             # Take the 7am-7pm data and save it to the dict
             indexes_to_keep = range(7,20)
             df_sliced = df.take(list(indexes_to_keep))
             weather_data_dict[date.strftime('%b %d, %Y')] = df_sliced
             print(weather_data_dict[date.strftime('%b %d, %Y')])
         # Demo: plot the solar radiation on a given day 10/10/2017, 7am-7pm
         # df_to_show = weather_data_dict['Oct 10, 2017']
         # plt.plot(range(7,20), df_to_show['solar_radiation'][7:20])
         return weather_data_dict
     
        
    def insolation(self, today):
        return today['solar_radiation'].tolist()

    def netInfS(self, today):
        sb_const = 5.67e-8 # Stefan-Boltzmann constant (W/m^2/K^4)
        data_air_temp = today['temperature'].tolist()
        data_tree_temp = [1]*12 # placeholder
        air_4 = [n ** 4 for n in data_air_temp]
        tree_4 = [n ** 4 for n in data_tree_temp]
        netInfS = [(air_i - tree_i)*sb_const for air_i, tree_i in zip(air_4, tree_4)]
        #print(data_air_temp)
        #print(netInfS)
        # netInfS = [0, -25, -35, -50, -75, -80, -80, -75, -50, -40, -25, -10, -5]  # IR_in-IR_out
        return netInfS
    
    def convectionS(self, today):
        h=1 # convective heat transfer coefficient (placeholder)
        t_air = today['temperature'].tolist()
        t_sfc = [1]*12  # placeholder
        def h_free(t_s, t_a):
            return 18.293 * (abs(t_s-t_a)/(t_a)**3)**0.25 * (t_a+97.77)/(179.02+t_a)**0.5
        h_free = [h_free(t_s, t_a) for t_s, t_a in zip(t_air, t_sfc)]
        
        def h_forced(t_s, t_a):
            R = 1 # edit later
            R = 90 # also edit later
            u = 1 
            # Next step: make windspeed variable
            return 3.458 * (t_a+t_s-0.74)*80.49 * (u/R/t_a)**0.5 * (1 - (theta/90)**3)
        
        #convectionS = h*(t_sfc[7:20] - t_air[7:20])
        #convectionS = [0, -30, -60, -75, -100, -125, -125, -100, -90, -65, -50, -30, -10]
        #return convectionS
        return 0
    
    
    ## to grab values
    def sourceTerms(self, today):
        totalS = []
    
        #  Make sure t = [0, 1, 2, ... , 12]
        #  because we need 13 on the t-axis for the 13 data points the Potter Anderson paper provides.
        t = np.linspace(0, 12, 13)
        
        # x is the number of points the linear interpolation will have
        # currently set to be 10x the size of t
        x = np.linspace(0, 12, 1000)
        
        # This stuff needs to be updated 
        #  interpolate the south data
        #  Also turn into an array so it sums how I want it to for totalS
        #insolationSI = (1-albedo)*np.array(np.interp(x, t, insolationS))
        # Scaled by (1-alpha)
        #netInfSI = np.array(np.interp(x, t, netInfS))  # 
        #convectionSI = np.array(np.interp(x, t, convectionS))
        
        # Sum all the lines up
        ins = self.insolation(today)
        inf = self.netInfS(today)
        conv = self.convectionS(today)
        #totalS = insolation() + netInfS() + convectionS()
        
        #return totalS[i]
        return 0
    
    


# Driver
# Define range of dates to study
daterange = []
years = [2017]
months = [10]
dates = range(10,11)
for year in years:
    for month in months:
        for date in dates:
            daterange.append(datetime(year, month, date))
num_dates = len(daterange)
names_dates = []
for d in range(num_dates):
    names_dates.append(daterange[d].strftime('%b %d, %Y'))
sources = HeatSource(daterange)



