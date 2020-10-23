#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 15 23:11:53 2020

@author: selinateng
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  8 10:46:29 2020

@author: selinateng
"""

import pysolar as sun
from datetime import *
from dateutil import relativedelta
import pytz
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

class HeatSources:
    
    # Coordinates for Tacoma, WA
    latitude = 47.2529
    longitude = -122.4443


    def __init__(self): # pass in the latitude, longitude, and dates
        # Define range of dates to study
        daterange = []
        years = [2017]
        months = [10]
        dates = range(10,11)
        for year in years:
            for month in months:
                for date in dates:
                    daterange.append(datetime(year, month, date)) #tzinfo = timezone(-timedelta(hours=-7))))
        num_dates = len(daterange)
        names_dates = []
        for d in range(num_dates):
            names_dates.append(daterange[d].strftime('%b %d, %Y'))


    def get_weather_data(self):
         """
         Retrieve weather data: temperature and windspeed
         Weather datasets available at
             https://www.kaggle.com/selfishgene/historical-hourly-weather-data
         Solar radiation data available from pysolar
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
             print(date)
             end_date = date + end_delta
             # Clip off the minutes and seconds since data is hourly
             end_date = end_date.replace(minute=0, second=0, microsecond=0)
             mask = (data_handle['datetime'] >= date) & (data_handle['datetime'] < end_date)
             df = data_handle.loc[mask]
             df.columns=['measure_date','temperature','windspeed']
             weather_data_dict[date.strftime('%b %d, %Y')]=df
             # Add solar data from pysolar
             solar_radiation = []
             for hour in range(24):
                 date_timezone = date.replace(tzinfo=pytz.timezone('America/Los_Angeles')) + timedelta(hours=hour)
                 #date_and_time = date(timezone(-timedelta(hours=-7))) + timedelta(hours=hour)
                 altitude_deg = sun.solar.get_altitude(latitude, longitude, date_timezone)
                 radiation = sun.radiation.get_radiation_direct(date_timezone, altitude_deg)
                 solar_radiation.append(radiation)
             df['solar_radiation'] = solar_radiation
             df.to_csv(r'Source Term Weather ' + date.strftime('%b-%d-%Y') + '.csv')
         # Demo: plot the solar radiation on a given day 10/10/2017
         df_to_show = weather_data_dict['Oct 10, 2017']
         plt.plot(range(24), df_to_show['solar_radiation'])
         
         return weather_data_dict
             

# Driver
sources = HeatSources()
weather = sources.get_data_temp_and_windspeed()
print(weather)


# =============================================================================
# # Visualize clear sky radiation data
# for n in range(num_dates):
#     plt.plot(range(24), radiation_data[n])
#     plt.xlabel('Time (Hr)')
#     plt.ylabel('Clear sky radiation (W/m^2)')
#     text_to_show = 'Clear Sky Radiation on ' + names_dates[n]
#     plt.title(text_to_show)
#     plt.savefig('solar_radiation_graph')
# =============================================================================