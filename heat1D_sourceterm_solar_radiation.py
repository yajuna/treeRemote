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
import pytz
import matplotlib.pyplot as plt

# Coordinates for Tacoma, WA
latitude = 47.2529
longitude = 122.4443

# Define range of dates to study
daterange = []
years = [2020]
months = [10]
dates = range(10,15)
for year in years:
    for month in months:
        for date in dates:
            daterange.append(datetime(year, month, date, tzinfo = timezone(-timedelta(hours=-7))))
num_dates = len(daterange)

# Retrieve clear sky radiation data for range of dates
radiation_data = []
for date in daterange:
    radiation_on_date = []
    for hour in range(24):
        date_and_time = date + timedelta(hours=hour)
        altitude_deg = sun.solar.get_altitude(latitude, longitude, date_and_time)
        radiation = sun.radiation.get_radiation_direct(date_and_time, altitude_deg)
        radiation_on_date.append(radiation)
    radiation_data.append(radiation_on_date)

# Visualize clear sky radiation data
for n in range(num_dates):
    plt.plot(range(24), radiation_data[n])
    plt.xlabel('Time (Hr)')
    plt.ylabel('Clear sky radiation (W/m^2)')
    text_to_show = 'Clear Sky Radiation on ' + daterange[n].strftime('%b %d, %Y')
    plt.title(text_to_show)
    plt.savefig('solar_radiation_graph')
    
    
    

    

        

# =============================================================================
# # Get the angle between the sun and a plane tangent to the earth where you are
# date = pytz.timezone('America/Vancouver').localize(datetime.datetime.now())
# print(get_altitude(latitude, longitude, date))
# 
# date = datetime.datetime(2007, 2, 18, 15, 13, 1, 130320, tzinfo=datetime.timezone.utc)
# print(get_altitude(latitude, longitude, date))
# 
# # Get azimuth of the sun
# date = datetime.datetime(2007, 2, 18, 15, 13, 1, 130320, tzinfo=datetime.timezone.utc)
# print(get_azimuth(latitude, longitude, date))
# 
# # Get estimate for clear-sky radiation
# latitude_deg = latitude # positive in the northern hemisphere
# longitude_deg = longitude # negative reckoning west from prime meridian in Greenwich, England
# date = datetime.datetime(2007, 2, 18, 15, 13, 1, 130320, tzinfo=datetime.timezone.utc)
# altitude_deg = get_altitude(latitude_deg, longitude_deg, date)
# radiation = radiation.get_radiation_direct(date, altitude_deg)
# print(radiation)
# =============================================================================

