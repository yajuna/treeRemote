#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Code to organize Protasio's data in a single CSV file
Data retrieved from: 
    https://admin.tago.io/public/dashboard/5fd231d5dc9ed600282fddb0/a6a2b07c-78c2-4ac0-b586-c594994cc614

Note that # rows varies slightly by day

It would be REALLY cool to add to this script so that it automatically gets
new data from the website up to the present day...
^^ See communication w/ Protasio on how he did this

Interpolation code from:
    https://stackoverflow.com/questions/51392012/convert-irregular-time-series-to-hourly-data-in-python-and-have-normal-distribut
    

Created on Fri Jan 15 11:39:32 2021
@author: selinateng
"""
import pandas as pd
from datetime import *
from dateutil import relativedelta
from scipy import interpolate


# Take in the data from Protasio's web tool and combine internal and external
# tree temperature data into one dataframe
data_core = '/Users/selin/Google Drive/The Numerical Tree (TNT)/External Temperature - Tree Bark (ºC)-1614190306194.csv'
data_bark = '/Users/selin/Google Drive/The Numerical Tree (TNT)/Internal Temperature (ºC)-1614190252423.csv'
output_filename = 'protasio_tree_data.csv'
df_core_temp = pd.read_csv(data_core)
df_bark_temp = pd.read_csv(data_bark)
df = df_core_temp.merge(df_bark_temp, left_on='Date and Time', right_on='Date and Time')
# Rename columns for simplicity
df.columns = ['datetime', 'core_temp', 'bark_temp']
print(df)
# Convert date and time info to datetime objects
df['datetime'] = pd.to_datetime(df['datetime']) #, format="%m/%d/%Y, %H:%M:%S %p"
print(type(df['datetime'][0]))
df['datetime'] = [datetime.date(i) for i in df['datetime']]
print(type(df['datetime'][0]))
# Convert temperatures to Kelvin
df['core_temp'] = df['core_temp'] + 273.15
df['bark_temp'] = df['bark_temp'] + 273.15
# Flip data into chronological order
df = df.iloc[::-1] 
# Interpolate to be hourly
#df_rs = df.set_index('datetime').resample('H').mean().interpolate('linear')
#print(df_rs)
# Save to CSV
df.to_csv(output_filename, index=True)


# Separate data from individual days, organized in a dictionary
data_dict = {}
# Here we specify dates to study
daterange = []
years = [2021]
months = [2]
dates = range(3,22)
for year in years:
    for month in months:
        for date in dates:
            daterange.append(datetime(year, month, date))
end_delta = relativedelta.relativedelta(days=1)

for date in daterange:
    end_date = date + end_delta
    end_date = end_date.replace(minute=0, second=0, microsecond=0)
    #print(date, end_date)
    mask = (df['datetime'] >= date) & (df['datetime'] < end_date)
    data_on_date = df.loc[mask]

    data_dict[date.strftime('%b %d, %Y')]=data_on_date

   
# Demo call
#print(data_dict['Feb 21, 2021'])




