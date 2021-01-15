#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Code to organize Protasio's data in a single CSV file
Data retrieved from: 
    https://admin.tago.io/public/dashboard/5fd231d5dc9ed600282fddb0/a6a2b07c-78c2-4ac0-b586-c594994cc614

Dates w/ 24hr data available so far: 1/9/2021 - 1/14/2021

It would be REALLY cool to add to this script so that it automatically gets
new data from the website up to the present day...

Created on Fri Jan 15 11:39:32 2021
@author: selinateng
"""
import pandas as pd
from datetime import *
from dateutil import relativedelta


# =============================================================================
# def dateparse(x): return pd.datetime.strptime(x, '%Y/%m/%d, %H:%M:%S AM')
# df_core = pd.read_csv('/Users/selinateng/Google Drive/protasio_core_temp.csv', parse_dates=[
#                       'Date and Time'], date_parser=dateparse, engine='python')
# df_bark = pd.read_csv('/Users/selinateng/Google Drive/protasio_bark_temp.csv', parse_dates=[
#                       'Date and Time'], date_parser=dateparse, engine='python')
# =============================================================================

# Take in the data from Protasio's web tool and combine internal and external
# tree temperature data into one dataframe
df_core_temp = pd.read_csv('/Users/selinateng/Google Drive/protasio_core_temp.csv')
df_bark_temp = pd.read_csv('/Users/selinateng/Google Drive/protasio_bark_temp.csv')
df = df_core_temp.merge(df_bark_temp, left_on='Date and Time', right_on='Date and Time')
# Rename columns for simplicity
df.columns = ['datetime', 'core_temp', 'bark_temp']
# Convert date and time info to datetime objects
df['datetime'] = pd.to_datetime(df['datetime'], format="%m/%d/%Y, %H:%M:%S %p")
# Convert temperatures to Kelvin
df['core_temp'] = df['core_temp'] + 273.15
df['bark_temp'] = df['bark_temp'] + 273.15


# Separate data from individual days, organized in a dictionary
data_dict = {}
# Here we specify dates to study
daterange = []
years = [2021]
months = [1]
dates = range(9,15)
for year in years:
    for month in months:
        for date in dates:
            daterange.append(datetime(year, month, date))
end_delta = relativedelta.relativedelta(days=1)
for date in daterange:
    end_date = date + end_delta
    mask = (df['datetime'] >= date) & (df['datetime'] < end_date)
    data_on_date = df.loc[mask]
    # Flip data rows so the start of the day is the top of the table
    # data_on_date = data_on_date.sort_values('datetime')
    # Fun fact: sort tool is not sophisticated enough to know that 12am is morning...
    # I'm going to leave it as is unless it creates issues later
    data_dict[date.strftime('%b %d, %Y')]=data_on_date
    print(date.strftime('%b %d, %Y'))
    print(data_dict[date.strftime('%b %d, %Y')])




