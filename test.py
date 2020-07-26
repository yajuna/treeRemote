#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 25 17:25:44 2020

@author: selinateng
"""

from pysolar.solar import *
from datetime import datetime, timedelta

def main():
    print('hello')
    
    
    test_pysolar()
    
    # test datetime_range
    dts = [dt.strftime('%Y-%m-%d T%H:%M Z') for dt in 
           datetime_range(datetime(2016, 9, 1, 7), datetime(2016, 9, 1, 9+12), 
           timedelta(minutes=15))]
    print(dts)
    
    # test pysolar getting radiation with a time interval using datetime
    my_test(dts)
    
def test_pysolar():
    latitude_deg = 42.206 # positive in the northern hemisphere
    longitude_deg = -71.382 # negative reckoning west from prime meridian in Greenwich, England
    date = dt.datetime(2007, 2, 18, 15, 13, 1, 130320, tzinfo=dt.timezone.utc)
    altitude_deg = get_altitude(latitude_deg, longitude_deg, date)
    print(radiation.get_radiation_direct(date, altitude_deg))

def datetime_range(start, end, delta):
    current = start
    while current < end:
        yield current
        current += delta

def my_test(dts):
    latitude_deg = 42.206 # positive in the northern hemisphere
    longitude_deg = -71.382 # negative reckoning west from prime meridian in Greenwich, England
    for date in dts:
        print(type(date))
        #altitude_deg = get_altitude(latitude_deg, longitude_deg, date)
        #print(radiation.get_radiation_direct(date, altitude_deg))
    
    
if __name__ == '__main__':
    main()
