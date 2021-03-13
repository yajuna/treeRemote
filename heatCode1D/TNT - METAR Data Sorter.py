# -*- coding: utf-8 -*-
"""
Created on Wed Feb 24 12:15:32 2021

Data retrieved from https://navlost.eu/metar/?icao=SBJP&ts0=2021-01-01&ts1=2021-02-19
Airport code: SBJP (Presidente Castro Airport)

@author: selin
"""

import pandas as pd
from metar import Metar
import pysolar
import pytz
from datetime import *
from dateutil import relativedelta
import json

def parse_metar(file_name):

    with open(file_name) as file:
        # Set up empty data arrays
        time = []
        windspeed = []
        temperature = []
        solar_radiation = []

        # Get the lines from the file
        lines = file.readlines()

        for line in lines:

            line = line.strip()
            info = Metar.Metar(line)

            # Retrieve time
            # I think it makes the data table look cleaner if we don't parse
            # the time measurement here, just store it as a string and parse it
            # in our other code later. But I left the parser code here for reference
            #def dateparse(x): return pd.datetime.strptime(x, '%Y-%m-%d %H:%M:%S')
            time.append(info.time)

            # Retrieve windspeed (m/s)
            w = info.wind("MPS")
            w_numeric = [int(s) for s in w.split() if s.isdigit()]
            wind = 0
            for i in w_numeric:
                wind += i
            windspeed.append(wind)

            # Retrieve temp (K)
            temp = float(info.temp.string('K')[:-2])
            temperature.append(temp)

        d = {'time': time, 'temperature (K)': temperature, 'windspeed (m/s)': windspeed}
        df = pd.DataFrame(data=d)
        
        # Retrieve solar radiation (W/m^3)
        latitude = -7.115
        longitude = -34.863
        zone = pytz.timezone('America/Sao_Paulo')
        for t in time:
            t_timezone_aware = zone.localize(t)
            altitude_deg = pysolar.solar.get_altitude(latitude, longitude, t_timezone_aware)
            radiation = pysolar.radiation.get_radiation_direct(t_timezone_aware, altitude_deg)
            solar_radiation.append(radiation)
        
        df['solar radiation (W/m^3'] = solar_radiation
        return df

def prep_file(inputfile):
    inputfile = 'Joao Pessoa weather report.txt'

    output_lines = []
    with open (inputfile) as file:
        lines = file.readlines()
        for line in lines:
            output_lines.append(line[0:6] + line[36:])
    outputfile = open('Joao Pessoa Weather METAR.txt', 'w')
    outputfile.write(''.join(output_lines))
     
 
def main():
    prep_file('Joao Pessoa weather report.txt')
    df = parse_metar('Joao Pessoa Weather METAR.txt')
    place = 'Joao Pessao'
    month = '2021'
    # Flip data into chronological order
    df = df.iloc[::-1] 
    df.to_csv(r'Weather in ' + place + ' ' + month +'.csv')
 
 
if __name__ == '__main__':
    main()