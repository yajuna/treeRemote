#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 26 11:24:33 2020

@author: selinateng


How to find i at a given time and angle (for all t and phi)?
Currently i is a filler constant 2pi

Potter & Andresen variables - 
    i       % angle of incidence = angle between direction of sunlight and local normal to tree surface
    sun     % angle of sunlight
    norm    % local normal to tree surface
    i = abs(norm - sun)
    %The local normal to the tree's surface can just be the angle phi
    i = abs(phi - sun)

The angle of sunlight would realistically depend on surrounding shade (other trees/objects that block out sun)
But maybe for simplicity we can treat it as the same as the azimuth angle
The azimuth angle is a horizontal coordinate, it describes the relative direction of the sun along the local horizon
There are multiple sign conventions so we need to be careful
Source: https://en.wikipedia.org/wiki/Solar_azimuth_angle

Can we use a mathematical function for azimuth angle or do we need data?
Do we already have some code that would do the job? -> No

Could code it but it's complicated? Needs equation of time to be accurate over extended periods.
Sun position (azimuth) equations:
    https://www.pveducation.org/pvcdrom/properties-of-sunlight/the-suns-position

Here is a polar path calculator that gives sun position for a certain day: 
    https://www2.pvlighthouse.com.au/calculators/solar%20path%20calculator/solar%20path%20calculator.aspx
    
*****Found an existing Python script for computing sun position on Github.*****
    https://github.com/s-bear/sun-position
    Use notes:
        Sign convention: the azimuth angle is 0 at North and positive towards the east
        Do we have a sign convention for phi? (-> no?)
        Input is timestamp, latitude, longitude, elevation (also more optional inputs available)
    Sample usage:
        now = datetime.utcnow()
        az,zen = sunpos(now,LAT,LON,0)
        
    So we also need to call datetime (a built-in Python function)
        Reference: https://docs.python.org/3/library/datetime.html
        Sample usage of sunpos calls datetime.utcnow(), Python says this datetime.now(tz=None) is preferred
        How do we get it to use timestamp data from our measurements and not the current time?
        
        
*****Found a collection of Python libraries for computing sun position + solar radiation*****
https://pysolar.readthedocs.io/en/latest/
We could get both the azimuth and direct solar radiation through this
Consider using Pysolar in replacement of the equation provided in Potter & Andresen for finding solar radiation? 

Not sure where to go from here.

Questions
    Should we use sunpos or pysolar as shortcut?
    Do we need azimuth angle for anything else?
    What is our sign convention for phi (where is it 0?)
    Need help understanding pseudocode to find i for all t and phi in heatPolar.py (unsure of how to find t and phi,
    and what units they are in, in the stepping method)
    For next week - should I work on implementing one of these options for getting azimuth angle so we can finish
    coding for source term i, or work on something else?

"""
    