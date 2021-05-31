# -*- coding: utf-8 -*-
"""
Created on Mon May 31 14:02:29 2021

@author: Mhock
"""


import heatMain as heat

gridPoints = 20
totalPoints = 5
heat.config['gridPoints'] = gridPoints

print(f"Using {gridPoints} gridpoints, the minimum and maximum are: \n")
for i in range(1, totalPoints+1):
    heat.config['at_point'] = int(i*(gridPoints-1)/totalPoints)
    heat.temp(heat.config)