#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 12 10:55:03 2021

this code computes the source term used for heatMain.py. It is 
to replace sourceS.py and sourceN.py

@author: yajuna
"""

import numpy as np
import csvReader as cR
import TNTvec as tntv

albedo = .3  # "For albedo, we assumed a value of 0.3" - Potter Anderson page 3
sigma = 5.67e-8

tempAir = np.interp(np.linspace(0,24,1000), np.linspace(0,24,cR.temp16np.size), cR.temp16np)


solar = (1 - albedo) * np.interp(np.linspace(0,24,1000), np.linspace(0,24,cR.solar16np.size), cR.solar16np)


blackbody = sigma * (tempAir**4 - tntv.bdryArray**4) # Tair^4 - Tsfc^4


convect = tntv.hArray * (tntv.bdryArray - tempAir) # Tsfc - Tair

sourceTerm = solar + blackbody + convect

## Change source term to negative
sourceTerm = - sourceTerm

## 

barkTemp = np.interp(np.linspace(0,24,1000), np.linspace(0,24,cR.barkTemp16np.size),cR.barkTemp16np)


