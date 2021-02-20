#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 18 15:24:43 2021

vetorize TNTiteration2.py: input Ta, Tbinit are vectors

Main function to call is BATCH (bark at tree, coefficient of heat)

EX. 
import TNTvec as tntV
h, Tb = tnt2.BATCH()

@author: yajuna
"""
import numpy as np
from functools import partial
from scipy.optimize import broyden1

import csvReader as cR


param = {"Ta": cR.Ta, "Va": 10, "qrads": 650, "Pr": 0.707, "Ka": 26.3e-3, "Kt": 0.11,
         "nu": 15.89e-6, "epsilon": 0.8, "sigma": 5.67e-8, "C": 0.193, "m": 0.618, "rb": 0.2,
         "L": 10, "DeltaT": 2, "DeltaR": 100 / 1000, "timeSteps": 1000}

def F(Tbinit, **param):
    Ta = param["Ta"]
    qrads = param["qrads"]
    Kt = param["Kt"]
    epsilon = param["epsilon"]
    sigma = param["sigma"]
    L = param["L"]
    DeltaT = param["DeltaT"]
    DeltaR = param["DeltaR"]
    
    rb = param["rb"]
    
    Pr = param["Pr"]
    Ka = param["Ka"]
    C = param["C"]
    m = param["m"]
    Va = param["Va"]
    nu = param["nu"]
    
    Re = Va * (2 * rb) / nu
    Nu = C * (Re ** m) * Pr ** (1 / 3)  # Nusselt number
    h = Nu * Ka / (rb * 2)  # W/m2-k heat transfer coefficient
    
    r1 = rb - DeltaR
    Rcond = np.log(rb / r1) / (2 * np.pi * L * Kt)  # K/W
    qcond = DeltaT / Rcond  # W
    qeq = qcond - qrads
    return Ta + qeq / (2 * np.pi * L * h * rb + 
                       2 * L * epsilon * rb * sigma * np.pi * (Ta ** 3 + Ta ** 2 * Tb + Ta * Tb ** 2 + Tb ** 3)) - Tb

def heatTransferCoeff(**param):
    rb = param["rb"]
    
    Pr = param["Pr"]
    Ka = param["Ka"]
    C = param["C"]
    m = param["m"]
    Va = param["Va"]
    nu = param["nu"]
    
    Re = Va * (2 * rb) / nu
    Nu = C * (Re ** m) * Pr ** (1 / 3)  # Nusselt number
    return Nu * Ka / (rb * 2)  # W/m2-k heat transfer coefficient

Tbinit = 22 + 273.2

def BATCH():
    h = heatTransferCoeff(**param)
    Tb = broyden1(partial(F, **param), Tbinit)
    Tbfinal = Tb - 273.2
    return h, Tbfinal
