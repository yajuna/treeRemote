#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 28 02:35:54 2021

Main function to call is BarkTempAndHeatCoeff

EX. 
import TNTiteration2 as tnt2
h, Tb = tnt2.BarkTempAndHeatCoeff()

@author: yajuna
"""
import numpy as np
from functools import partial
from scipy.optimize import broyden1

param = {"Ta": 29 + 273.2, "Va": 10, "qrads": 650, "Pr": 0.707, "Ka": 26.3e-3, "Kt": 0.11,
         "nu": 15.89e-6, "epsilon": 0.8, "sigma": 5.67e-8, "C": 0.193, "m": 0.618, "rb": 0.2,
         "L": 10, "DeltaT": 2, "DeltaR": 100 / 1000}


def F(Tb, **param):
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

def BarkTempAndHeatCoeff():
    h = heatTransferCoeff(**param)
    Tb = broyden1(partial(F, **param), Tbinit)
    Tbfinal = Tb - 273.2
    return h, Tbfinal

"""
%Tb = unknown - this is what we are looking for
% Variables determined from the weather data look up
Ta = 29+273.2; %K Ambient air temperature
Va = 10; %m/s wind velocity
% Variables determined from the pysolar
qrads = 650; %W/m2 incoming radiation from the sun 
% Constants that we look up
Pr = 0.707; % Prandtl number for air @300K
Ka = 26.3*10^(-3); %W/m-K Conductivity of air @300K
Kt = 0.11; %W/m-K Conductivity estimate of the tree
nu = 15.89*10^(-6); %m2/s viscosity of the air @300K
epsilon = 0.8; %emissivity of the tree trunk
sigma = 5.67*10^(-8); %W/m2-K^4 Stephan Boltzmann constant
% From Table 7.2 in Bergman
C = 0.193;
m = 0.618;
% Values we estimate or match with Brazil experimental data
rb = 0.2; %m radius of the tree at the bark
L = 10; %m height of the tree
%Things we need to estimate
DeltaT = 2; %oC or K. This is the temperature difference between the
%bark of the tree and the sensor located a known distance in. Ideally this
%is from the experimental data in Brazil. I estimated this number from
%Protasio's slides
% We could also estimate this value using the last iteration of the FEA model
% at a known node
r1 = rb-(100/1000); %m radius of the tree at the point the second temperature is measured
% Estimate h using correlations
Re = Va*(2*rb)/nu %Reynolds number
Nu = C*(Re^m)*(Pr^(1/3)); %Nusselt number
h = Nu*Ka/(rb*2) %W/m2-k heat transfer coefficient
% Calculate the conduction term
Rcond = log(rb/r1)/(2*pi*L*Kt); % K/W
qcond = DeltaT/Rcond; %W
%% Start the iteration/root solver
% all this code needs to live inside the root solver, which will be more
% efficient than a loop
Tb = 22+273.2; %initial guess of the temperature

for i=1:10
% Estimate hrad using iteration
hr = epsilon*sigma*(Tb+Ta)*(Tb^2 + Ta^2); %W/m2-K approximation of the radiation term
% Resistors for convection and conduction
Rconv = 1/(h*2*pi*rb*L); %K/W
Rrad = 1/(hr*2*pi*rb*L); %K/W
Req = 1/((1/Rrad)+(1/Rconv)); %K/W
% Conservation of energy at the Tb node
%qcond = qrads + qeq
qeq = qcond-qrads; %W
DeltaTeq = qeq*Req; %K
Tb2=(DeltaTeq + Ta); %K
"""





