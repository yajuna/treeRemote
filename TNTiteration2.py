#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan  9 15:15:18 2021

Translation of Heather's TNTiteration file to Python. Original Matlab code attached in the end

To do:
    1. incorporate nonlinear solver from Michael and Selina to replace the for loop
    2. incorporate the solver into main 1D heat solver (heat1dK.py). The Tb would replace 
    g1(t)

Tb = K unknown tree bark temperature
Ta = K ambient air temperature
Va = m/s wind velocity
qrads = W/m2 incoming radiation from the sun
rb = %m radius of the tree at the bark
L = %m height of the tree
DeltaT = 2; %oC or K. This is the temperature difference between the
bark of the tree and the sensor located a known distance in. Ideally this
is from the experimental data in Brazil. Heather estimated this number from Protasio's slides
Pr = 0.707; % Prandtl number for air @300K
Ka = 26.3*10^(-3); %W/m-K Conductivity of air @300K
Kt = 0.11; %W/m-K Conductivity estimate of the tree
nu = 15.89*10^(-6); %m2/s viscosity of the air @300K
epsilon = 0.8; %emissivity of the tree trunk
sigma = 5.67*10^(-8); %W/m2-K^4 Stephan Boltzmann constant
% From Table 7.2 in Bergman
C = 0.193;
m = 0.618; 
rb = m radius of the tree at the bark
L = m height of the three
DeltaT = C or K This is the temperature difference between the bark of the tree and the sensor located a known distance in. Ideally this
is from the experimental data in Brazil. 
We could also estimate this value using the last iteration of the FEA model at a known node

@author: yajun
"""
import numpy as np
import scipy.optimize

def F(Tb):
    epsilon = 0.8
    sigma = 5.67*10**-8
    Ta = 29+273.2
    rb = 0.2
    L = 10
    Kt = 0.11
    DeltaT = 100/1000
    r1 = rb - 100 / 1000
    Rcond = np.log(rb / r1) / (2 * np.pi * L * Kt)  # K/W
    qcond = DeltaT / Rcond  # W
    qrads = 650
    qeq = qcond-qrads
    return qeq/(2*L*epsilon*pi*rb*sigma*Ta**3 + 2*L*epsilon*pi*rb*sigma*Ta**2*Tb + 2*L*epsilon*pi*rb*sigma*Ta*Tb**2 + 2*L*epsilon*pi*rb*sigma*Tb**3 + 2*L*hr*pi*rb) - Tb


param = {"Ta": 29 + 273.2, "Tb": 22 + 273.2, "Va": 10, "qrads": 650, "Pr": 0.707, "Ka": 26.3*10**-3, "Kt": 0.11, 
         "nu": 15.89*10**-6, "epsilon": 0.8, "sigma": 5.67*10**-8, "C": 0.193, "m": 0.618, "rb": 0.2,
         "L": 10, "DeltaT": 2, "DeltaR": 100 / 1000}
def outerBdry(**param):
    Ta = param["Ta"]
    Tb = param["Tb"]
    Va = param["Va"]
    qrads = param["qrads"]
    Pr = param["Pr"]
    Ka = param["Ka"]
    Kt = param["Kt"]
    nu = param["nu"]
    epsilon = param["epsilon"]
    sigma = param["sigma"]
    C = param["C"]
    m = param["m"]
    rb = param["rb"]
    L = param["L"]
    DeltaT = param["DeltaT"]
    DeltaR = param["DeltaR"]

# r1 = m radius of tree at which 2nd temp is measured    
    r1 = rb - DeltaR

#%% estimate h using correlations
    Re = Va * (2 * rb)/nu # Reynolds number
    Nu = C * (Re**m) * Pr**(1/3) # Nusselt number
    h = Nu * Ka / (rb * 2) # W/m2-k heat transfer coefficient
#%% calculate the conduction term
    Rcond = np.log(rb/r1) / (2 * np.pi * L * Kt) # K/W
    qcond = DeltaT / Rcond # W
    qeq = qcond - qrads
#%% iteration step. Replace with root solver
    Tbfinal = scipy.optimize.root(F, Tb)
    hr = epsilon * sigma * (Tb + Ta) * (Tb**2 + Ta**2)
    print("hr = ")
    print(hr)
    return Tbfinal


sltn = outerBdry(param)





"""
% Simple script to test the idea of iteration as a method to solve for
% missing tree bark temperature.
clear all
close all
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
C=0.193;
m=0.618;
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
% this new value of the temperature is used for the calculations
Tb = Tb2;
end
Tbfinal=Tb-273.2 %convert back to C
%% code output to check your numbers
% % h =
%
%
24.6058
%
%
% Tbfinal =
%
%
27.3045
"""