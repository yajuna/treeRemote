#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 20 21:09:10 2021

convert temperature difference (core minus bark) to voltage

How to calculate the current (E) using the temperature difference. We will be 
using the Seebeck effect, so the calculation is super simple. E=-S*DeltaT. S 
is the Seebeck coefficient for wood, but the experimental project team already 
found this for us. So S=-0.259 mV/K if the outside temperature is colder than 
the tree (like night), and S=0.207 mV/K if the outside temperature is warmer 
than the tree [Cite Orlando Paper]. The DeltaT will be the difference between 
the core temperature and the bark temperature (in either C or K). The current 
you calculate E, will be in mV, and we should have a different value of E for 
each point in time we do the calculation for. 

run voltage

@author: yajun
"""
import numpy as np
from scipy import sparse
from matplotlib import pyplot as plt

from heatMain import *
config['output'] = 'temp'
config['visualization dimension'] = 2
soln = temp(config)

point_pair = [1,-5]
soln_diff = soln[:,point_pair[0]]-soln[:,point_pair[1]]

import source as stree # source term at bdry 

s1 = 0.207 # outside warmer
s2 = -0.259 # outside colder

E = soln[:,0]
for j in range(E.size):
    if stree.tempAir[j] > soln[j,-1]:
        E[j] = -s1 * soln_diff[j]
    else:
        E[j] = -s2 * soln_diff[j]

print(E.max(), E.min())   

t = np.linspace(0, 24, E.size, endpoint=False) 

plt.plot(t, E, '.r-')
message = f"Voltage generated at {point_pair}"
plt.title(message)
plt.axis([0,24,-0.01,0.01])
plt.xlabel('Time (hrs)')
plt.ylabel('Voltage (mV)')
#plt.savefig('/home/yajun/Documents/treePower/figs/' + 'StempK' + str(at_point) + '.eps', format='eps', dpi=300,bbox_inches='tight')
plt.show()    
    
    























