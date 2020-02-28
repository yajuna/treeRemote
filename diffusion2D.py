#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 17:11:10 2020

Created for 2D diffusion equation. Equations based on Potter and Andresen (2002)

code ref: https://scipython.com/book/chapter-7-matplotlib/examples/the-two-dimensional-diffusion-equation/

to do: add source terms, add parameters, write in polar coordinates

code written in Python2.7, need to convert to Python 3.*

@author: yajun
"""

import numpy as np
import scipy.sparse as spsp
import matplotlib.pyplot as plt
import matplotlib.cm as cm
############################## Problem Setup ###################################

# Note: When you import functions from this script, code that is inside an
# "if __name__ == '__main__'" block will not be executed.  When you *run* the
# script from an IPython console, code inside these blocks will run.

if __name__ == '__main__':
	# Setup the problem paramters
    config = dict()
	
	# Configure the tree
    config['k'] = 1 # thermal conductivity
    config['c'] = 1 # specific heat unit J*kg^-1*K^-1	
	
    # Physical domain parameters
    config['r_limits'] = [0.0, 20.0]
    config['nr'] = 201
    config['dr'] = (config['r_limits'][1] - config['r_limits'][0]) / (config['nr']-1)
    
    config['phi_limits'] = [0, 2*np.pi]
    config['nphi'] = 201
    config['dphi'] = (config['phi_limits'][1] - config['phi_limits'][0]) / (config['nphi']-1)
	
	# Load the source term
#	from source_tools import source1
#	f = source(config)
	
	# Set CFL safety constant
    config['alpha'] = 1.0/6.0
	
	# Define time step parameters
    config['T'] = 3 #seconds
    config['dt'] = config['alpha'] * config['dx']
    config['nt'] = int(config['T']/config['dt'])




