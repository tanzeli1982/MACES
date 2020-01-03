#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  3 14:03:30 2020

Draw the dynamics of water flow velocity on the TAI platform

@author: Zeli Tan
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.io import netcdf

# read the simulated water flow velocity
filename = '/Users/tanz151/Python_maces/src/out_hydro_2001-01-01_2001-01-06.nc'
try:
    nc = netcdf.netcdf_file(filename,'r')
    xv = np.array(nc.variables['x'][:][0])
    Uw = np.array(nc.variables['U'][:][0])
finally:
    nc.close()
    
# plot
plt.clf()
fig, ax = plt.subplots(figsize=(8,10))

plt.style.use('default')

