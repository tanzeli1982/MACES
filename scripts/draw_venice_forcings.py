#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 21 06:44:46 2020

Draw the forcing data of 

@author: Zeli Tan
"""

import numpy as np
from scipy.io import netcdf
import matplotlib.pyplot as plt
#from matplotlib.ticker import AutoMinorLocator

# read data
filename = r'/Users/tanz151/Documents/Projects/TAI_BGC/Data/Hydrodynamics_obs/' + \
    r'VeniceLagoon/force_h.nc'
try:
    nc = netcdf.netcdf_file(filename,'r')
    h_var = 100 * np.array(nc.variables['h'][0:1008]) # cm
finally:
    nc.close()
h_tt = np.arange(1008)
    
filename = r'/Users/tanz151/Documents/Projects/TAI_BGC/Data/Hydrodynamics_obs/' + \
    r'VeniceLagoon/force_U10.nc'
try:
    nc = netcdf.netcdf_file(filename,'r')
    U10_var = np.array(nc.variables['U10'][0:672])
finally:
    nc.close()
U10_tt = np.arange(672)
    
# plot
plt.clf()
fig, axes = plt.subplots(2, 1, figsize=(6,7))

plt.style.use('default')

ax = axes[0]
ax.plot(h_tt, h_var, color='black', linestyle='-', linewidth=2, alpha=0.9)
ax.set_xlim(0, 1009)
ax.set_ylim(-75, 75)
ax.xaxis.set_ticks(np.arange(0,1009,144))
ax.yaxis.set_ticks(np.linspace(-75,75,7))
ax.set_xticklabels(['1/1','1/2','1/3','1/4','1/5','1/6','1/7','1/8'])
#ax.xaxis.set_minor_locator(AutoMinorLocator(2))
ylabel = 'Water level ($\mathregular{cm}$ $\mathregular{asl}$)'
ax.set_ylabel(ylabel, fontsize=12, fontname='Times New Roman', color='black')
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(12) for label in labels]
[label.set_color('black') for label in labels]

ax = axes[1]
ax.plot(U10_tt, U10_var, color='black', linestyle='-', linewidth=2, alpha=0.9)
ax.set_xlim(0, 673)
ax.set_ylim(0, 15)
ax.xaxis.set_ticks(np.arange(0,673,96))
ax.yaxis.set_ticks(np.linspace(0,15,6))
ax.set_xticklabels(['1/1','1/2','1/3','1/4','1/5','1/6','1/7','1/8'])
#ax.xaxis.set_minor_locator(AutoMinorLocator(2))
ylabel = 'Wind speed ($\mathregular{m}$ $\mathregular{{s}^{-1}}$)'
ax.set_ylabel(ylabel, fontsize=12, fontname='Times New Roman', color='black')
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(12) for label in labels]
[label.set_color('black') for label in labels]

plt.tight_layout()
fig.savefig('venice_forcing.png', dpi=300)
#fig.savefig('venice_forcing.pdf', dpi=600)
plt.show()
