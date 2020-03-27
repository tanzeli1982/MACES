#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 21 06:44:46 2020

Draw the forcing data of 

@author: Zeli Tan
"""

import numpy as np
import matplotlib.pyplot as plt
#from matplotlib.ticker import AutoMinorLocator
from netCDF4 import Dataset

day0 = 456
day1 = 459
# read data
filename = r'/Users/tanz151/Documents/Projects/TAI_BGC/Data/Hydrodynamics_obs/' + \
    r'VeniceLagoon/force_h.nc'
try:
    nc = Dataset(filename,'r')
    h_var = 100 * np.array(nc.variables['h'][day0*144:day1*144+1]) # cm
finally:
    nc.close()
nt_h = np.size(h_var)
tt_h = np.arange(nt_h)
    
filename = r'/Users/tanz151/Documents/Projects/TAI_BGC/Data/Hydrodynamics_obs/' + \
    r'VeniceLagoon/force_U10.nc'
try:
    nc = Dataset(filename,'r')
    U10_var = np.array(nc.variables['U10'][day0*96:day1*96+1])
    U10_dir_var = np.array(nc.variables['U10_dir'][day0*96:day1*96+1])
finally:
    nc.close()
nt_U10 = np.size(U10_var)
tt_U10 = np.arange(nt_U10)
    
# plot
plt.clf()
fig, axes = plt.subplots(3, 1, figsize=(5,8))

plt.style.use('default')

ax = axes[0]
ax.plot(tt_h, h_var, color='black', linestyle='-', linewidth=2, alpha=0.9)
ax.set_xlim(0, nt_h-1)
#ax.set_ylim(-75, 75)
ax.xaxis.set_ticks(np.arange(0,nt_h,144))
#ax.yaxis.set_ticks(np.linspace(-75,75,7))
ax.set_xticklabels(['4/2','4/3','4/4','4/5'])
#ax.xaxis.set_minor_locator(AutoMinorLocator(2))
ylabel = 'Water level ($\mathregular{cm}$ $\mathregular{asl}$)'
ax.set_ylabel(ylabel, fontsize=12, fontname='Times New Roman', color='black')
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(12) for label in labels]
[label.set_color('black') for label in labels]

ax = axes[1]
ax.plot(tt_U10, U10_var, color='black', linestyle='-', linewidth=2, alpha=0.9)
ax.set_xlim(0, nt_U10-1)
#ax.set_ylim(0, 15)
ax.xaxis.set_ticks(np.arange(0,nt_U10,96))
#ax.yaxis.set_ticks(np.linspace(0,15,6))
ax.set_xticklabels(['4/2','4/3','4/4','4/5'])
#ax.xaxis.set_minor_locator(AutoMinorLocator(2))
ylabel = 'Wind speed ($\mathregular{m}$ $\mathregular{{s}^{-1}}$)'
ax.set_ylabel(ylabel, fontsize=12, fontname='Times New Roman', color='black')
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(12) for label in labels]
[label.set_color('black') for label in labels]

ax = axes[2]
ax.plot(tt_U10, U10_dir_var, color='black', linestyle='-', linewidth=2, alpha=0.9)
ax.set_xlim(0, nt_U10-1)
#ax.set_ylim(0, 15)
ax.xaxis.set_ticks(np.arange(0,nt_U10,96))
#ax.yaxis.set_ticks(np.linspace(0,15,6))
ax.set_xticklabels(['4/2','4/3','4/4','4/5'])
#ax.xaxis.set_minor_locator(AutoMinorLocator(2))
ylabel = 'Wind direction (degree)'
ax.set_ylabel(ylabel, fontsize=12, fontname='Times New Roman', color='black')
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(12) for label in labels]
[label.set_color('black') for label in labels]

plt.tight_layout()
fig.savefig('venice_forcing2.png', dpi=300)
#fig.savefig('venice_forcing2.pdf', dpi=600)
plt.show()
