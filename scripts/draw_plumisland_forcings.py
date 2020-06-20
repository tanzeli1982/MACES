#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 21 06:44:46 2020

Draw the forcing data of Plum Island

@author: Zeli Tan
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from netCDF4 import Dataset
from datetime import date

t0s = []
t1s = []
t0s.append( (date(2017,7,17) - date(2012,1,1)).days )
t1s.append( (date(2017,7,24) - date(2012,1,1)).days )
t0s.append( (date(2017,10,6) - date(2012,1,1)).days )
t1s.append( (date(2017,10,13) - date(2012,1,1)).days )

h_vars = []
nt_h = []
# read data
filename = r'/Users/tanz151/Documents/Projects/TAI_BGC/Data/Hydrodynamics_obs/' + \
    r'PlumIsland/force_h.nc'
try:
    nc = Dataset(filename,'r')
    for ii, t0 in enumerate(t0s):
        t1 = t1s[ii]
        h_var = 100 * np.array(nc.variables['h'][96*t0:96*t1,0])
        h_vars.append(h_var)
        nt_h.append(len(h_var))
finally:
    nc.close()
    
U10_vars = []
nt_U10 = []
filename = r'/Users/tanz151/Documents/Projects/TAI_BGC/Data/Hydrodynamics_obs/' + \
    r'PlumIsland/force_U10.nc'
try:
    nc = Dataset(filename,'r')
    for ii, t0 in enumerate(t0s):
        t1 = t1s[ii]
        U10_var = np.array(nc.variables['U10'][96*t0:96*t1,0])
        U10_vars.append(U10_var)
        nt_U10.append(len(U10_var))
finally:
    nc.close()
    
# plot
plt.clf()
fig, axes = plt.subplots(2, 2, figsize=(8,6))

plt.style.use('default')

ax = axes[0][0]
ax.plot(np.arange(nt_h[0]), h_vars[0], color='black', linestyle='-', 
        linewidth=2, alpha=0.9)
ax.set_xlim(0, nt_h[0])
ax.set_ylim(-200, 200)
ax.xaxis.set_ticks(np.arange(0,nt_h[0]+1,96))
ax.yaxis.set_ticks(np.linspace(-200,200,5))
ax.set_xticklabels(['7/17','7/18','7/19','7/20','7/21','7/22','7/23','7/24'])
ax.xaxis.set_minor_locator(AutoMinorLocator(4))
#ax.yaxis.set_minor_locator(AutoMinorLocator(2))
ylabel = 'Water level ($\mathregular{cm}$ $\mathregular{asl}$)'
ax.set_ylabel(ylabel, fontsize=12, fontname='Times New Roman', color='black')
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(12) for label in labels]
[label.set_color('black') for label in labels]
ax.tick_params(which='major', direction='in', colors='xkcd:black', length=6, pad=8)
ax.tick_params(which='minor', direction='in', colors='xkcd:black')

ax = axes[0][1]
ax.plot(np.arange(nt_h[1]), h_vars[1], color='black', linestyle='-', 
        linewidth=2, alpha=0.9)
ax.set_xlim(0, nt_h[1])
ax.set_ylim(-200, 200)
ax.xaxis.set_ticks(np.arange(0,nt_h[1]+1,96))
ax.yaxis.set_ticks(np.linspace(-200,200,5))
ax.set_xticklabels(['10/6','10/7','10/8','10/9','10/10','10/11','10/12','10/13'])
ax.xaxis.set_minor_locator(AutoMinorLocator(4))
#ax.yaxis.set_minor_locator(AutoMinorLocator(2))
#ax.set_xlabel('Time', fontsize=12, fontname='Times New Roman', color='black')
#ylabel = 'Water level ($\mathregular{cm}$ $\mathregular{AHD}$)'
#ax.set_ylabel(ylabel, fontsize=12, fontname='Times New Roman', color='black')
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(12) for label in labels]
[label.set_color('black') for label in labels]
ax.tick_params(which='major', direction='in', colors='xkcd:black', length=6, pad=8)
ax.tick_params(which='minor', direction='in', colors='xkcd:black')

ax = axes[1][0]
ax.plot(np.arange(nt_U10[0]), U10_vars[0], color='black', linestyle='-', 
        linewidth=2, alpha=0.9)
ax.set_xlim(0, nt_U10[0])
ax.set_ylim(0, 9)
ax.xaxis.set_ticks(np.arange(0,nt_U10[0]+1,96))
ax.yaxis.set_ticks(np.linspace(0,8,5))
ax.set_xticklabels(['7/17','7/18','7/19','7/20','7/21','7/22','7/23','7/24'])
ax.xaxis.set_minor_locator(AutoMinorLocator(4))
#ax.yaxis.set_minor_locator(AutoMinorLocator(2))
ax.set_xlabel('Time', fontsize=12, fontname='Times New Roman', color='black')
ylabel = 'Wind speed ($\mathregular{m}$ $\mathregular{{s}^{-1}}$)'
ax.set_ylabel(ylabel, fontsize=12, fontname='Times New Roman', color='black')
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(12) for label in labels]
[label.set_color('black') for label in labels]
ax.tick_params(which='major', direction='in', colors='xkcd:black', length=6, pad=8)
ax.tick_params(which='minor', direction='in', colors='xkcd:black')

ax = axes[1][1]
ax.plot(np.arange(nt_U10[1]), U10_vars[1], color='black', linestyle='-', 
        linewidth=2, alpha=0.9)
ax.set_xlim(0, nt_U10[1])
ax.set_ylim(0, 9)
ax.xaxis.set_ticks(np.arange(0,nt_U10[1]+1,96))
ax.yaxis.set_ticks(np.linspace(0,8,5))
ax.set_xticklabels(['10/6','10/7','10/8','10/9','10/10','10/11','10/12','10/13'])
ax.xaxis.set_minor_locator(AutoMinorLocator(4))
#ax.yaxis.set_minor_locator(AutoMinorLocator(2))
ax.set_xlabel('Time', fontsize=12, fontname='Times New Roman', color='black')
#ylabel = 'Wind speed ($\mathregular{m}$ $\mathregular{{s}^{-1}}$)'
#ax.set_ylabel(ylabel, fontsize=12, fontname='Times New Roman', color='black')
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(12) for label in labels]
[label.set_color('black') for label in labels]
ax.tick_params(which='major', direction='in', colors='xkcd:black', length=6, pad=8)
ax.tick_params(which='minor', direction='in', colors='xkcd:black')

plt.tight_layout()
fig.savefig('F7.png', dpi=300)
fig.savefig('F7.pdf', dpi=600)
plt.show()
