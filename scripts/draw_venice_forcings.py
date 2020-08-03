#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 21 06:44:46 2020

Draw the forcing data of Venice Lagoon

@author: Zeli Tan
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from netCDF4 import Dataset

t1 = {'d0': 343, 'd1': 345}
t2 = {'d0': 456, 'd1': 459}
# read data
filename = r'/Users/tanz151/Documents/Projects/TAI_BGC/Data/Hydrodynamics_obs/' + \
    r'VeniceLagoon/force_h.nc'
try:
    nc = Dataset(filename,'r')
    h_var1 = 100 * np.array(nc.variables['h'][t1['d0']*144:t1['d1']*144+1]) # cm
    h_var2 = 100 * np.array(nc.variables['h'][t2['d0']*144:t2['d1']*144+1]) # cm
finally:
    nc.close()
nt1_h = np.size(h_var1)
tt1_h = np.arange(nt1_h)
nt2_h = np.size(h_var2)
tt2_h = np.arange(nt2_h)
print(np.max(h_var1), np.max(h_var2))
    
filename = r'/Users/tanz151/Documents/Projects/TAI_BGC/Data/Hydrodynamics_obs/' + \
    r'VeniceLagoon/force_U10.nc'
try:
    nc = Dataset(filename,'r')
    U10_var1 = np.array(nc.variables['U10'][t1['d0']*96:t1['d1']*96+1])
    U10_var2 = np.array(nc.variables['U10'][t2['d0']*96:t2['d1']*96+1])
    U10_dir_var1 = np.array(nc.variables['U10_dir'][t1['d0']*96:t1['d1']*96+1])
    U10_dir_var2 = np.array(nc.variables['U10_dir'][t1['d0']*96:t1['d1']*96+1])
finally:
    nc.close()
nt1_U10 = np.size(U10_var1)
tt1_U10 = np.arange(nt1_U10)
nt2_U10 = np.size(U10_var2)
tt2_U10 = np.arange(nt2_U10)
print(np.max(U10_var1),np.max(U10_var2))
    
# plot
plt.clf()
fig, axes = plt.subplots(2, 2, figsize=(8,6))

plt.style.use('default')

ax = axes[0][0]
ax.plot(tt1_h, h_var1, color='black', linestyle='-', linewidth=2, alpha=0.9)
ax.set_xlim(0, nt1_h-1)
ax.set_ylim(-52, 65)
ax.xaxis.set_ticks(np.arange(0,nt1_h,144))
ax.yaxis.set_ticks(np.linspace(-40,60,6))
ax.set_xticklabels(['12/10','12/11','12/12'])
ax.xaxis.set_minor_locator(AutoMinorLocator(4))
ylabel = 'Water level ($\mathregular{cm}$ $\mathregular{asl}$)'
ax.set_ylabel(ylabel, fontsize=12, fontname='Times New Roman', color='black')
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(12) for label in labels]
[label.set_color('black') for label in labels]
ax.tick_params(which='major', direction='in', colors='xkcd:black', length=6, pad=8)
ax.tick_params(which='minor', direction='in', colors='xkcd:black')

ax = axes[1][0]
ax.plot(tt1_U10, U10_var1, color='black', linestyle='-', linewidth=2, alpha=0.9)
ax.set_xlim(0, nt1_U10-1)
ax.set_ylim(1, 18)
ax.xaxis.set_ticks(np.arange(0,nt1_U10,96))
ax.yaxis.set_ticks(np.linspace(3,18,6))
ax.set_xticklabels(['12/10','12/11','12/12'])
ax.xaxis.set_minor_locator(AutoMinorLocator(4))
ax.set_xlabel('Time', fontsize=12, fontname='Times New Roman', color='black')
ylabel = 'Wind speed ($\mathregular{m}$ $\mathregular{{s}^{-1}}$)'
ax.set_ylabel(ylabel, fontsize=12, fontname='Times New Roman', color='black')
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(12) for label in labels]
[label.set_color('black') for label in labels]
ax.tick_params(which='major', direction='in', colors='xkcd:black', length=6, pad=8)
ax.tick_params(which='minor', direction='in', colors='xkcd:black')

ax = axes[0][1]
ax.plot(tt2_h, h_var2, color='black', linestyle='-', linewidth=2, alpha=0.9)
ax.set_xlim(0, nt1_h-1)
ax.set_ylim(-52, 65)
ax.xaxis.set_ticks(np.arange(0,nt2_h,144))
ax.yaxis.set_ticks(np.linspace(-40,60,6))
ax.set_xticklabels(['4/2','4/3','4/4','4/5'])
ax.xaxis.set_minor_locator(AutoMinorLocator(4))
#ylabel = 'Water level ($\mathregular{cm}$ $\mathregular{asl}$)'
#ax.set_ylabel(ylabel, fontsize=12, fontname='Times New Roman', color='black')
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(12) for label in labels]
[label.set_color('black') for label in labels]
ax.tick_params(which='major', direction='in', colors='xkcd:black', length=6, pad=8)
ax.tick_params(which='minor', direction='in', colors='xkcd:black')

ax = axes[1][1]
ax.plot(tt2_U10, U10_var2, color='black', linestyle='-', linewidth=2, alpha=0.9)
ax.set_xlim(0, nt2_U10-1)
ax.set_ylim(1, 18)
ax.xaxis.set_ticks(np.arange(0,nt2_U10,96))
ax.yaxis.set_ticks(np.linspace(3,18,6))
ax.set_xticklabels(['4/2','4/3','4/4','4/5'])
ax.xaxis.set_minor_locator(AutoMinorLocator(4))
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
fig.savefig('F3.png', dpi=300)
fig.savefig('F3.pdf', dpi=600)
plt.show()
