#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  4 13:12:39 2020

Plot the transect of the study sites.

@author: Zeli Tan
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from netCDF4 import Dataset
from matplotlib.ticker import AutoMinorLocator

# read the transect of Venice Lagoon
filename = '/Users/tanz151/Python_maces/scripts/maces_ecogeom_2002-01-01_2004-01-01_466.nc'
try:
    nc = Dataset(filename,'r')
    x_venice = np.array(nc.variables['x'][:])
    zh_venice = np.array(nc.variables['zh'][0,:])
finally:
    nc.close()
index = np.argmin(np.abs(zh_venice))
x_venice = x_venice - x_venice[index]
    
# read the transect of Plum Island
filename = '/Users/tanz151/Python_maces/scripts/maces_ecogeom_2017-01-01_2019-01-01_4097.nc'
try:
    nc = Dataset(filename,'r')
    x_plum = np.array(nc.variables['x'][:])
    zh_plum = np.array(nc.variables['zh'][0,:])
finally:
    nc.close()
index = np.argmin(np.abs(zh_plum))
x_plum = x_plum - x_plum[index]
    
# read the transect of Venice Lagoon
filename = '/Users/tanz151/Python_maces/scripts/maces_ecogeom_2004-01-01_2005-01-01_11099.nc'
try:
    nc = Dataset(filename,'r')
    x_hunter = np.array(nc.variables['x'][:])
    zh_hunter = np.array(nc.variables['zh'][0,:])
finally:
    nc.close()
index = np.argmin(np.abs(zh_hunter))
x_hunter = x_hunter - x_hunter[index]
    
plt.clf()
fig = plt.figure(figsize=(8,6.5))

gs = gridspec.GridSpec(nrows=2, ncols=2)

plt.style.use('default')

# Venice Lagoon
ax = fig.add_subplot(gs[0,0])
ax.plot(1e-3*x_venice, zh_venice, color='black', ls='-', lw=3, alpha=1)
ax.plot(1e-3*x_venice, np.zeros_like(x_venice), color='black', ls=':', lw=2, alpha=0.8)
ax.plot([0,0], [-6,2], color='black', ls=':', lw=1, alpha=0.8)
ax.xaxis.set_ticks(np.arange(-20,25,5))
ax.yaxis.set_ticks(np.arange(-6,4,2))
ax.set_xlim(1e-3*x_venice[0], 1e-3*x_venice[-1])
ax.set_ylim(-6, 2)
ax.xaxis.set_minor_locator(AutoMinorLocator(5))
ax.tick_params(which='major', direction='in', length=6)
ax.tick_params(which='minor', direction='in')
ax.set_xlabel('Distance (km)', fontsize=12, fontname='Times New Roman', 
              fontweight='bold')
ax.set_ylabel('Elevation (masl)', fontsize=12, fontname='Times New Roman', 
              fontweight='bold')
ax.set_title('Venice Lagoon', fontsize=14, fontname='Times New Roman', 
             fontweight='bold')
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(12) for label in labels]
[label.set_color('black') for label in labels]
[label.set_fontweight('bold') for label in labels]

# Plum Island Estuary
ax = fig.add_subplot(gs[0,1])
ax.plot(x_plum, zh_plum, color='black', ls='-', lw=3, alpha=1)
ax.plot(x_plum, np.zeros_like(x_plum), color='black', ls=':', lw=2, alpha=0.8)
ax.plot([0,0], [-5,25], color='black', ls=':', lw=1, alpha=0.8)
ax.xaxis.set_ticks(np.arange(-800,1600,400))
ax.yaxis.set_ticks(np.arange(-5,25,5))
ax.set_xlim(x_plum[0], x_plum[-1])
ax.set_ylim(-5, 20)
ax.xaxis.set_minor_locator(AutoMinorLocator(5))
ax.tick_params(which='major', direction='in', length=6)
ax.tick_params(which='minor', direction='in')
ax.set_xlabel('Distance (m)', fontsize=12, fontname='Times New Roman', 
              fontweight='bold')
ax.set_ylabel('Elevation (masl)', fontsize=12, fontname='Times New Roman', 
              fontweight='bold')
ax.set_title('Plum Island Estuary', fontsize=14, fontname='Times New Roman', 
             fontweight='bold')
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(12) for label in labels]
[label.set_color('black') for label in labels]
[label.set_fontweight('bold') for label in labels]

# Hunter Estuary
ax = fig.add_subplot(gs[1,0])
ax.plot(x_hunter, zh_hunter, color='black', ls='-', lw=3, alpha=1)
ax.plot(x_hunter, np.zeros_like(x_hunter), color='black', ls=':', lw=2, alpha=0.8)
ax.plot([0,0], [-2,3], color='black', ls=':', lw=1, alpha=0.8)
ax.xaxis.set_ticks(np.arange(-25,200,25))
ax.yaxis.set_ticks(np.arange(-2,4,1))
ax.set_xlim(x_hunter[0], x_hunter[-1])
ax.set_ylim(-2, 3)
ax.xaxis.set_minor_locator(AutoMinorLocator(5))
ax.tick_params(which='major', direction='in', length=6)
ax.tick_params(which='minor', direction='in')
ax.set_xlabel('Distance (m)', fontsize=12, fontname='Times New Roman', 
              fontweight='bold')
ax.set_ylabel('Elevation (mAHD)', fontsize=12, fontname='Times New Roman', 
              fontweight='bold')
ax.set_title('Hunter Estuary', fontsize=14, fontname='Times New Roman', 
             fontweight='bold')
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(12) for label in labels]
[label.set_color('black') for label in labels]
[label.set_fontweight('bold') for label in labels]

plt.tight_layout()
fig.savefig('F3.png', dpi=300)
fig.savefig('F3.jpg', dpi=600)
plt.show()