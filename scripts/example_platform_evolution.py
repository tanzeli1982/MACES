#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 21 22:14:31 2019

Compare the evolution of platform elevation by different mineral accretion
models and sea-level rise scenarios

@author: Zeli Tan
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.io import netcdf
from matplotlib.ticker import AutoMinorLocator

rdir = '/qfs/projects/taim/TAIMOD/tests/'
# read data
models = ['F06MOD_NULLMOD', 'T03MOD_NULLMOD', 'KM12MOD_NULLMOD', 
          'VDK05MOD_NULLMOD']
filename = rdir + 'T03MOD_NULLMOD/low_rslr/466_ecogeom.nc'
try:
    nc = netcdf.netcdf_file(filename, 'r')
    x = 1e-3 * np.array(nc.variables['x'][:])   # km
    zh_yr0 = np.array(nc.variables['zh'][:])[0]
finally:
    nc.close()
platform_zh_high = {}   # high sea-level rise
platform_zh_low = {}    # no sea-level rise
for model in models:
    # no sea-level rise
    filename = rdir + model + '/low_rslr/466_ecogeom.nc'
    try:
        nc = netcdf.netcdf_file(filename, 'r')
        zh_yr20 = np.array(nc.variables['zh'][:])[-2]
    finally:
        nc.close()
    platform_zh_low[model] = zh_yr20
    # rapid sea-level rise
    filename = rdir + model + '/high_rslr/466_ecogeom.nc'
    try:
        nc = netcdf.netcdf_file(filename, 'r')
        zh_yr20 = np.array(nc.variables['zh'][:])[-2]
    finally:
        nc.close()
    platform_zh_high[model] = zh_yr20
    
# plot
plt.clf()
fig, axes = plt.subplots(1, 2, figsize=(8,3))

plt.style.use('default')

linestyles = ['-','--',':','-.']
colors = ['C0','C1','C2','C3']

# no sea-level rise
ax = axes[0]
handles = []
h, = ax.plot(x, zh_yr0, color='black', ls='-', lw=1.5, alpha=0.8)
handles.append(h)
for ii, key in enumerate(models):
    h, = ax.plot(x, platform_zh_low[key], color=colors[ii], ls=linestyles[ii], 
                 lw=1, alpha=0.8)
    handles.append(h)
ax.legend(handles, ['initial platform','F06MOD','T03MOD','KM12MOD','VDK05MOD'],
          numpoints=1, loc=0, ncol=1, framealpha=0.0,
          prop={'family':'Times New Roman', 'size':'medium'})
ax.set_xlim(15, 30)
ax.set_ylim(-1, 1)
ax.xaxis.set_ticks(np.arange(15,35,5))
ax.yaxis.set_ticks(np.arange(-1,1.5,0.5))
ax.xaxis.set_minor_locator(AutoMinorLocator(5))
ax.yaxis.set_minor_locator(AutoMinorLocator(4))
ax.set_xlabel('Coordinate ($\mathregular{km}$)', fontsize=12, 
              fontname='Times New Roman', color='black')
ax.set_ylabel('Elevation ($\mathregular{msl}$)', fontsize=12, 
              fontname='Times New Roman', color='black')
ax.text(0.05, 0.9, 'a', transform=ax.transAxes, fontsize=16,
        fontname='Times New Roman', fontweight='bold')

ax.tick_params(which='major', direction='in', colors='xkcd:black', length=6, pad=8)
ax.tick_params(which='minor', direction='in', colors='xkcd:black')

# high sea-level rise
ax = axes[1]
handles = []
h, = ax.plot(x, zh_yr0, color='black', ls='-', lw=1.5, alpha=0.8)
handles.append(h)
for ii, key in enumerate(models):
    h, = ax.plot(x, platform_zh_high[key], color=colors[ii], ls=linestyles[ii], 
                 lw=1, alpha=0.8)
    handles.append(h)
ax.legend(handles, ['initial platform','F06MOD','T03MOD','KM12MOD','VDK05MOD'],
          numpoints=1, loc=0, ncol=1, framealpha=0.0,
          prop={'family':'Times New Roman', 'size':'medium'})
ax.set_xlim(15, 30)
ax.set_ylim(-1, 1)
ax.xaxis.set_ticks(np.arange(15,35,5))
ax.yaxis.set_ticks(np.arange(-1,1.5,0.5))
ax.xaxis.set_minor_locator(AutoMinorLocator(5))
ax.yaxis.set_minor_locator(AutoMinorLocator(4))
ax.set_xlabel('Coordinate ($\mathregular{km}$)', fontsize=12, 
              fontname='Times New Roman', color='black')
ax.set_ylabel('Elevation ($\mathregular{msl}$)', fontsize=12, 
              fontname='Times New Roman', color='black')
ax.text(0.05, 0.9, 'b', transform=ax.transAxes, fontsize=16,
        fontname='Times New Roman', fontweight='bold')

ax.tick_params(which='major', direction='in', colors='xkcd:black', length=6, pad=8)
ax.tick_params(which='minor', direction='in', colors='xkcd:black')

# change font
for ax in axes.flatten():
    labels = ax.get_xticklabels() + ax.get_yticklabels()
    [label.set_fontname('Times New Roman') for label in labels]
    [label.set_fontsize(12) for label in labels]
    [label.set_color('black') for label in labels]

plt.tight_layout()
fig.savefig('platform_evolution.png', dpi=300)
#fig.savefig('platform_evolution.pdf', dpi=600)
plt.show()