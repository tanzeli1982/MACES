#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 21 23:57:20 2019

Compare organic matter accretion under different sea-level rise

@author: Zeli Tan
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.io import netcdf
from matplotlib.ticker import AutoMinorLocator

rdir = '/qfs/projects/taim/TAIMOD/tests/'
Lcoast = 113.059584     # coast line length (km)
# read data
filename = rdir + 'T03MOD_NULLMOD/low_rslr/466_ecogeom.nc'
try:
    nc = netcdf.netcdf_file(filename, 'r')
    x = np.array(nc.variables['x'][:])
finally:
    nc.close()
nx = np.size(x)
dx = np.zeros(nx, dtype=np.float64)
for ii in range(nx):
    if ii==0:
        dx[ii] = 0.5*(x[ii+1]-x[ii])
    elif ii==nx-1:
        dx[ii] = 0.5*(x[ii]-x[ii-1])
    else:
        dx[ii] = 0.5*(x[ii+1]-x[ii-1])
models = ['F06MOD', 'T03MOD', 'KM12MOD', 'VDK05MOD']
omac_methods = ['M12MOD', 'DA07MOD', 'KM12MOD']
minvals = [0.0, 0.0, 0.01]
carbon_burial_high = {}   # high sea-level rise
carbon_burial_low = {}    # no sea-level rise
for ii, method in enumerate(omac_methods):
    nmod = len(models)
    carbon_burial_low[method] = np.NaN * np.ones(nmod)
    for jj, model in enumerate(models):
        filename = rdir + model + '_' + method + '/low_rslr/466_ecogeom.nc'
        try:
            nc = netcdf.netcdf_file(filename, 'r')
            depOM = np.sum(np.array(nc.variables['DepOM'][:])[-367:-2], axis=0)
        finally:
            nc.close()
        carbon_burial_low[method][jj] = max(minvals[ii], \
            1e-9 * 8.64e4 * np.sum(depOM * dx * 1e3 * Lcoast))   # Tg
    carbon_burial_high[method] = np.NaN * np.ones(nmod)
    for jj, model in enumerate(models):
        filename = rdir + model + '_' + method + '/high_rslr/466_ecogeom.nc'
        try:
            nc = netcdf.netcdf_file(filename, 'r')
            depOM = np.sum(np.array(nc.variables['DepOM'][:])[-367:-2], axis=0)
        finally:
            nc.close()
        carbon_burial_high[method][jj] = max(minvals[ii], \
            1e-9 * 8.64e4 * np.sum(depOM * dx * 1e3 * Lcoast))   # Tg 

# plot
plt.clf()
fig, axes = plt.subplots(1, 3, figsize=(16,5))

plt.style.use('default')

xm = np.arange(len(models))
width = 0.2 # the width of the bars

# M12MOD
ax = axes[0]
h1 = ax.bar(xm-width/2, 1e3*carbon_burial_low['M12MOD'], width, color='C0')
h2 = ax.bar(xm+width/2, 1e3*carbon_burial_high['M12MOD'], width, color='C3')
handles = [h1, h2]
ax.legend(handles, ['No sea-level rise','Rapid sea-level rise'],
          numpoints=1, loc=0, ncol=1, framealpha=0.0,
          prop={'family':'Times New Roman', 'size':'x-large'})
#ax.set_xlim(-1, len(models))
ax.set_ylim(0, 12)
ax.xaxis.set_ticks(xm)
ax.yaxis.set_ticks(np.arange(0,15,3))
ax.xaxis.set_ticklabels(models)
ax.yaxis.set_minor_locator(AutoMinorLocator(4))
ylabel = 'OM deposition ($\mathregular{Gg}$ $\mathregular{yr^{-1}}$)'
ax.set_ylabel(ylabel, fontsize=14, fontname='Times New Roman', color='black')
ax.text(0.05, 0.95, 'a', transform=ax.transAxes, fontsize=16,
        fontname='Times New Roman', fontweight='bold')
ax.tick_params(which='major', direction='in', colors='xkcd:black', length=8, pad=8)
ax.tick_params(which='minor', direction='in', colors='xkcd:black')

# DA07MOD
ax = axes[1]
h1 = ax.bar(xm-width/2, carbon_burial_low['DA07MOD'], width, color='C0')
h2 = ax.bar(xm+width/2, carbon_burial_high['DA07MOD'], width, color='C3')
handles = [h1, h2]
ax.legend(handles, ['No sea-level rise','Rapid sea-level rise'],
          numpoints=1, loc=0, ncol=1, framealpha=0.0,
          prop={'family':'Times New Roman', 'size':'x-large'})
#ax.set_xlim(-1, len(models))
ax.set_ylim(0, 2)
ax.xaxis.set_ticks(xm)
ax.yaxis.set_ticks(np.arange(0,2.4,0.4))
ax.xaxis.set_ticklabels(models)
ax.yaxis.set_minor_locator(AutoMinorLocator(4))
ylabel = 'OM deposition ($\mathregular{Tg}$ $\mathregular{yr^{-1}}$)'
ax.set_ylabel(ylabel, fontsize=14, fontname='Times New Roman', color='black')
ax.text(0.05, 0.95, 'b', transform=ax.transAxes, fontsize=16,
        fontname='Times New Roman', fontweight='bold')
ax.tick_params(which='major', direction='in', colors='xkcd:black', length=8, pad=8)
ax.tick_params(which='minor', direction='in', colors='xkcd:black')

# KM12MOD
ax = axes[2]
h1 = ax.bar(xm-width/2, carbon_burial_low['KM12MOD'], width, color='C0')
h2 = ax.bar(xm+width/2, carbon_burial_high['KM12MOD'], width, color='C3')
handles = [h1, h2]
ax.legend(handles, ['No sea-level rise','Rapid sea-level rise'],
          numpoints=1, loc=0, ncol=1, framealpha=0.0,
          prop={'family':'Times New Roman', 'size':'x-large'})
#ax.set_xlim(-1, len(models))
ax.set_ylim(0, 2.0)
ax.xaxis.set_ticks(xm)
ax.yaxis.set_ticks(np.arange(0,2.4,0.4))
ax.xaxis.set_ticklabels(models)
ax.yaxis.set_minor_locator(AutoMinorLocator(4))
ylabel = 'OM deposition ($\mathregular{Tg}$ $\mathregular{yr^{-1}}$)'
ax.set_ylabel(ylabel, fontsize=14, fontname='Times New Roman', color='black')
ax.text(0.05, 0.95, 'c', transform=ax.transAxes, fontsize=16,
        fontname='Times New Roman', fontweight='bold')
ax.tick_params(which='major', direction='in', colors='xkcd:black', length=8, pad=8)
ax.tick_params(which='minor', direction='in', colors='xkcd:black')

# change font
for ax in axes.flatten():
    labels = ax.get_xticklabels() + ax.get_yticklabels()
    [label.set_fontname('Times New Roman') for label in labels]
    [label.set_fontsize(12) for label in labels]
    [label.set_color('black') for label in labels]

plt.tight_layout()
fig.savefig('carbon_deposition.png', dpi=300)
#fig.savefig('carbon_deposition.pdf', dpi=600)
plt.show()