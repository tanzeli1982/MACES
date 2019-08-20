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

rdir = '/Users/tanz151/Downloads/tests/'
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
carbon_burial_diff = {}   # high sea-level rise - no sea-level rise
for ii, method in enumerate(omac_methods):
    nmod = len(models)
    carbon_burial_low = np.NaN * np.ones(nmod)
    for jj, model in enumerate(models):
        filename = rdir + model + '_' + method + '/low_rslr/466_ecogeom.nc'
        try:
            nc = netcdf.netcdf_file(filename, 'r')
            depOM = np.sum(np.array(nc.variables['DepOM'][:])[-367:-2], axis=0)
        finally:
            nc.close()
        carbon_burial_low[jj] = max(0.0, 1e-9 * 8.64e4 * np.sum(depOM * dx * 1e3 * Lcoast))   # Tg
    carbon_burial_high = np.NaN * np.ones(nmod)
    for jj, model in enumerate(models):
        filename = rdir + model + '_' + method + '/high_rslr/466_ecogeom.nc'
        try:
            nc = netcdf.netcdf_file(filename, 'r')
            depOM = np.sum(np.array(nc.variables['DepOM'][:])[-367:-2], axis=0)
        finally:
            nc.close()
        carbon_burial_high[jj] = max(0.0, 1e-9 * 8.64e4 * np.sum(depOM * dx * 1e3 * Lcoast))   # Tg 
    carbon_burial_diff[method] = carbon_burial_high - carbon_burial_low

# plot
plt.clf()
fig, ax = plt.subplots(figsize=(7.5,5))

plt.style.use('default')

xm = np.arange(len(models))
width = 0.2 # the width of the bars

h1 = ax.bar(xm-width, 10*carbon_burial_diff['M12MOD'], width, color='C0')
h2 = ax.bar(xm, carbon_burial_diff['DA07MOD'], width, color='C4')
h3 = ax.bar(xm+width, carbon_burial_diff['KM12MOD'], width, color='C3')
handles = [h1, h2, h3]
ax.legend(handles, ['OM1','OM2', 'OM3'],
          numpoints=1, loc=0, ncol=1, framealpha=0.0,
          prop={'family':'Times New Roman', 'size':'x-large'})
#ax.set_xlim(-1, len(models))
ax.set_ylim(-0.1, 2.5)
ax.xaxis.set_ticks(xm)
ax.yaxis.set_ticks([0,0.5,1,1.5,2,2.5])
ax.xaxis.set_ticklabels(models)
ax.yaxis.set_minor_locator(AutoMinorLocator(5))
ylabel = r'$\mathregular{\Delta}$OM deposition ($\mathregular{Tg}$ $\mathregular{yr^{-1}}$)'
ax.set_ylabel(ylabel, fontsize=14, fontname='Times New Roman', color='black')
ax.tick_params(which='major', direction='in', colors='xkcd:black', length=8, pad=8)
ax.tick_params(which='minor', direction='in', colors='xkcd:black')

# change font
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(12) for label in labels]
[label.set_color('black') for label in labels]

plt.tight_layout()
fig.savefig('carbon_deposition.png', dpi=300)
#fig.savefig('carbon_deposition.pdf', dpi=600)
plt.show()