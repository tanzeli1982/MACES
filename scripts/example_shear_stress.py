#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  7 16:25:48 2020

Draw the dynamics of bottom shear stress on the TAI platform

@author: Zeli Tan
"""

import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from matplotlib.ticker import AutoMinorLocator

# read the platform elevation
filename = '/Users/tanz151/Python_maces/src/out_ecogeom_2001-01-01_2001-01-06.nc'
try:
    nc = Dataset(filename,'r')
    zh = np.array(nc.variables['zh'][:][0][0])
finally:
    nc.close()
    
# read the simulated water flow velocity
filename = '/Users/tanz151/Python_maces/src/out_hydro_2001-01-01_2001-01-06.nc'
try:
    nc = Dataset(filename,'r')
    xv = np.array(nc.variables['x'][:][0])
    tau = np.array(nc.variables['tau'][:][0])
finally:
    nc.close()
xv = xv[100:308]
tau = tau[:,100:308]
zh = zh[100:308]
indx = np.argmin(np.abs(zh))

xv = xv - xv[indx]

ntime = np.shape(tau)[0]
tt = -np.arange(ntime)

# plot
plt.clf()
fig, ax = plt.subplots(figsize=(7.5,9.5))

plt.style.use('default')

cf = ax.contourf(xv, tt, tau, 10, cmap='hot_r')
ax.plot([0,0], [tt[-1],tt[0]], color='black', ls='--', lw=1, 
        alpha=0.8)
ax.set_xlim([xv[0], xv[-1]])
ax.set_ylim([tt[-1], tt[0]])
ax.xaxis.set_ticks(np.arange(-7500,12500,2500))
ax.yaxis.set_ticks(np.arange(-100,20,20))
ax.yaxis.set_ticklabels(['100','80','60','40','20','0'])
ax.xaxis.set_minor_locator(AutoMinorLocator(5))
ax.yaxis.set_minor_locator(AutoMinorLocator(5))
ax.set_xlabel('Coordinate ($\mathregular{m}$)', fontsize=12, 
              fontname='Times New Roman', color='black')
ax.set_ylabel('Time (hours)', fontsize=12, 
              fontname='Times New Roman', color='black')
ax.tick_params(which='major', direction='in', length=6)
ax.tick_params(which='minor', direction='in', length=2)
cbar = fig.colorbar(cf, ax=ax, orientation='horizontal', pad=0.07)
ylabel = 'Bottom shear stress ($\mathregular{Pa}$)'
cbar.set_label(ylabel, fontsize=12, fontname='Times New Roman', 
               labelpad=0)
#cbar.ax.xaxis.set_minor_locator(AutoMinorLocator(2))
labels = cbar.ax.get_xticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(12) for label in labels]

labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(12) for label in labels]

plt.tight_layout()
fig.savefig('tau.png', dpi=300)
#fig.savefig('F1.pdf', dpi=600)
plt.show()