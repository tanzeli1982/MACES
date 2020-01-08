#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  7 15:38:28 2020

Draw the dynamics of significant wave height on the TAI platform

@author: Zeli Tan
"""

import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset

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
    Hwav = np.array(nc.variables['Hwav'][:][0])
finally:
    nc.close()
xv = xv[100:308]
Hwav = Hwav[:,100:308]
zh = zh[100:308]
indx = np.argmin(np.abs(zh))

ntime = np.shape(Hwav)[0]
tt = -np.arange(ntime)

# plot
plt.clf()
fig, ax = plt.subplots(figsize=(7.5,10))

plt.style.use('default')

cf = ax.contourf(xv, tt, Hwav, 10, cmap='hot_r')
ax.plot([xv[indx],xv[indx]], [tt[-1],tt[0]], color='black', ls='--', lw=1, 
        alpha=0.8)
ax.set_xlim([xv[0], xv[-1]])
ax.set_ylim([tt[-1], tt[0]])
ax.yaxis.set_ticks(np.arange(-100,20,20))
ax.yaxis.set_ticklabels(['100','80','60','40','20','0'])
ax.set_xlabel('Coordinate ($\mathregular{m}$)', fontsize=12, 
              fontname='Times New Roman', color='black')
ax.set_ylabel('Time (hours)', fontsize=12, 
              fontname='Times New Roman', color='black')
cbar = fig.colorbar(cf, ax=ax, orientation='horizontal', pad=0.08)
ylabel = 'Significant wave height ($\mathregular{m}$)'
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
fig.savefig('wave.png', dpi=300)
#fig.savefig('F1.pdf', dpi=600)
plt.show()