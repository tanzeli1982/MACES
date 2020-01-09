#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  7 16:21:16 2020

Draw the dynamics of surface water elevation on the TAI platform

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
    
# read the simulated water depth
filename = '/Users/tanz151/Python_maces/src/out_hydro_2001-01-01_2001-01-06.nc'
try:
    nc = Dataset(filename,'r')
    xv = np.array(nc.variables['x'][:][0])
    Dp = np.array(nc.variables['h'][:][0])
finally:
    nc.close()
    
xv = xv[100:308]
Dp = Dp[:,100:308]
zh = zh[100:308]

zv = np.linspace(-3,2,1001)
nz = len(zv)
nx = len(xv)

indx = np.argmin(np.abs(zh))
xv = xv - xv[indx]

ntime = np.shape(Dp)[0]

mask = -1 * np.ones((nz,nx))
for ii in range(nx):
    indice = np.logical_and(zv>zh[ii], zv<=zh[ii]+Dp[0,ii])
    mask[indice,ii] = 1
Zpos = np.ma.masked_less(mask, 0)

dx = 0.5 * (xv[1] - xv[0])
dz = 0.5 * (zv[1] - zv[0])
extent = [xv[0]-dx, xv[-1]+dx, zv[0]-dz, zv[-1]+dz]

# plot
plt.clf()
fig, ax = plt.subplots(figsize=(7.5,5))

plt.style.use('default')

ax.plot(xv, zh, color='black', ls='-', lw=2, alpha=1.0)
ax.imshow(Zpos, cmap='Blues', vmin=-1, vmax=1, interpolation='none', 
          origin='lower', extent=extent, aspect='auto')
ax.set_xlim([xv[0], xv[-1]])
ax.set_ylim([-3, 2])
ax.xaxis.set_ticks(np.arange(-7500,12500,2500))
ax.yaxis.set_ticks(np.arange(-3,3,1))
ax.xaxis.set_minor_locator(AutoMinorLocator(5))
ax.yaxis.set_minor_locator(AutoMinorLocator(5))
ax.set_xlabel('Coordinate ($\mathregular{m}$)', fontsize=12, 
              fontname='Times New Roman', color='black')
ax.set_ylabel('Elevation (m.a.s.l.)', fontsize=12, 
              fontname='Times New Roman', color='black')
ax.tick_params(which='major', direction='in', length=6)
ax.tick_params(which='minor', direction='in', length=2)


labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(12) for label in labels]

plt.tight_layout()
fig.savefig('depth.png', dpi=300)
#fig.savefig('F1.pdf', dpi=600)
plt.show()
