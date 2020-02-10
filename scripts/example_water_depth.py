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
from scipy import interpolate
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

indx = np.argmin(np.abs(zh))
xv = xv - xv[indx]
ntime = np.shape(Dp)[0]

zwater = np.linspace(-3,2,1001)
xwater = np.linspace(xv[0],xv[-1],1001)
nz = len(zwater)
nx = len(xwater)

f = interpolate.interp1d(xv, zh)
zh_water = f(xwater)

# plot
plt.clf()
fig, axes = plt.subplots(4, 3, figsize=(8,10))

plt.style.use('default')

for mm in range(4):
    for nn in range(3):
        ax = axes[mm][nn]
        
        f = interpolate.interp1d(xv, Dp[mm*3+nn])
        Dp_water = f(xwater)

        mask = -1 * np.ones((nz,nx))
        for ii in range(nx):
            indice = np.logical_and(zwater>zh_water[ii], zwater<=zh_water[ii]+Dp_water[ii])
            mask[indice,ii] = 1
        Zpos = np.ma.masked_less(mask, 0)

        dx = 0.5 * (xwater[1] - xwater[0])
        dz = 0.5 * (zwater[1] - zwater[0])
        extent = [xwater[0]-dx, xwater[-1]+dx, zwater[0]-dz, zwater[-1]+dz]

        ax.plot(xv, zh, color='black', ls='-', lw=2, alpha=1.0)
        ax.imshow(Zpos, cmap='Blues', vmin=-1, vmax=1, interpolation='none', 
                  origin='lower', extent=extent, aspect='auto', alpha=0.8)
        ax.set_xlim([xv[0], xv[-1]])
        ax.set_ylim([-3, 2])
        ax.xaxis.set_ticks(np.arange(-5000,15000,5000))
        ax.yaxis.set_ticks(np.arange(-3,3,1))
        ax.xaxis.set_minor_locator(AutoMinorLocator(5))
        ax.yaxis.set_minor_locator(AutoMinorLocator(5))
        #ax.set_xlabel('Coordinate ($\mathregular{m}$)', fontsize=10, 
        #              fontname='Times New Roman', color='black')
        ax.set_ylabel('Elevation (m.a.s.l.)', fontsize=10, 
                      fontname='Times New Roman', color='black')
        ax.set_title('Time step '+ '{:02d}'.format(mm*3+nn+1), color='k', 
                     fontsize=10, fontname='Times New Roman')
        ax.tick_params(which='major', direction='in', length=6)
        ax.tick_params(which='minor', direction='in', length=2)

        labels = ax.get_xticklabels() + ax.get_yticklabels()
        [label.set_fontname('Times New Roman') for label in labels]
        [label.set_fontsize(10) for label in labels]

plt.tight_layout()
fig.savefig('depth.png', dpi=300)
#fig.savefig('F1.pdf', dpi=600)
plt.show()
