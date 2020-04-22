#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  1 18:20:59 2020

@author: Zeli Tan
"""

import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset

day0 = 9
day1 = 11
# read simulation outputs
filename = '/Users/tanz151/Python_maces/src/maces_ecogeom_2002-12-01_2002-12-13_466.nc'
try:
    nc = Dataset(filename,'r')
    x = np.array(nc.variables['x'][:])
    zh = np.array(nc.variables['zh'][0,:])
finally:
    nc.close()
    
filename = '/Users/tanz151/Python_maces/src/maces_hydro_2002-12-01_2002-12-13_466.nc'
try:
    nc = Dataset(filename,'r')
    Esed = 3.6e6*np.array(nc.variables['Esed'][day0*24:day1*24,:])
    Dsed = 3.6e6*np.array(nc.variables['Dsed'][day0*24:day1*24,:])
finally:
    nc.close()

t0 = 0

plt.clf()
fig, axes = plt.subplots(4, 3, figsize=(8,10))

plt.style.use('default')

for ii in range(4):
    for jj in range(3):
        ax = axes[ii][jj]
        ax.plot(1e-3*x, Esed[t0], color='black', linestyle='-', linewidth=2, alpha=0.9)
        ax.plot(1e-3*x, Dsed[t0], color='C3', linestyle='-', marker='.', markersize=5)
        ax.set_xlim(10, 20)
        ax.set_ylim(0,200)
        ax.xaxis.set_ticks(np.linspace(10,20,5))
        ax.yaxis.set_ticks(np.linspace(0,200,5))
        if jj>0:
            ax.set_yticklabels([])
        if ii<3:
            ax.set_xticklabels([])
        labels = ax.get_xticklabels() + ax.get_yticklabels()
        [label.set_fontname('Times New Roman') for label in labels]
        [label.set_fontsize(12) for label in labels]
        [label.set_color('black') for label in labels]
        t0 = t0 + 1
    
plt.tight_layout()
fig.savefig('venice_Esed_Dsed.png', dpi=300)
#fig.savefig('venice_Esed_Dsed.pdf', dpi=600)
plt.show()