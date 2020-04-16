#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 31 14:35:12 2020

Draw the dynamics of sediment detachment and deposition.

@author: Zeli Tan
"""

import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset

z_1BF = -1.1    # 1BF is -1.1 m and 2BF is -2.1 m
z_2BF = -2.1
day0 = 9
day1 = 11

# read simulation outputs
filename = '/Users/tanz151/Python_maces/src/out_ecogeom_2002-12-01_2002-12-13.nc'
try:
    nc = Dataset(filename,'r')
    x = np.array(nc.variables['x'][0])
    zh = np.array(nc.variables['zh'][0,0,:])
finally:
    nc.close()
    
# find the site index
index1 = np.argmin(np.abs(zh - z_1BF))
index2 = np.argmin(np.abs(zh - z_2BF))
print(x[index1],x[index2])

filename = '/Users/tanz151/Python_maces/src/out_hydro_2002-12-01_2002-12-13.nc'
try:
    nc = Dataset(filename,'r')
    Esed_1BF = np.array(nc.variables['Esed'][0,day0*24:day1*24+1,index1])
    Dsed_1BF = np.array(nc.variables['Dsed'][0,day0*24:day1*24+1,index1])
    Esed_2BF = np.array(nc.variables['Esed'][0,day0*24:day1*24+1,index2])
    Dsed_2BF = np.array(nc.variables['Dsed'][0,day0*24:day1*24+1,index2])
finally:
    nc.close()
Esed_1BF = 3.6e6*np.reshape(Esed_1BF,(24*(day1-day0)+1))
Dsed_1BF = 3.6e6*np.reshape(Dsed_1BF,(24*(day1-day0)+1))
Esed_2BF = 3.6e6*np.reshape(Esed_2BF,(24*(day1-day0)+1))
Dsed_2BF = 3.6e6*np.reshape(Dsed_2BF,(24*(day1-day0)+1))

nt_model = np.size(Esed_1BF)
tt_model = np.arange(nt_model)

# plot water level, significant wave height, suspended sediment
plt.clf()
fig, axes = plt.subplots(2, 1, figsize=(6,6))

plt.style.use('default')

# 1BF
ax = axes[0]
ax.plot(tt_model, Esed_1BF, color='black', linestyle='-', linewidth=2, alpha=0.9)
ax.plot(tt_model, Dsed_1BF, color='C3', linestyle='-', marker='.', markersize=5)
ax.set_xlim(0, nt_model-1)
#ax.set_ylim(-100, 100)
ax.xaxis.set_ticks(np.arange(0,nt_model,24))
#ax.yaxis.set_ticks(np.linspace(-100,100,5))
ax.set_xticklabels(['12/10','12/11','12/12'])
#ax.xaxis.set_minor_locator(AutoMinorLocator(24))
ylabel = 'Sediment'
ax.set_ylabel(ylabel, fontsize=12, fontname='Times New Roman', color='black')
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(12) for label in labels]
[label.set_color('black') for label in labels]

# 2BF
ax = axes[1]
ax.plot(tt_model, Esed_2BF, color='black', linestyle='-', linewidth=2, alpha=0.9)
ax.plot(tt_model, Dsed_2BF, color='C3', linestyle='-', marker='.', markersize=5)
ax.set_xlim(0, nt_model-1)
#ax.set_ylim(-100, 100)
ax.xaxis.set_ticks(np.arange(0,nt_model,24))
#ax.yaxis.set_ticks(np.linspace(-100,100,5))
ax.set_xticklabels(['12/10','12/11','12/12'])
#ax.xaxis.set_minor_locator(AutoMinorLocator(24))
ylabel = 'Sediment'
ax.set_ylabel(ylabel, fontsize=12, fontname='Times New Roman', color='black')
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(12) for label in labels]
[label.set_color('black') for label in labels]

plt.tight_layout()
fig.savefig('venice_sed_sim.png', dpi=300)
#fig.savefig('venice_sed_sim.pdf', dpi=600)
plt.show()