#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 31 14:35:12 2020

Draw the dynamics of sediment detachment and deposition.

@author: Zeli Tan
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from matplotlib.ticker import AutoMinorLocator

z_1BF = -1.1    # 1BF is -1.1 m and 2BF is -2.1 m
z_2BF = -2.1
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
    
# find the site index
index1 = np.argmin(np.abs(zh - z_1BF))
index2 = np.argmin(np.abs(zh - z_2BF))

filename = '/Users/tanz151/Python_maces/src/maces_hydro_2002-12-01_2002-12-13_466.nc'
try:
    nc = Dataset(filename,'r')
    sed_1BF = np.array(nc.variables['TSM'][day0*24:day1*24+1,index1])
    sed_2BF = np.array(nc.variables['TSM'][day0*24:day1*24+1,index2])
finally:
    nc.close()
sed_1BF = 1e3 * np.reshape(sed_1BF,(24*(day1-day0)+1))    # mg/L
sed_2BF = 1e3 * np.reshape(sed_2BF,(24*(day1-day0)+1))    # mg/L

nt_model = np.size(sed_1BF)
tt_model = np.arange(nt_model)

# read data
filename = r'/Users/tanz151/Documents/Projects/TAI_BGC/Data/Hydrodynamics_obs/' + \
    'VeniceLagoon/1BF_OBS.xls'
df = pd.read_excel(filename, sheet_name='1BF', header=None, skiprows=range(3), 
                   usecols='A,B,F,O,Q')
df.columns = ['Time','Hmo','Hmax','hw','Turbidity']
sed_obs_1BF = np.array(df['Turbidity'])[5334:5526]  # mg/l

filename = r'/Users/tanz151/Documents/Projects/TAI_BGC/Data/Hydrodynamics_obs/' + \
    'VeniceLagoon/2BF_OBS.xls'
df = pd.read_excel(filename, sheet_name='2BF', header=None, skiprows=range(3), 
                   usecols='A,B,O,Q')
df.columns = ['Time','Hmo','hw','Turbidity']
sed_obs_2BF = np.array(df['Turbidity'])[5319:5511]

nt_obs = np.size(sed_obs_1BF)
tt_obs = np.arange(nt_obs)/4

# plot water level, significant wave height, suspended sediment
plt.clf()
fig, axes = plt.subplots(1, 2, figsize=(8,3))

plt.style.use('default')

# 1BF
ax = axes[0]
ax.plot(tt_model, sed_1BF, color='black', linestyle='-', linewidth=2, alpha=0.9)
ax.plot(tt_obs, sed_obs_1BF, color='C3', linestyle='-', marker='.', markersize=5)
ax.set_xlim(0, nt_model-1)
ax.set_ylim(0, 150)
ax.xaxis.set_ticks(np.arange(0,nt_model,24))
ax.yaxis.set_ticks(np.linspace(0,150,6))
ax.set_xticklabels(['12/10','12/11','12/12'])
ax.xaxis.set_minor_locator(AutoMinorLocator(4))
ylabel = 'Suspended sediment ($\mathregular{mg}$ $\mathregular{l^{-1}}$)'
ax.set_ylabel(ylabel, fontsize=12, fontname='Times New Roman', color='black')
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(12) for label in labels]
[label.set_color('black') for label in labels]

# 2BF
ax = axes[1]
ax.plot(tt_model, sed_2BF, color='black', linestyle='-', linewidth=2, alpha=0.9)
ax.plot(tt_obs, sed_obs_2BF, color='C3', linestyle='-', marker='.', markersize=5)
ax.set_xlim(0, nt_model-1)
ax.set_ylim(0, 150)
ax.xaxis.set_ticks(np.arange(0,nt_model,24))
ax.yaxis.set_ticks(np.linspace(0,150,6))
ax.set_xticklabels(['12/10','12/11','12/12'])
ax.xaxis.set_minor_locator(AutoMinorLocator(4))
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(12) for label in labels]
[label.set_color('black') for label in labels]

plt.tight_layout()
fig.savefig('F5.png', dpi=300)
#fig.savefig('F5.pdf', dpi=600)
plt.show()