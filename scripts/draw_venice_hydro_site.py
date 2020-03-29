#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 23 10:51:30 2020

Draw simulated hydrodynamics at the specified landscape locations

@author: Zeli Tan
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from netCDF4 import Dataset

z_1BF = -1.1    # 1BF is -1.1 m and 2BF is -2.1 m
z_2BF = -2.1
day0 = 343
day1 = 345

# read simulation outputs
filename = '/Users/tanz151/Python_maces/src/out_ecogeom_2002-01-01_2004-01-01.nc'
try:
    nc = Dataset(filename,'r')
    x = np.array(nc.variables['x'][0])
    zh = np.array(nc.variables['zh'][0,0,:])
finally:
    nc.close()
    
# find the site index
index1 = np.argmin(np.abs(zh - z_1BF))
index2 = np.argmin(np.abs(zh - z_2BF))

filename = '/Users/tanz151/Python_maces/src/out_hydro_2002-01-01_2004-01-01.nc'
try:
    nc = Dataset(filename,'r')
    h_1BF = np.array(nc.variables['h'][0,day0*24:day1*24+1,index1])
    Hwav_1BF = np.array(nc.variables['Hwav'][0,day0*24:day1*24+1,index1])
    sed_1BF = np.array(nc.variables['TSM'][0,day0*24:day1*24+1,index1])
    tau_1BF = np.array(nc.variables['tau'][0,day0*24:day1*24+1,index1])
    h_2BF = np.array(nc.variables['h'][0,day0*24:day1*24+1,index2])
    Hwav_2BF = np.array(nc.variables['Hwav'][0,day0*24:day1*24+1,index2])
    sed_2BF = np.array(nc.variables['TSM'][0,day0*24:day1*24+1,index2])
    tau_2BF = np.array(nc.variables['tau'][0,day0*24:day1*24+1,index2])
finally:
    nc.close()
h_1BF = 100 * np.reshape(h_1BF,(24*(day1-day0)+1))    # cm
Hwav_1BF = 100 * np.reshape(Hwav_1BF,(24*(day1-day0)+1))  # cm
sed_1BF = 1e3 * np.reshape(sed_1BF,(24*(day1-day0)+1))    # mg/l
tau_1BF = np.reshape(tau_1BF,(24*(day1-day0)+1))
h_2BF = 100 * np.reshape(h_2BF,(24*(day1-day0)+1))    # cm
Hwav_2BF = 100 * np.reshape(Hwav_2BF,(24*(day1-day0)+1))  # cm
sed_2BF = 1e3 * np.reshape(sed_2BF,(24*(day1-day0)+1))    # mg/l
tau_2BF = np.reshape(tau_2BF,(24*(day1-day0)+1))

nt_model = np.size(h_1BF)
tt_model = np.arange(nt_model)

# read data
filename = r'/Users/tanz151/Documents/Projects/TAI_BGC/Data/Hydrodynamics_obs/' + \
    'VeniceLagoon/1BF_OBS.xls'
df = pd.read_excel(filename, sheet_name='1BF', header=None, skiprows=range(3), 
                   usecols='A,B,F,O,Q')
df.columns = ['Time','Hmo','Hmax','hw','Turbidity']
Hwav_obs_1BF = 100 * np.array(df['Hmo'])[5334:5526] # cm
sed_obs_1BF = np.array(df['Turbidity'])[5334:5526]  # mg/l

filename = '/Users/tanz151/Documents/Projects/TAI_BGC/Data/Hydrodynamics_obs/' + \
    'VeniceLagoon/WaterLevelClose1BF.xls'
df = pd.read_excel(filename, sheet_name='Valori orari 2002', header=None, 
                   skiprows=range(4), usecols='A:C')
df.columns = ['Date','Hour','hw']
h_obs_1BF = 100 * np.array(df['hw'])[8231:8303]
nt_obs2 = np.size(h_obs_1BF)
tt_obs2 = np.arange(nt_obs2)

filename = r'/Users/tanz151/Documents/Projects/TAI_BGC/Data/Hydrodynamics_obs/' + \
    'VeniceLagoon/2BF_OBS.xls'
df = pd.read_excel(filename, sheet_name='2BF', header=None, skiprows=range(3), 
                   usecols='A,B,O,Q')
df.columns = ['Time','Hmo','hw','Turbidity']
Hwav_obs_2BF = 100 * np.array(df['Hmo'])[5319:5511]
h_obs_2BF = 100 * np.array(df['hw'])[5319:5511]
sed_obs_2BF = np.array(df['Turbidity'])[5319:5511]

h_obs_1BF = h_obs_1BF + 110
h_obs_2BF = h_obs_2BF + 210

nt_obs = np.size(Hwav_obs_1BF)
tt_obs = np.arange(nt_obs)/4

# plot water level, significant wave height, suspended sediment
plt.clf()
fig, axes = plt.subplots(4, 2, figsize=(8,10))

# 1BF
ax = axes[0][0]
ax.plot(tt_model, h_1BF, color='black', linestyle='-', linewidth=2, alpha=0.9)
ax.plot(tt_obs2, h_obs_1BF, color='C3', linestyle='-', marker='.', markersize=5)
ax.set_xlim(0, nt_model-1)
#ax.set_ylim(-100, 100)
ax.xaxis.set_ticks(np.arange(0,nt_model,24))
#ax.yaxis.set_ticks(np.linspace(-100,100,5))
ax.set_xticklabels(['12/10','12/11','12/12'])
#ax.xaxis.set_minor_locator(AutoMinorLocator(24))
ylabel = 'Water depth ($\mathregular{cm}$)'
ax.set_ylabel(ylabel, fontsize=12, fontname='Times New Roman', color='black')
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(12) for label in labels]
[label.set_color('black') for label in labels]

ax = axes[1][0]
ax.plot(tt_model, Hwav_1BF, color='black', linestyle='-', linewidth=2, alpha=0.9)
ax.plot(tt_obs, Hwav_obs_1BF, color='C3', linestyle='-', marker='.', markersize=5)
ax.set_xlim(0, nt_model-1)
#ax.set_ylim(0, 100)
ax.xaxis.set_ticks(np.arange(0,nt_model,24))
#ax.yaxis.set_ticks(np.linspace(0,100,6))
ax.set_xticklabels(['12/10','12/11','12/12'])
#ax.xaxis.set_minor_locator(AutoMinorLocator(24))
ylabel = 'Significant wave height ($\mathregular{cm}$)'
ax.set_ylabel(ylabel, fontsize=12, fontname='Times New Roman', color='black')
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(12) for label in labels]
[label.set_color('black') for label in labels]

ax = axes[2][0]
ax.plot(tt_model, sed_1BF, color='black', linestyle='-', linewidth=2, alpha=0.9)
ax.plot(tt_obs, sed_obs_1BF, color='C3', linestyle='-', marker='.', markersize=5)
ax.set_xlim(0, nt_model-1)
#ax.set_ylim(0, 150)
ax.xaxis.set_ticks(np.arange(0,nt_model,24))
#ax.yaxis.set_ticks(np.linspace(0,150,6))
ax.set_xticklabels(['12/10','12/11','12/12'])
#ax.xaxis.set_minor_locator(AutoMinorLocator(24))
ylabel = 'Suspended sediment ($\mathregular{mg}$ $\mathregular{{l}^{-1}}$)'
ax.set_ylabel(ylabel, fontsize=12, fontname='Times New Roman', color='black')
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(12) for label in labels]
[label.set_color('black') for label in labels]

ax = axes[3][0]
ax.plot(tt_model, tau_1BF, color='black', linestyle='-', linewidth=2, alpha=0.9)
ax.set_xlim(0, nt_model-1)
#ax.set_ylim(0, 1)
ax.xaxis.set_ticks(np.arange(0,nt_model,24))
#ax.yaxis.set_ticks(np.linspace(0,1,6))
ax.set_xticklabels(['12/10','12/11','12/12'])
#ax.xaxis.set_minor_locator(AutoMinorLocator(24))
ylabel = 'Bottom shear stress ($\mathregular{Pa}$)'
ax.set_ylabel(ylabel, fontsize=12, fontname='Times New Roman', color='black')
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(12) for label in labels]
[label.set_color('black') for label in labels]

# 2BF
ax = axes[0][1]
ax.plot(tt_model, h_2BF, color='black', linestyle='-', linewidth=2, alpha=0.9)
ax.plot(tt_obs, h_obs_2BF, color='C3', linestyle='-', marker='.', markersize=5)
ax.set_xlim(0, nt_model-1)
#ax.set_ylim(-100, 100)
ax.xaxis.set_ticks(np.arange(0,nt_model+1,24))
#ax.yaxis.set_ticks(np.linspace(-100,100,5))
ax.set_xticklabels(['12/10','12/11','12/12'])
#ax.xaxis.set_minor_locator(AutoMinorLocator(24))
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(12) for label in labels]
[label.set_color('black') for label in labels]

ax = axes[1][1]
ax.plot(tt_model, Hwav_2BF, color='black', linestyle='-', linewidth=2, alpha=0.9)
ax.plot(tt_obs, Hwav_obs_2BF, color='C3', linestyle='-', marker='.', markersize=5)
ax.set_xlim(0, nt_model-1)
#ax.set_ylim(0, 100)
ax.xaxis.set_ticks(np.arange(0,nt_model+1,24))
#ax.yaxis.set_ticks(np.linspace(0,100,6))
ax.set_xticklabels(['12/10','12/11','12/12'])
#ax.xaxis.set_minor_locator(AutoMinorLocator(24))
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(12) for label in labels]
[label.set_color('black') for label in labels]

ax = axes[2][1]
ax.plot(tt_model, sed_2BF, color='black', linestyle='-', linewidth=2, alpha=0.9)
ax.plot(tt_obs, sed_obs_2BF, color='C3', linestyle='-', marker='.', markersize=5)
ax.set_xlim(0, nt_model-1)
#ax.set_ylim(0, 150)
ax.xaxis.set_ticks(np.arange(0,nt_model+1,24))
#ax.yaxis.set_ticks(np.linspace(0,150,6))
ax.set_xticklabels(['12/10','12/11','12/12'])
#ax.xaxis.set_minor_locator(AutoMinorLocator(24))
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(12) for label in labels]
[label.set_color('black') for label in labels]

ax = axes[3][1]
ax.plot(tt_model, tau_2BF, color='black', linestyle='-', linewidth=2, alpha=0.9)
ax.set_xlim(0, nt_model-1)
#ax.set_ylim(0, 1)
ax.xaxis.set_ticks(np.arange(0,nt_model,24))
#ax.yaxis.set_ticks(np.linspace(0,1,6))
ax.set_xticklabels(['12/10','12/11','12/12'])
#ax.xaxis.set_minor_locator(AutoMinorLocator(24))
ylabel = 'Bottom shear stress ($\mathregular{Pa}$)'
ax.set_ylabel(ylabel, fontsize=12, fontname='Times New Roman', color='black')
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(12) for label in labels]
[label.set_color('black') for label in labels]

plt.tight_layout()
fig.savefig('venice_hydro_sim.png', dpi=300)
#fig.savefig('venice_hydro_sim.pdf', dpi=600)
plt.show()