#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 31 14:35:12 2020

Draw the dynamics of sediment detachment and deposition.

@author: Zeli Tan
"""

import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from matplotlib.ticker import AutoMinorLocator
from datetime import date

z_channel = -1.446      # channel elevation
z_marsh = 1.686         # marsh elevation
day0 = (date(2017,7,19) - date(2017,7,17)).days
day1 = (date(2017,7,23) - date(2017,7,17)).days

models = ['F06MOD', 'T03MOD', 'KM12MOD', 'F07MOD', 'VDK05MOD', 'DA07MOD', 'M12MOD']

# read simulation outputs
filename = '/Users/tanz151/Python_maces/src/maces_ecogeom_2017-07-17_2017-08-01_4097.nc'
try:
    nc = Dataset(filename,'r')
    x = np.array(nc.variables['x'][:])
    zh = np.array(nc.variables['zh'][0,:])
finally:
    nc.close()
    
# find the site index
index1 = np.argmin(np.abs(zh - z_channel))
index2 = np.argmin(np.abs(zh - z_marsh))

filename = '/Users/tanz151/Python_maces/src/maces_hydro_2017-07-17_2017-08-01_4097.nc'
try:
    nc = Dataset(filename,'r')
    sed_sim_c = np.array(nc.variables['TSM'][day0*24-1:day1*24-1,index1])
    sed_sim_m = np.array(nc.variables['TSM'][day0*24-1:day1*24-1,index2])
    tau_sim_c = np.array(nc.variables['tau'][day0*24-1:day1*24-1,index1])
finally:
    nc.close()
sed_sim_c = 1e3 * np.reshape(sed_sim_c,(24*(day1-day0)))    # mg/L
sed_sim_m = 1e3 * np.reshape(sed_sim_m,(24*(day1-day0)))    # mg/L

nt_model = np.size(sed_sim_c)
tt_model = np.arange(nt_model)

warnings.filterwarnings('ignore')

# read data
filename = r'/Users/tanz151/Documents/Projects/TAI_BGC/Data/Hydrodynamics_obs/' + \
    'PlumIsland/LawsPoint/EST-RO-TC-WE-TurbiditySuspSed.xlsx'
df = pd.read_excel(filename, sheet_name='EST-RO-TC-WE-TurbiditySuspSed', 
                   header=0, usecols='A,B,C,E,F')
df.columns = ['Date','Time','Site','SSC','Depth']
sed_c = np.array(df['SSC'])[20350:20734]
sed_m = np.array(df['SSC'])[87632:88016]
sed_obs_c = np.NaN * np.ones(24*(day1-day0))
sed_obs_m = np.NaN * np.ones(24*(day1-day0))
nt_obs = np.size(sed_obs_c)
for ii in range(nt_obs):
    sed_obs_c[ii] = np.nanmean(sed_c[4*ii:4*(ii+1)])
    sed_obs_m[ii] = np.nanmean(sed_m[4*ii:4*(ii+1)])

nt_obs = np.size(sed_obs_c)
tt_obs = np.arange(nt_obs)

# plot water level, significant wave height, suspended sediment
plt.clf()
fig, axes = plt.subplots(2, 1, figsize=(6,8))

plt.style.use('default')

# channel
ax = axes[0]
ax.plot(tt_model, sed_sim_c, color='black', linestyle='-', linewidth=2, alpha=0.9)
ax.plot(tt_obs, sed_obs_c, color='C3', linestyle='-', marker='.', markersize=5)
ax.set_xlim(0, nt_model)
ax.set_ylim(0, 50)
ax.xaxis.set_ticks(np.arange(0,nt_model+1,24))
ax.yaxis.set_ticks(np.linspace(0,50,6))
ax.set_xticklabels(['7/19','7/20','7/21','7/22','7/23'])
ax.xaxis.set_minor_locator(AutoMinorLocator(4))
ylabel = 'Suspended sediment ($\mathregular{mg}$ $\mathregular{l^{-1}}$)'
ax.set_ylabel(ylabel, fontsize=12, fontname='Times New Roman', color='black')
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(12) for label in labels]
[label.set_color('black') for label in labels]

# marsh
ax = axes[1]
ax.plot(tt_model, sed_sim_m, color='black', linestyle='-', linewidth=2, alpha=0.9)
ax.plot(tt_obs, sed_obs_m, color='C3', linestyle='-', marker='.', markersize=5)
ax.set_xlim(0, nt_model)
ax.set_ylim(0, 10)
ax.xaxis.set_ticks(np.arange(0,nt_model+1,24))
ax.yaxis.set_ticks(np.linspace(0,10,6))
ax.set_xticklabels(['7/19','7/20','7/21','7/22','7/23'])
ax.xaxis.set_minor_locator(AutoMinorLocator(4))
ax.set_xlabel('Time', fontsize=12, fontname='Times New Roman', color='black')
ylabel = 'Suspended sediment ($\mathregular{mg}$ $\mathregular{l^{-1}}$)'
ax.set_ylabel(ylabel, fontsize=12, fontname='Times New Roman', color='black')
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(12) for label in labels]
[label.set_color('black') for label in labels]

plt.tight_layout()
fig.savefig('F9.png', dpi=300)
#fig.savefig('F9.pdf', dpi=600)
plt.show()