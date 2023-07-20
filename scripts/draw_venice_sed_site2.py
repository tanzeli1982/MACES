#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  1 18:20:59 2020

@author: Zeli Tan
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from netCDF4 import Dataset
from matplotlib.ticker import AutoMinorLocator
from datetime import date

z_1BF = -1.1    # 1BF is -1.1 m and 2BF is -2.1 m
z_2BF = -2.1
day0 = (date(2002,12,10) - date(2002,12,1)).days
day1 = (date(2002,12,12) - date(2002,12,1)).days

models = ['F06', 'T03', 'KM12', 'M12', 'F07', 'VDK05', 'DA07']
min_accr_sim = {}
sed_sim = {}

# read simulation outputs
rdir = '/Users/tanz151/Documents/Projects/TAI_BGC/Drafts/Outputs/VeniceLagoon/'
filename = rdir + 'maces_ecogeom_2002-12-01_2002-12-13_F06_466.nc'
try:
    nc = Dataset(filename,'r')
    x = np.array(nc.variables['x'][:])
    zh = np.array(nc.variables['zh'][0,:])
finally:
    nc.close()
    
index0 = np.argmin(np.abs(zh))
x = x - x[index0]
# find the site index
index1 = np.argmin(np.abs(zh - z_1BF))
index2 = np.argmin(np.abs(zh - z_2BF))
print(index1, index2)

for model in models:
    # hydrodynamics
    filename = rdir + 'maces_hydro_2002-12-01_2002-12-13_' + model + '_466.nc'
    try:
        nc = Dataset(filename,'r')
        sed_1BF = np.array(nc.variables['TSM'][day0*24:day1*24+1,index1])
        sed_2BF = np.array(nc.variables['TSM'][day0*24:day1*24+1,index2])
    finally:
        nc.close()
    sed_1BF = 1e3 * np.reshape(sed_1BF,(24*(day1-day0)+1))    # mg/L
    sed_2BF = 1e3 * np.reshape(sed_2BF,(24*(day1-day0)+1))    # mg/L
    sed_sim[model] = sed_1BF
    # eco-geomorphology
    filename = rdir + 'maces_ecogeom_2002-12-01_2002-12-13_' + model + '_466.nc'
    try:
        nc = Dataset(filename,'r')
        min_accr = 8.64e7*np.mean(np.array(nc.variables['Dsed'][day0:day1,:]),axis=0) - \
            8.64e7*np.mean(np.array(nc.variables['Esed'][day0:day1,:]),axis=0)   # g/m2/day
    finally:
        nc.close()
    min_accr_sim[model] = min_accr

nt_model = np.size(sed_sim['M12'])
tt_model = np.arange(nt_model)

# read data
filename = rdir + '1BF_OBS.xls'
df = pd.read_excel(filename, sheet_name='1BF', header=None, skiprows=range(3), 
                   usecols='A,B,F,O,Q')
df.columns = ['Time','Hmo','Hmax','hw','Turbidity']
sed_obs_1BF = np.array(df['Turbidity'])[5334:5526]  # mg/l

filename = rdir + '2BF_OBS.xls'
df = pd.read_excel(filename, sheet_name='2BF', header=None, skiprows=range(3), 
                   usecols='A,B,O,Q')
df.columns = ['Time','Hmo','hw','Turbidity']
sed_obs_2BF = np.array(df['Turbidity'])[5319:5511]

nt_obs = np.size(sed_obs_1BF)
tt_obs = np.arange(nt_obs)/4

# plotting
plt.clf()
fig = plt.figure(figsize=(8,5))

plt.style.use('default')
gs = gridspec.GridSpec(nrows=1, ncols=1)

#colors = ["#41bbc5", "#672d7e", "#b4d170", "#463df6", "#36f459", "#bf209e", 
#          "#256b33"]
colors = ['#7b85d4', '#f37738', '#83c995', '#d7369e', '#c4c9d8', '#859795',
          '#e9d043', '#ad5b50', '#e377c2']
linestyles = ['-', '--', '-.', ':', '-', '--', '-.']

# comparison of observed and simulated suspended sediment
print(np.nanmin(sed_obs_1BF), np.nanmax(sed_obs_1BF))

ax = fig.add_subplot(gs[0,0])
ax.plot(tt_obs, sed_obs_1BF, color='black', linestyle='-', linewidth=2, 
        marker='.', markersize=10)
handles = []
for key in sed_sim:
    indx = len(handles)
    h, = ax.plot(tt_model, sed_sim[key], color=colors[indx], 
                 linestyle=linestyles[indx], linewidth=2, alpha=1)
    handles.append(h)
    rmse = 0.0
    count = 0
    for ii, tt in enumerate(tt_obs):
        if tt in tt_model and np.isfinite(sed_obs_1BF[ii]):
            index = np.where(tt_model==tt)
            rmse = rmse + (sed_sim[key][index] - sed_obs_1BF[ii])**2
            count = count + 1
    rmse = np.sqrt( rmse / count )
    nrmse = rmse / np.nanmean(sed_obs_1BF)
    print(key, ': ', rmse, nrmse)
legend = ax.legend(handles, list(sed_sim.keys()), numpoints=1, loc=1, 
                   prop={'family':'Times New Roman', 'size':'large', 
                         'weight': 'bold'}, framealpha=0.0)
ax.set_xlim(0, nt_model)
ax.set_ylim(0, 120)
ax.xaxis.set_ticks(np.arange(0,nt_model+1,24))
ax.yaxis.set_ticks(np.linspace(0,120,5))
ax.set_xticklabels(['12/10','12/11','12/12'])
ax.xaxis.set_minor_locator(AutoMinorLocator(4))
ax.set_xlabel('Time', fontsize=16, fontname='Times New Roman', fontweight='bold')
ylabel = 'Suspended sediment ($\mathregular{mg}$ $\mathregular{l^{-1}}$)'
ax.set_ylabel(ylabel, fontsize=16, fontname='Times New Roman', fontweight='bold')
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(14) for label in labels]
[label.set_color('black') for label in labels]
[label.set_fontweight('bold') for label in labels]
#ax.text(0.05, 0.93, 'a', transform=ax.transAxes, fontsize=20,
#        fontname='Times New Roman', fontweight='bold')
ax.tick_params(which='major', direction='in', colors='xkcd:black', length=6, pad=8)
ax.tick_params(which='minor', direction='in', colors='xkcd:black')

## comparison of simulated sedimentation
#ax = axes[1]
#handles = []
#for key in min_accr_sim:
#    indx = len(handles)
#    h, = ax.plot(1e-3*x, min_accr_sim[key], color=colors[indx], 
#                 linestyle=linestyles[indx], linewidth=2, alpha=1)
#    handles.append(h)
#ax.set_xlim(-1, 2)
##ax.set_ylim(0, 150)
#ax.xaxis.set_ticks(np.linspace(-1,2,7))
##ax.yaxis.set_ticks(np.linspace(0,150,6))
#ax.xaxis.set_minor_locator(AutoMinorLocator(5))
#ax.set_xlabel('Distance ($\mathregular{km}$)', fontsize=12, 
#              fontname='Times New Roman', color='black')
#ylabel = 'Net sedimentation ($\mathregular{g}$ $\mathregular{m^{-2}}$ $\mathregular{day^{-1}}$)'
#ax.set_ylabel(ylabel, fontsize=12, fontname='Times New Roman', color='black')
#labels = ax.get_xticklabels() + ax.get_yticklabels()
#[label.set_fontname('Times New Roman') for label in labels]
#[label.set_fontsize(12) for label in labels]
#[label.set_color('black') for label in labels]
#ax.text(0.05, 0.93, 'b', transform=ax.transAxes, fontsize=20,
#        fontname='Times New Roman', fontweight='bold')
#ax.tick_params(which='major', direction='in', colors='xkcd:black', length=6, pad=8)
#ax.tick_params(which='minor', direction='in', colors='xkcd:black')
    
plt.tight_layout()
fig.savefig('F7.png', dpi=300)
fig.savefig('F7.jpg', dpi=600)
plt.show()
