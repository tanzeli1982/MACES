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

models = ['F06', 'T03', 'KM12']#, 'M12', 'F07', 'VDK05', 'DA07']
min_accr_sim = {}
sed_sim_c = {}
sed_sim_m = {}

# read simulation outputs
rdir = '/Users/tanz151/Documents/Projects/TAI_BGC/Drafts/Outputs/PlumIsland/'
filename = rdir + 'maces_ecogeom_2017-07-17_2017-08-01_F06_4097.nc'
try:
    nc = Dataset(filename,'r')
    x = 1e-3 * np.array(nc.variables['x'][:])
    zh = np.array(nc.variables['zh'][0,:])
finally:
    nc.close()
index0 = np.argmin(np.abs(zh))
x = x - x[index0]

for model in models:
    filename = rdir + 'maces_ecogeom_2017-07-17_2017-08-01_' + model + '_4097.nc'
    try:
        nc = Dataset(filename,'r')
        min_accr = 8.64e7*np.mean(np.array(nc.variables['Dsed'][day0:day1,:]),axis=0) - \
            8.64e7*np.mean(np.array(nc.variables['Esed'][day0:day1,:]),axis=0)   # g/m2/day
    finally:
        nc.close()
    min_accr_sim[model] = min_accr
    
# find the site index
index1 = np.argmin(np.abs(zh - z_channel))
index2 = np.argmin(np.abs(zh - z_marsh))

for model in models:
    filename = rdir + 'maces_hydro_2017-07-17_2017-08-01_' + model + '_4097.nc'
    try:
        nc = Dataset(filename,'r')
        sed_c = np.array(nc.variables['TSM'][day0*24:day1*24,index1])
        sed_m = np.array(nc.variables['TSM'][day0*24:day1*24,index2])
    finally:
        nc.close()
    sed_sim_c[model] = 1e3 * np.reshape(sed_c,(24*(day1-day0)))    # mg/L
    sed_sim_m[model] = 1e3 * np.reshape(sed_m,(24*(day1-day0)))    # mg/L

nt_model = np.size(sed_sim_c['F06'])
tt_model = np.arange(nt_model)

warnings.filterwarnings('ignore')

# read data
filename = rdir + 'EST-RO-TC-WE-TurbiditySuspSed.xlsx'
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
fig = plt.figure(figsize=(8,3.5))

plt.style.use('default')

#colors = ["#aee39a", "#643176", "#4be32e", "#e72fc2", "#518413", "#7540fc", 
#          "#b3e61c"]
colors = ['#7b85d4', '#f37738', '#83c995', '#d7369e', '#c4c9d8', '#859795',
          '#e9d043', '#ad5b50', '#e377c2']
linestyles = ['-', '--', '-.', ':', '-', '--', '-.']

# channel
ax = plt.subplot2grid((1,2),(0,0))
ax.plot(tt_obs, sed_obs_c, color='black', linestyle='-', linewidth=2, marker='.', 
        markersize=8)
handles = []
for key in sed_sim_c:
    #if key in ['F06', 'T03', 'KM12']:
    #    sed_sim_c[key][:] = 2.68165
    indx = len(handles)
    h, = ax.plot(tt_model, sed_sim_c[key], color=colors[indx], 
                 linestyle=linestyles[indx], linewidth=2, alpha=1)
    handles.append(h)
ax.set_xlim(0, nt_model)
ax.set_ylim(0, 50)
ax.xaxis.set_ticks(np.arange(0,nt_model+1,24))
ax.yaxis.set_ticks(np.linspace(0,50,6))
ax.set_xticklabels(['7/19','7/20','7/21','7/22','7/23'])
ax.xaxis.set_minor_locator(AutoMinorLocator(4))
ax.set_xlabel('Time', fontsize=12, fontname='Times New Roman', color='black')
ylabel = 'Suspended sediment ($\mathregular{mg}$ $\mathregular{l^{-1}}$)'
ax.set_ylabel(ylabel, fontsize=12, fontname='Times New Roman', color='black')
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(12) for label in labels]
[label.set_color('black') for label in labels]
ax.text(0.05, 0.93, 'a', transform=ax.transAxes, fontsize=16,
        fontname='Times New Roman', fontweight='bold')
ax.tick_params(which='major', direction='in', colors='xkcd:black', length=6, pad=8)
ax.tick_params(which='minor', direction='in', colors='xkcd:black')

# marsh
ax = plt.subplot2grid((1,2),(0,1))
ax.plot(tt_obs, sed_obs_m, color='black', linestyle='-', linewidth=2, marker='.', 
        markersize=8)
handles = []
for key in sed_sim_m:
    indx = len(handles)
    h, = ax.plot(tt_model, sed_sim_m[key], color=colors[indx], 
                 linestyle=linestyles[indx], linewidth=2, alpha=1)
    handles.append(h)
legend = ax.legend(handles, list(sed_sim_c.keys()), numpoints=1, loc=1, 
                   prop={'family':'Times New Roman', 'size':'large'}, 
                   framealpha=0.0)
ax.set_xlim(0, nt_model)
ax.set_ylim(0, 20)
ax.xaxis.set_ticks(np.arange(0,nt_model+1,24))
ax.yaxis.set_ticks(np.linspace(0,20,6))
ax.set_xticklabels(['7/19','7/20','7/21','7/22','7/23'])
ax.xaxis.set_minor_locator(AutoMinorLocator(4))
ax.set_xlabel('Time', fontsize=12, fontname='Times New Roman', color='black')
ylabel = 'Suspended sediment ($\mathregular{mg}$ $\mathregular{l^{-1}}$)'
#ax.set_ylabel(ylabel, fontsize=12, fontname='Times New Roman', color='black')
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(12) for label in labels]
[label.set_color('black') for label in labels]
ax.text(0.05, 0.93, 'b', transform=ax.transAxes, fontsize=16,
        fontname='Times New Roman', fontweight='bold')
ax.tick_params(which='major', direction='in', colors='xkcd:black', length=6, pad=8)
ax.tick_params(which='minor', direction='in', colors='xkcd:black')

## sedimentation
#ax = plt.subplot2grid((2,2), (1,0), colspan=2)
#handles = []
#for key in min_accr_sim:
#    indx = len(handles)
#    h, = ax.plot(x, min_accr_sim[key], color=colors[indx], 
#                 linestyle=linestyles[indx], linewidth=2, alpha=1)
#    handles.append(h)
#ax.set_xlim(-0.2, 0.2)
##ax.set_ylim(0, 150)
#ax.xaxis.set_ticks(np.linspace(-0.2,0.2,5))
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
#ax.text(0.03, 0.93, 'c', transform=ax.transAxes, fontsize=20,
#        fontname='Times New Roman', fontweight='bold')
#ax.tick_params(which='major', direction='in', colors='xkcd:black', length=6, pad=8)
#ax.tick_params(which='minor', direction='in', colors='xkcd:black')

plt.tight_layout()
fig.savefig('F8.png', dpi=300)
#fig.savefig('F8.jpg', dpi=600)
plt.show()