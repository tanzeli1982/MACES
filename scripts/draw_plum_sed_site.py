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
import matplotlib.gridspec as gridspec
from netCDF4 import Dataset
from matplotlib.ticker import AutoMinorLocator
from datetime import date

z_channel = -1.446      # channel elevation
z_marsh = 1.686         # marsh elevation
day0 = (date(2017,7,19) - date(2017,7,17)).days
day1 = (date(2017,7,23) - date(2017,7,17)).days

models = ['F06', 'T03', 'KM12', 'M12', 'F07', 'VDK05', 'DA07']

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

min_accr_sim = {}
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
print(index1, index2)

sed_sim_c = {}
sed_sim_m = {}
U_sim_c = {}
U_sim_m = {}
Uwav_sim_c = {}
Uwav_sim_m = {}
tau_sim_c = {}
tau_sim_m = {}
for model in models:
    filename = rdir + 'maces_hydro_2017-07-17_2017-08-01_' + model + '_4097.nc'
    try:
        nc = Dataset(filename,'r')
        sed_c = np.array(nc.variables['TSM'][day0*24:day1*24,index1])
        sed_m = np.array(nc.variables['TSM'][day0*24:day1*24,index2])
        U_c = np.array(nc.variables['U'][day0*24:day1*24,index1])
        U_m = np.array(nc.variables['U'][day0*24:day1*24,index2])
        Uwav_c = np.array(nc.variables['Uwav'][day0*24:day1*24,index1])
        Uwav_m = np.array(nc.variables['Uwav'][day0*24:day1*24,index2])
        tau_c = np.array(nc.variables['tau'][day0*24:day1*24,index1])
        tau_m = np.array(nc.variables['tau'][day0*24:day1*24,index2])
    finally:
        nc.close()
    sed_sim_c[model] = 1e3 * np.reshape(sed_c,(24*(day1-day0)))    # mg/L
    sed_sim_m[model] = 1e3 * np.reshape(sed_m,(24*(day1-day0)))    # mg/L
    U_sim_c[model] = 100*np.reshape(U_c,(24*(day1-day0)))
    U_sim_m[model] = 100*np.reshape(U_m,(24*(day1-day0)))
    Uwav_sim_c[model] = 100*np.reshape(Uwav_c,(24*(day1-day0)))
    Uwav_sim_m[model] = 100*np.reshape(Uwav_m,(24*(day1-day0)))
    tau_sim_c[model] = 1e2*np.reshape(tau_c,(24*(day1-day0)))
    tau_sim_m[model] = 1e2*np.reshape(tau_m,(24*(day1-day0)))

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

gs = gridspec.GridSpec(nrows=1, ncols=2)

plt.style.use('default')

#colors = ["#aee39a", "#643176", "#4be32e", "#e72fc2", "#518413", "#7540fc", 
#          "#b3e61c"]
colors = ['#7b85d4', '#f37738', '#83c995', '#d7369e', '#c4c9d8', '#859795',
          '#e9d043', '#ad5b50', '#e377c2']
linestyles = ['-', '--', '-.', ':', '-', '--', '-.']

# channel
ax = fig.add_subplot(gs[0,0])
ax.plot(tt_obs, sed_obs_c, color='black', linestyle='-', linewidth=2, marker='.', 
        markersize=10)
handles = []
for key in sed_sim_c:
    #if key in ['F06', 'T03', 'KM12']:
    #    sed_sim_c[key][:] = 2.68165
    indx = len(handles)
    h, = ax.plot(tt_model, sed_sim_c[key], color=colors[indx], 
                 linestyle=linestyles[indx], linewidth=2, alpha=1)
    handles.append(h)
    rmse = 0.0
    count = 0
    for ii, tt in enumerate(tt_obs):
        if tt in tt_model and np.isfinite(sed_obs_c[ii]):
            index = np.where(tt_model==tt)
            rmse = rmse + (sed_sim_c[key][index] - sed_obs_c[ii])**2
            count = count + 1
    rmse = np.sqrt( rmse / count )
    nrmse = rmse / np.nanmean(sed_obs_c)
    print(key, ': ', rmse, nrmse)
#ax.plot(tt_model, U_sim_c['M12'], color='gray', linestyle='-', linewidth=1, alpha=0.9)
#ax.plot(tt_model, Uwav_sim_c['M12'], color='gray', linestyle=':', linewidth=1, alpha=0.9)
#ax.plot(tt_model, tau_sim_c['DA07'], color='gray', linestyle='-', linewidth=1, alpha=0.9)
ax.set_xlim(0, nt_model)
ax.set_ylim(0, 40)
ax.xaxis.set_ticks(np.arange(0,nt_model+1,24))
ax.yaxis.set_ticks(np.linspace(0,40,5))
ax.set_xticklabels(['7/19','7/20','7/21','7/22','7/23'])
ax.xaxis.set_minor_locator(AutoMinorLocator(4))
ax.set_xlabel('Time', fontsize=12, fontname='Times New Roman', color='black', 
              fontweight='bold')
ylabel = 'Suspended sediment ($\mathregular{mg}$ $\mathregular{l^{-1}}$)'
ax.set_ylabel(ylabel, fontsize=12, fontname='Times New Roman', color='black', 
              fontweight='bold')
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(12) for label in labels]
[label.set_color('black') for label in labels]
[label.set_fontweight('bold') for label in labels]
ax.text(0.05, 0.90, 'a', transform=ax.transAxes, fontsize=16,
        fontname='Times New Roman', fontweight='bold')
ax.tick_params(which='major', direction='in', colors='xkcd:black', length=6, pad=8)
ax.tick_params(which='minor', direction='in', colors='xkcd:black')

# marsh
ax = fig.add_subplot(gs[0,1])
ax.plot(tt_obs, sed_obs_m, color='black', linestyle='-', linewidth=2, marker='.', 
        markersize=10)
handles = []
for key in sed_sim_m:
    indx = len(handles)
    h, = ax.plot(tt_model, sed_sim_m[key], color=colors[indx], 
                 linestyle=linestyles[indx], linewidth=2, alpha=1)
    handles.append(h)
    rmse = 0.0
    count = 0
    for ii, tt in enumerate(tt_obs):
        if tt in tt_model and np.isfinite(sed_obs_m[ii]):
            index = np.where(tt_model==tt)
            rmse = rmse + (sed_sim_m[key][index] - sed_obs_m[ii])**2
            count = count + 1
    rmse = np.sqrt( rmse / count )
    nrmse = rmse / np.nanmean(sed_obs_m)
    print(key, ': ', rmse, nrmse)
#ax.plot(tt_model, U_sim_m['M12'], color='gray', linestyle='-', linewidth=1, alpha=0.9)
#ax.plot(tt_model, Uwav_sim_m['M12'], color='gray', linestyle=':', linewidth=1, alpha=0.9)
#ax.plot(tt_model, tau_sim_m['DA07'], color='gray', linestyle='-', linewidth=1, alpha=0.9)
legend = ax.legend(handles, list(sed_sim_c.keys()), numpoints=1, 
                   loc='upper right', ncol=2,
                   prop={'family':'Times New Roman', 'size':'large', 'weight': 'bold'}, 
                   framealpha=0.0)
ax.set_xlim(0, nt_model)
ax.set_ylim(0, 40)
ax.xaxis.set_ticks(np.arange(0,nt_model+1,24))
ax.yaxis.set_ticks(np.linspace(0,40,5))
ax.set_xticklabels(['7/19','7/20','7/21','7/22','7/23'])
ax.xaxis.set_minor_locator(AutoMinorLocator(4))
ax.set_xlabel('Time', fontsize=12, fontname='Times New Roman', color='black', 
              fontweight='bold')
ylabel = 'Suspended sediment ($\mathregular{mg}$ $\mathregular{l^{-1}}$)'
#ax.set_ylabel(ylabel, fontsize=12, fontname='Times New Roman', color='black')
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(12) for label in labels]
[label.set_color('black') for label in labels]
[label.set_fontweight('bold') for label in labels]
ax.text(0.05, 0.90, 'b', transform=ax.transAxes, fontsize=16,
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
fig.savefig('F8.jpg', dpi=600)
plt.show()