#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 19 22:37:31 2020

Compare simulated and benchmark estimate of suspended sediment

@author: Zeli Tan
"""

import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from matplotlib.ticker import AutoMinorLocator
from datetime import date

z1 = -0.22  # channel
z2 = 0.05   # mangrove edge
z3 = 0.38   # mangrove interior
z4 = 0.65   # salt marsh edge

models = ['F06', 'T03', 'KM12', 'F07', 'VDK05', 'DA07', 'M12']

day0 = (date(2004,9,28) - date(2004,9,25)).days
day1 = (date(2004,10,1) - date(2004,9,25)).days
# read simulated water depth and suspended sediment
rdir = '/Users/tanz151/Documents/Projects/TAI_BGC/Data/Hydrodynamics_obs/HunterEstuary/Outputs/'
filename = rdir + 'maces_ecogeom_2004-09-25_2004-10-06_11099.F06.nc'
try:
    nc = Dataset(filename,'r')
    x = np.array(nc.variables['x'][:])
    zh = np.array(nc.variables['zh'][0,:])
    #Bag = 1e3 * np.mean(np.array(nc.variables['Bag'][day0:day1,:]),axis=0)
    #min_accr = 8.64e7*np.mean(np.array(nc.variables['Dsed'][day0:day1,:]),axis=0) - \
    #    8.64e7*np.mean(np.array(nc.variables['Esed'][day0:day1,:]),axis=0)   # g/m2/day
finally:
    nc.close()
    
# find the site index
index1 = np.argmin(np.abs(zh - z1))
index2 = np.argmin(np.abs(zh - z2))
index3 = np.argmin(np.abs(zh - z3))
index4 = np.argmin(np.abs(zh - z4))

sed1_sim = {}
sed2_sim = {}
sed3_sim = {}
sed4_sim = {}
for model in models:
    filename = rdir + 'maces_hydro_2004-09-25_2004-10-06_11099.' + model + '.nc'
    try:
        nc = Dataset(filename, 'r')
        sed1 = 1e3*np.array(nc.variables['TSM'][day0*24:day1*24,index1])
        sed2 = 1e3*np.array(nc.variables['TSM'][day0*24:day1*24,index2])
        sed3 = 1e3*np.array(nc.variables['TSM'][day0*24:day1*24,index3])
        sed4 = 1e3*np.array(nc.variables['TSM'][day0*24:day1*24,index4])
    finally:
        nc.close()
    sed1_sim[model] = sed1
    sed2_sim[model] = sed2
    sed3_sim[model] = sed3
    sed4_sim[model] = sed4
nt_model = np.size(sed1_sim['M12'])
tt_model = np.arange(nt_model)

# read benchmark estimates
filename = '/Users/tanz151/Documents/Projects/TAI_BGC/Data/Hydrodynamics_obs/' + \
    'HunterEstuary/SSC_benchmark_2004-09-28_2004-09-30.nc'
try:
    nc = Dataset(filename,'r')
    sed1_obs = np.array(nc.variables['SSC'][:,26,41])
    sed2_obs = np.array(nc.variables['SSC'][:,26,36])
    sed3_obs = np.array(nc.variables['SSC'][:,26,33])
    sed4_obs = np.array(nc.variables['SSC'][:,27,27])
finally:
    nc.close()
nt_obs = np.size(sed1_obs)
tt_obs = np.arange(nt_obs) / 4.

# plot
plt.clf()
fig, axes = plt.subplots(2, 2, figsize=(8,7))

plt.style.use('default')

colors = ['#7b85d4', '#f37738', '#83c995', '#d7369e', '#c4c9d8', '#859795',
          '#e9d043', '#ad5b50', '#e377c2']
linestyles = ['-', '--', '-.', ':', '-', '--', '-.']

# channel
ax = axes[0][0]
ax.plot(tt_obs, sed1_obs, color='black', linestyle='-', linewidth=2, marker='.', 
        markersize=8)
handles = []
for key in sed1_sim:
    indx = len(handles)
    h, = ax.plot(tt_model, sed1_sim[key], color=colors[indx], 
                 linestyle=linestyles[indx], linewidth=2, alpha=1)
    handles.append(h)
ax.set_xlim(0, nt_model)
#ax.set_ylim(0, 100)
ax.xaxis.set_ticks(np.arange(0,nt_model+1,24))
#ax.yaxis.set_ticks(np.linspace(0,100,6))
ax.set_xticklabels(['9/28','9/29','9/30','10/1'])
ax.xaxis.set_minor_locator(AutoMinorLocator(4))
#ax.set_xlabel('Time', fontsize=11, fontname='Times New Roman', color='black')
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

# mangrove edge
ax = axes[0][1]
ax.plot(tt_obs, sed2_obs, color='black', linestyle='-', linewidth=2, marker='.', 
        markersize=8)
handles = []
for key in sed2_sim:
    indx = len(handles)
    h, = ax.plot(tt_model, sed2_sim[key], color=colors[indx], 
                 linestyle=linestyles[indx], linewidth=2, alpha=1)
    handles.append(h)
ax.set_xlim(0, nt_model)
#ax.set_ylim(0, 100)
ax.xaxis.set_ticks(np.arange(0,nt_model+1,24))
#ax.yaxis.set_ticks(np.linspace(0,100,6))
ax.set_xticklabels(['9/28','9/29','9/30','10/1'])
ax.xaxis.set_minor_locator(AutoMinorLocator(4))
#ax.set_xlabel('Time', fontsize=11, fontname='Times New Roman', color='black')
#ylabel = 'Suspended sediment ($\mathregular{mg}$ $\mathregular{l^{-1}}$)'
#ax.set_ylabel(ylabel, fontsize=12, fontname='Times New Roman', color='black')
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(12) for label in labels]
[label.set_color('black') for label in labels]
ax.text(0.05, 0.93, 'b', transform=ax.transAxes, fontsize=16,
        fontname='Times New Roman', fontweight='bold')
ax.tick_params(which='major', direction='in', colors='xkcd:black', length=6, pad=8)
ax.tick_params(which='minor', direction='in', colors='xkcd:black')

# mangrove interior
ax = axes[1][0]
ax.plot(tt_obs, sed3_obs, color='black', linestyle='-', linewidth=2, marker='.', 
        markersize=8)
handles = []
for key in sed3_sim:
    indx = len(handles)
    h, = ax.plot(tt_model, sed3_sim[key], color=colors[indx], 
                 linestyle=linestyles[indx], linewidth=2, alpha=1)
    handles.append(h)
ax.set_xlim(0, nt_model)
#ax.set_ylim(0, 50)
ax.xaxis.set_ticks(np.arange(0,nt_model+1,24))
#ax.yaxis.set_ticks(np.linspace(0,50,6))
ax.set_xticklabels(['9/28','9/29','9/30','10/1'])
ax.xaxis.set_minor_locator(AutoMinorLocator(4))
ax.set_xlabel('Time', fontsize=12, fontname='Times New Roman', color='black')
ylabel = 'Suspended sediment ($\mathregular{mg}$ $\mathregular{l^{-1}}$)'
ax.set_ylabel(ylabel, fontsize=12, fontname='Times New Roman', color='black')
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(12) for label in labels]
[label.set_color('black') for label in labels]
ax.text(0.05, 0.93, 'c', transform=ax.transAxes, fontsize=16,
        fontname='Times New Roman', fontweight='bold')
ax.tick_params(which='major', direction='in', colors='xkcd:black', length=6, pad=8)
ax.tick_params(which='minor', direction='in', colors='xkcd:black')

# marsh edge
ax = axes[1][1]
ax.plot(tt_obs, sed4_obs, color='black', linestyle='-', linewidth=2, marker='.', 
        markersize=8)
handles = []
for key in sed4_sim:
    indx = len(handles)
    h, = ax.plot(tt_model, sed4_sim[key], color=colors[indx], 
                 linestyle=linestyles[indx], linewidth=2, alpha=1)
    handles.append(h)
legend = ax.legend(handles, list(sed4_sim.keys()), numpoints=1, loc=1, 
                   prop={'family':'Times New Roman', 'size':'large'}, 
                   framealpha=0.0)
ax.set_xlim(0, nt_model)
#ax.set_ylim(0, 10)
ax.xaxis.set_ticks(np.arange(0,nt_model+1,24))
#ax.yaxis.set_ticks(np.linspace(0,10,6))
ax.set_xticklabels(['9/28','9/29','9/30','10/1'])
ax.xaxis.set_minor_locator(AutoMinorLocator(4))
ax.set_xlabel('Time', fontsize=12, fontname='Times New Roman', color='black')
#ylabel = 'Suspended sediment ($\mathregular{mg}$ $\mathregular{l^{-1}}$)'
#ax.set_ylabel(ylabel, fontsize=12, fontname='Times New Roman', color='black')
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(12) for label in labels]
[label.set_color('black') for label in labels]
ax.text(0.05, 0.93, 'd', transform=ax.transAxes, fontsize=16,
        fontname='Times New Roman', fontweight='bold')
ax.tick_params(which='major', direction='in', colors='xkcd:black', length=6, pad=8)
ax.tick_params(which='minor', direction='in', colors='xkcd:black')

plt.tight_layout()
fig.savefig('F12.png', dpi=300)
fig.savefig('F12.pdf', dpi=600)
plt.show()