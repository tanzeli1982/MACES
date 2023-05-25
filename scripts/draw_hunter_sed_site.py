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

models = ['F06', 'T03', 'KM12', 'M12', 'F07', 'VDK05', 'DA07']
casemodel = 'M12'

day0 = (date(2004,9,28) - date(2004,9,25)).days
day1 = (date(2004,10,1) - date(2004,9,25)).days
# read simulated water depth and suspended sediment
rdir = '/Users/tanz151/Documents/Projects/TAI_BGC/Drafts/Outputs/HunterEstuary/'
filename = rdir + 'maces_ecogeom_2004-09-25_2004-10-06_F06_11099.nc'
try:
    nc = Dataset(filename,'r')
    x = np.array(nc.variables['x'][:])
    zh = np.array(nc.variables['zh'][0,:])
finally:
    nc.close()
    
# find the site index
index1 = np.argmin(np.abs(zh - z1))
index2 = np.argmin(np.abs(zh - z2))
index3 = np.argmin(np.abs(zh - z3))
index4 = np.argmin(np.abs(zh - z4))

print(index1, index2, index3, index4)

sed1_sim = {}
sed2_sim = {}
sed3_sim = {}
sed4_sim = {}
tau1_sim = {}
tau2_sim = {}
tau3_sim = {}
tau4_sim = {}
for model in models:
    filename = rdir + 'maces_hydro_2004-09-25_2004-10-06_' + model + '_11099.nc'
    try:
        nc = Dataset(filename, 'r')
        sed1 = 1e3*np.array(nc.variables['TSM'][day0*24:day1*24,index1])
        sed2 = 1e3*np.array(nc.variables['TSM'][day0*24:day1*24,index2])
        sed3 = 1e3*np.array(nc.variables['TSM'][day0*24:day1*24,index3])
        sed4 = 1e3*np.array(nc.variables['TSM'][day0*24:day1*24,index4])
        tau1 = 1e2*np.array(nc.variables['tau'][day0*24:day1*24,index1])
        tau2 = 1e2*np.array(nc.variables['tau'][day0*24:day1*24,index2])
        tau3 = 1e2*np.array(nc.variables['tau'][day0*24:day1*24,index3])
        tau4 = 1e2*np.array(nc.variables['tau'][day0*24:day1*24,index4])
    finally:
        nc.close()
    sed1_sim[model] = sed1
    sed2_sim[model] = sed2
    sed3_sim[model] = sed3
    sed4_sim[model] = sed4
    tau1_sim[model] = tau1
    tau2_sim[model] = tau2
    tau3_sim[model] = tau3
    tau4_sim[model] = tau4
nt_model = np.size(sed1_sim['F06'])
tt_model = np.arange(nt_model)

Bag_sim = {}
for model in models:
    filename = rdir + 'maces_ecogeom_2004-09-25_2004-10-06_' + model + '_11099.nc'
    try:
        nc = Dataset(filename,'r')
        Bag = 1e3 * np.mean(np.array(nc.variables['Bag'][day0:day1,:]),axis=0)
    finally:
        nc.close()
    Bag_sim[model] = Bag

# read benchmark estimates
filename = rdir + 'SSC_benchmark_2004-09-28_2004-09-30.nc'
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
        markersize=10)
handles = []
for key in models:
    indx = len(handles)
    h, = ax.plot(tt_model, sed1_sim[key], color=colors[indx], 
                 linestyle=linestyles[indx], linewidth=2, alpha=1)
    handles.append(h)
ax.plot(tt_model, tau1_sim[casemodel], color='gray', linestyle='-', linewidth=1, alpha=0.9)
ax.set_xlim(0, nt_model)
ax.set_ylim(0, 40)
ax.xaxis.set_ticks(np.arange(0,nt_model+1,24))
ax.yaxis.set_ticks(np.linspace(0,40,5))
ax.set_xticklabels(['9/28','9/29','9/30','10/1'])
ax.xaxis.set_minor_locator(AutoMinorLocator(4))
#ax.set_xlabel('Time', fontsize=11, fontname='Times New Roman', color='black')
ylabel = 'Suspended sediment ($\mathregular{mg}$ $\mathregular{l^{-1}}$)'
ax.set_ylabel(ylabel, fontsize=12, fontname='Times New Roman', color='black')
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(12) for label in labels]
[label.set_color('black') for label in labels]
ax.text(0.05, 0.9, 'a', transform=ax.transAxes, fontsize=16,
        fontname='Times New Roman', fontweight='bold')
ax.tick_params(which='major', direction='in', colors='xkcd:black', length=6, pad=8)
ax.tick_params(which='minor', direction='in', colors='xkcd:black')

rmse_sed1 = {}
for key in models:
    rmse_sed1[key] = 0.0
    count = 0
    for ii, tt in enumerate(tt_obs):
        if tt in tt_model and np.isfinite(sed1_obs[ii]):
            index = np.where(tt_model==tt)
            if sed1_sim[key][index]<1e-3 and sed1_obs[ii]>10:
                continue
            rmse_sed1[key] = rmse_sed1[key] + (sed1_sim[key][index] - sed1_obs[ii])**2
            count = count + 1
    rmse_sed1[key] = np.sqrt( rmse_sed1[key] / count )

# mangrove edge
ax = axes[0][1]
ax.plot(tt_obs, sed2_obs, color='black', linestyle='-', linewidth=2, marker='.', 
        markersize=10)
handles = []
for key in sed2_sim:
    indx = len(handles)
    h, = ax.plot(tt_model, sed2_sim[key], color=colors[indx], 
                 linestyle=linestyles[indx], linewidth=2, alpha=1)
    handles.append(h)
ax.plot(tt_model, tau2_sim[casemodel], color='gray', linestyle='-', linewidth=1, alpha=0.9)
ax.set_xlim(0, nt_model)
ax.set_ylim(0, 40)
ax.xaxis.set_ticks(np.arange(0,nt_model+1,24))
ax.yaxis.set_ticks(np.linspace(0,40,5))
ax.set_xticklabels(['9/28','9/29','9/30','10/1'])
ax.xaxis.set_minor_locator(AutoMinorLocator(4))
#ax.set_xlabel('Time', fontsize=11, fontname='Times New Roman', color='black')
#ylabel = 'Suspended sediment ($\mathregular{mg}$ $\mathregular{l^{-1}}$)'
#ax.set_ylabel(ylabel, fontsize=12, fontname='Times New Roman', color='black')
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(12) for label in labels]
[label.set_color('black') for label in labels]
ax.text(0.05, 0.9, 'b', transform=ax.transAxes, fontsize=16,
        fontname='Times New Roman', fontweight='bold')
ax.tick_params(which='major', direction='in', colors='xkcd:black', length=6, pad=8)
ax.tick_params(which='minor', direction='in', colors='xkcd:black')

rmse_sed2 = {}
for key in models:
    rmse_sed2[key] = 0.0
    count = 0
    for ii, tt in enumerate(tt_obs):
        if tt in tt_model and np.isfinite(sed2_obs[ii]):
            index = np.where(tt_model==tt)
            if sed2_sim[key][index]<1e-3 and sed2_obs[ii]>10:
                continue
            rmse_sed2[key] = rmse_sed2[key] + (sed1_sim[key][index] - sed1_obs[ii])**2
            count = count + 1
    rmse_sed2[key] = np.sqrt( rmse_sed2[key] / count )

# mangrove interior
ax = axes[1][0]
ax.plot(tt_obs, sed3_obs, color='black', linestyle='-', linewidth=2, marker='.', 
        markersize=10)
handles = []
for key in sed3_sim:
    indx = len(handles)
    h, = ax.plot(tt_model, sed3_sim[key], color=colors[indx], 
                 linestyle=linestyles[indx], linewidth=2, alpha=1)
    handles.append(h)
ax.plot(tt_model, tau3_sim[casemodel], color='gray', linestyle='-', linewidth=1, alpha=0.9)
ax.set_xlim(0, nt_model)
ax.set_ylim(0, 40)
ax.xaxis.set_ticks(np.arange(0,nt_model+1,24))
ax.yaxis.set_ticks(np.linspace(0,40,5))
ax.set_xticklabels(['9/28','9/29','9/30','10/1'])
ax.xaxis.set_minor_locator(AutoMinorLocator(4))
ax.set_xlabel('Time', fontsize=12, fontname='Times New Roman', color='black')
ylabel = 'Suspended sediment ($\mathregular{mg}$ $\mathregular{l^{-1}}$)'
ax.set_ylabel(ylabel, fontsize=12, fontname='Times New Roman', color='black')
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(12) for label in labels]
[label.set_color('black') for label in labels]
ax.text(0.05, 0.9, 'c', transform=ax.transAxes, fontsize=16,
        fontname='Times New Roman', fontweight='bold')
ax.tick_params(which='major', direction='in', colors='xkcd:black', length=6, pad=8)
ax.tick_params(which='minor', direction='in', colors='xkcd:black')

rmse_sed3 = {}
for key in models:
    rmse_sed3[key] = 0.0
    count = 0
    for ii, tt in enumerate(tt_obs):
        if tt in tt_model and np.isfinite(sed3_obs[ii]):
            index = np.where(tt_model==tt)
            if sed3_sim[key][index]<1e-3 and sed3_obs[ii]>10:
                continue
            rmse_sed3[key] = rmse_sed3[key] + (sed3_sim[key][index] - sed3_obs[ii])**2
            count = count + 1
    rmse_sed3[key] = np.sqrt( rmse_sed3[key] / count )

# marsh edge
ax = axes[1][1]
ax.plot(tt_obs, sed4_obs, color='black', linestyle='-', linewidth=2, marker='.', 
        markersize=10)
handles = []
for key in sed4_sim:
    indx = len(handles)
    h, = ax.plot(tt_model, sed4_sim[key], color=colors[indx], 
                 linestyle=linestyles[indx], linewidth=2, alpha=1)
    handles.append(h)
ax.plot(tt_model, tau4_sim[casemodel], color='gray', linestyle='-', linewidth=1, alpha=0.9)
legend = ax.legend(handles, models, numpoints=1, loc='upper right', 
                   prop={'family':'Times New Roman', 'size':'large'}, 
                   framealpha=0.0, ncol=2)
ax.set_xlim(0, nt_model)
ax.set_ylim(0, 40)
ax.xaxis.set_ticks(np.arange(0,nt_model+1,24))
ax.yaxis.set_ticks(np.linspace(0,40,5))
ax.set_xticklabels(['9/28','9/29','9/30','10/1'])
ax.xaxis.set_minor_locator(AutoMinorLocator(4))
ax.set_xlabel('Time', fontsize=12, fontname='Times New Roman', color='black')
#ylabel = 'Suspended sediment ($\mathregular{mg}$ $\mathregular{l^{-1}}$)'
#ax.set_ylabel(ylabel, fontsize=12, fontname='Times New Roman', color='black')
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(12) for label in labels]
[label.set_color('black') for label in labels]
ax.text(0.05, 0.9, 'd', transform=ax.transAxes, fontsize=16,
        fontname='Times New Roman', fontweight='bold')
ax.tick_params(which='major', direction='in', colors='xkcd:black', length=6, pad=8)
ax.tick_params(which='minor', direction='in', colors='xkcd:black')

rmse_sed4 = {}
for key in models:
    rmse_sed4[key] = 0.0
    count = 0
    for ii, tt in enumerate(tt_obs):
        if tt in tt_model and np.isfinite(sed4_obs[ii]):
            index = np.where(tt_model==tt)
            if sed4_sim[key][index]<1e-3 and sed4_obs[ii]>10:
                continue
            rmse_sed4[key] = rmse_sed4[key] + (sed4_sim[key][index] - sed4_obs[ii])**2
            count = count + 1
    rmse_sed4[key] = np.sqrt( rmse_sed4[key] / count )
    
print('RMSE ' + casemodel + ': ', rmse_sed1[casemodel], rmse_sed2[casemodel], 
      rmse_sed3[casemodel], rmse_sed4[casemodel])

plt.tight_layout()
fig.savefig('F9.png', dpi=300)
#fig.savefig('F9.jpg', dpi=600)
plt.show()