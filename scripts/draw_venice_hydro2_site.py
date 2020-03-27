#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 26 16:10:42 2020

Compare the simulated bottom shear stress with benchmark estimates

@author: Zeli Tan
"""

import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset

z_1BF = -1.1    # 1BF is -1.1 m and 2BF is -2.1 m
z_2BF = -2.1
day0 = 456
day1 = 459

# read simulation outputs
filename = '/Users/tanz151/Python_maces/src/out_ecogeom_2002-01-01_2004-01-01.nc'
try:
    nc = Dataset(filename,'r')
    x = np.array(nc.variables['x'][:][0])
    zh = np.array(nc.variables['zh'][:][0,0,:])
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
    U_1BF = np.array(nc.variables['U'][0,day0*24:day1*24+1,index1])
    h_2BF = np.array(nc.variables['h'][0,day0*24:day1*24+1,index2])
    Hwav_2BF = np.array(nc.variables['Hwav'][0,day0*24:day1*24+1,index2])
    sed_2BF = np.array(nc.variables['TSM'][0,day0*24:day1*24+1,index2])
    tau_2BF = np.array(nc.variables['tau'][0,day0*24:day1*24+1,index2])
    U_2BF = np.array(nc.variables['U'][0,day0*24:day1*24+1,index2])
finally:
    nc.close()
h_1BF = 100 * np.reshape(h_1BF,(24*(day1-day0)+1))    # cm
Hwav_1BF = 100 * np.reshape(Hwav_1BF,(24*(day1-day0)+1))  # cm
sed_1BF = 1e3 * np.reshape(sed_1BF,(24*(day1-day0)+1))    # mg/l
tau_1BF = np.reshape(tau_1BF,(24*(day1-day0)+1))
U_1BF = np.reshape(U_1BF,(24*(day1-day0)+1))
h_2BF = 100 * np.reshape(h_2BF,(24*(day1-day0)+1))    # cm
Hwav_2BF = 100 * np.reshape(Hwav_2BF,(24*(day1-day0)+1))  # cm
sed_2BF = 1e3 * np.reshape(sed_2BF,(24*(day1-day0)+1))    # mg/l
tau_2BF = np.reshape(tau_2BF,(24*(day1-day0)+1))
U_2BF = np.reshape(U_2BF,(24*(day1-day0)+1))

nt_model = np.size(tau_1BF)
tt_model = np.arange(nt_model)

# read data
filename = r'/Users/tanz151/Documents/Projects/TAI_BGC/Data/Hydrodynamics_obs/' + \
    r'VeniceLagoon/2-4Apr03-TauMax1BF&2BF.txt'
tau_3d_1BF = []
tau_3d_2BF = []
try:
    f = open(filename, 'r')
    f.readline()    # skip header
    for line in f:
        line = line.strip()
        columns = line.split()
        tau_3d_1BF.append(float(columns[1]))
        tau_3d_2BF.append(float(columns[2]))
finally:
    f.close()
tau_3d_1BF = np.array(tau_3d_1BF)
tau_3d_2BF = np.array(tau_3d_2BF)

nt_obs = np.size(tau_3d_1BF)
tt_obs = np.arange(nt_obs)/2

# plot
plt.clf()
fig, axes = plt.subplots(4, 2, figsize=(8,10))

plt.style.use('default')

# 1BF
ax = axes[0][0]
ax.plot(tt_model, h_1BF, color='black', linestyle='-', linewidth=2, alpha=0.9)
ax.set_xlim(0, nt_model-1)
#ax.set_ylim(-100, 100)
ax.xaxis.set_ticks(np.arange(0,nt_model,24))
#ax.yaxis.set_ticks(np.linspace(-100,100,5))
ax.set_xticklabels(['4/2','4/3','4/4','4/5'])
#ax.xaxis.set_minor_locator(AutoMinorLocator(24))
ylabel = 'Water depth ($\mathregular{cm}$)'
ax.set_ylabel(ylabel, fontsize=12, fontname='Times New Roman', color='black')
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(12) for label in labels]
[label.set_color('black') for label in labels]

ax = axes[1][0]
ax.plot(tt_model, Hwav_1BF, color='black', linestyle='-', linewidth=2, alpha=0.9)
ax.set_xlim(0, nt_model-1)
#ax.set_ylim(0, 100)
ax.xaxis.set_ticks(np.arange(0,nt_model,24))
#ax.yaxis.set_ticks(np.linspace(0,100,6))
ax.set_xticklabels(['4/2','4/3','4/4','4/5'])
#ax.xaxis.set_minor_locator(AutoMinorLocator(24))
ylabel = 'Significant wave height ($\mathregular{cm}$)'
ax.set_ylabel(ylabel, fontsize=12, fontname='Times New Roman', color='black')
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(12) for label in labels]
[label.set_color('black') for label in labels]

ax = axes[2][0]
ax.plot(tt_model, sed_1BF, color='black', linestyle='-', linewidth=2, alpha=0.9)
ax.set_xlim(0, nt_model-1)
#ax.set_ylim(0, 150)
ax.xaxis.set_ticks(np.arange(0,nt_model,24))
#ax.yaxis.set_ticks(np.linspace(0,150,6))
ax.set_xticklabels(['4/2','4/3','4/4','4/5'])
#ax.xaxis.set_minor_locator(AutoMinorLocator(24))
ylabel = 'Suspended sediment ($\mathregular{mg}$ $\mathregular{{l}^{-1}}$)'
ax.set_ylabel(ylabel, fontsize=12, fontname='Times New Roman', color='black')
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(12) for label in labels]
[label.set_color('black') for label in labels]

ax = axes[3][0]
ax.plot(tt_model, tau_1BF, color='black', linestyle='-', linewidth=2, alpha=0.9)
ax.plot(tt_obs, tau_3d_1BF, color='C3', linestyle='-', marker='.', markersize=5)
ax.set_xlim(0, nt_model-1)
#ax.set_ylim(0, 1)
ax.xaxis.set_ticks(np.arange(0,nt_model,24))
#ax.yaxis.set_ticks(np.linspace(0,1,6))
ax.set_xticklabels(['4/2','4/3','4/4','4/5'])
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
ax.set_xlim(0, nt_model-1)
#ax.set_ylim(-100, 100)
ax.xaxis.set_ticks(np.arange(0,nt_model+1,24))
#ax.yaxis.set_ticks(np.linspace(-100,100,5))
ax.set_xticklabels(['4/2','4/3','4/4','4/5'])
#ax.xaxis.set_minor_locator(AutoMinorLocator(24))
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(12) for label in labels]
[label.set_color('black') for label in labels]

ax = axes[1][1]
ax.plot(tt_model, Hwav_2BF, color='black', linestyle='-', linewidth=2, alpha=0.9)
ax.set_xlim(0, nt_model-1)
#ax.set_ylim(0, 100)
ax.xaxis.set_ticks(np.arange(0,nt_model+1,24))
#ax.yaxis.set_ticks(np.linspace(0,100,6))
ax.set_xticklabels(['4/2','4/3','4/4','4/5'])
#ax.xaxis.set_minor_locator(AutoMinorLocator(24))
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(12) for label in labels]
[label.set_color('black') for label in labels]

ax = axes[2][1]
ax.plot(tt_model, sed_2BF, color='black', linestyle='-', linewidth=2, alpha=0.9)
ax.set_xlim(0, nt_model-1)
#ax.set_ylim(0, 150)
ax.xaxis.set_ticks(np.arange(0,nt_model+1,24))
#ax.yaxis.set_ticks(np.linspace(0,150,6))
ax.set_xticklabels(['4/2','4/3','4/4','4/5'])
#ax.xaxis.set_minor_locator(AutoMinorLocator(24))
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(12) for label in labels]
[label.set_color('black') for label in labels]

ax = axes[3][1]
ax.plot(tt_model, tau_2BF, color='black', linestyle='-', linewidth=2, alpha=0.9)
ax.plot(tt_obs, tau_3d_2BF, color='C3', linestyle='-', marker='.', markersize=5)
ax.set_xlim(0, nt_model-1)
#ax.set_ylim(0, 1)
ax.xaxis.set_ticks(np.arange(0,nt_model,24))
#ax.yaxis.set_ticks(np.linspace(0,1,6))
ax.set_xticklabels(['4/2','4/3','4/4','4/5'])
#ax.xaxis.set_minor_locator(AutoMinorLocator(24))
ylabel = 'Bottom shear stress ($\mathregular{Pa}$)'
ax.set_ylabel(ylabel, fontsize=12, fontname='Times New Roman', color='black')
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(12) for label in labels]
[label.set_color('black') for label in labels]

plt.tight_layout()
fig.savefig('venice_hydro2_sim.png', dpi=300)
#fig.savefig('venice_hydro2_sim.pdf', dpi=600)
plt.show()