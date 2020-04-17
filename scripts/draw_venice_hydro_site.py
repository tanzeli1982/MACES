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
from matplotlib.ticker import AutoMinorLocator

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

filename = '/Users/tanz151/Python_maces/src/out_hydro_2002-12-01_2002-12-13.nc'
try:
    nc = Dataset(filename,'r')
    h1_1BF = np.array(nc.variables['h'][0,day0*24:day1*24+1,index1])
    Hwav1_1BF = np.array(nc.variables['Hwav'][0,day0*24:day1*24+1,index1])
    tau1_1BF = np.array(nc.variables['tau'][0,day0*24:day1*24+1,index1])
    h1_2BF = np.array(nc.variables['h'][0,day0*24:day1*24+1,index2])
    Hwav1_2BF = np.array(nc.variables['Hwav'][0,day0*24:day1*24+1,index2])
    tau1_2BF = np.array(nc.variables['tau'][0,day0*24:day1*24+1,index2])
finally:
    nc.close()
h1_1BF = 100 * np.reshape(h1_1BF,(24*(day1-day0)+1))    # cm
Hwav1_1BF = 100 * np.reshape(Hwav1_1BF,(24*(day1-day0)+1))  # cm
tau1_1BF = np.reshape(tau1_1BF,(24*(day1-day0)+1))
h1_2BF = 100 * np.reshape(h1_2BF,(24*(day1-day0)+1))    # cm
Hwav1_2BF = 100 * np.reshape(Hwav1_2BF,(24*(day1-day0)+1))  # cm
tau1_2BF = np.reshape(tau1_2BF,(24*(day1-day0)+1))

nt1_model = np.size(h1_1BF)
tt1_model = np.arange(nt1_model)

day0 = 12
day1 = 15
# read simulation outputs
filename = '/Users/tanz151/Python_maces/src/out_ecogeom_2003-03-21_2003-04-06.nc'
try:
    nc = Dataset(filename,'r')
    x = np.array(nc.variables['x'][0])
    zh = np.array(nc.variables['zh'][0,0,:])
finally:
    nc.close()
    
# find the site index
index1 = np.argmin(np.abs(zh - z_1BF))
index2 = np.argmin(np.abs(zh - z_2BF))

filename = '/Users/tanz151/Python_maces/src/out_hydro_2003-03-21_2003-04-06.nc'
try:
    nc = Dataset(filename,'r')
    h2_1BF = np.array(nc.variables['h'][0,day0*24:day1*24+1,index1])
    Hwav2_1BF = np.array(nc.variables['Hwav'][0,day0*24:day1*24+1,index1])
    tau2_1BF = np.array(nc.variables['tau'][0,day0*24:day1*24+1,index1])
    h2_2BF = np.array(nc.variables['h'][0,day0*24:day1*24+1,index2])
    Hwav2_2BF = np.array(nc.variables['Hwav'][0,day0*24:day1*24+1,index2])
    tau2_2BF = np.array(nc.variables['tau'][0,day0*24:day1*24+1,index2])
finally:
    nc.close()
h2_1BF = 100 * np.reshape(h2_1BF,(24*(day1-day0)+1))    # cm
Hwav2_1BF = 100 * np.reshape(Hwav2_1BF,(24*(day1-day0)+1))  # cm
tau2_1BF = np.reshape(tau2_1BF,(24*(day1-day0)+1))
h2_2BF = 100 * np.reshape(h2_2BF,(24*(day1-day0)+1))    # cm
Hwav2_2BF = 100 * np.reshape(Hwav2_2BF,(24*(day1-day0)+1))  # cm
tau2_2BF = np.reshape(tau2_2BF,(24*(day1-day0)+1))

nt2_model = np.size(h2_1BF)
tt2_model = np.arange(nt2_model)

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
nt1_obs2 = np.size(h_obs_1BF)
tt1_obs2 = np.arange(nt1_obs2)

filename = r'/Users/tanz151/Documents/Projects/TAI_BGC/Data/Hydrodynamics_obs/' + \
    'VeniceLagoon/2BF_OBS.xls'
df = pd.read_excel(filename, sheet_name='2BF', header=None, skiprows=range(3), 
                   usecols='A,B,O,Q')
df.columns = ['Time','Hmo','hw','Turbidity']
Hwav_obs_2BF = 100 * np.array(df['Hmo'])[5319:5511]
h_obs_2BF = 100 * np.array(df['hw'])[5319:5511]
sed_obs_2BF = np.array(df['Turbidity'])[5319:5511]

h_obs_1BF = h_obs_1BF - 100*z_1BF
h_obs_2BF = h_obs_2BF - 100*z_2BF

nt1_obs = np.size(Hwav_obs_1BF)
tt1_obs = np.arange(nt1_obs)/4

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

nt2_obs = np.size(tau_3d_1BF)
tt2_obs = np.arange(nt2_obs)/2

# plot water level, significant wave height, suspended sediment
plt.clf()
fig, axes = plt.subplots(6, 2, figsize=(8,10))

plt.style.use('default')

# 1BF
ax = axes[0][0]
ax.plot(tt1_model, h1_1BF, color='black', linestyle='-', linewidth=2, alpha=0.9)
ax.plot(tt1_obs2, h_obs_1BF, color='C3', linestyle='-', marker='.', markersize=5)
ax.set_xlim(0, nt1_model-1)
ax.set_ylim(75, 160)
ax.xaxis.set_ticks(np.arange(0,nt1_model,24))
ax.yaxis.set_ticks(np.linspace(80,160,5))
ax.set_xticklabels(['12/10','12/11','12/12'])
ax.xaxis.set_minor_locator(AutoMinorLocator(4))
ax.set_title('1BF', fontsize=14, fontname='Times New Roman', color='black')
ylabel = 'Water depth\n($\mathregular{cm}$)'
ax.set_ylabel(ylabel, fontsize=11, fontname='Times New Roman', color='black')
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(11) for label in labels]
[label.set_color('black') for label in labels]
ax.tick_params(which='major', direction='in', colors='xkcd:black', length=6, pad=8)
ax.tick_params(which='minor', direction='in', colors='xkcd:black')

ax = axes[1][0]
ax.plot(tt1_model, Hwav1_1BF, color='black', linestyle='-', linewidth=2, alpha=0.9)
ax.plot(tt1_obs, Hwav_obs_1BF, color='C3', linestyle='-', marker='.', markersize=5)
ax.set_xlim(0, nt1_model-1)
ax.set_ylim(0, 50)
ax.xaxis.set_ticks(np.arange(0,nt1_model,24))
ax.yaxis.set_ticks(np.linspace(0,50,6))
ax.set_xticklabels(['12/10','12/11','12/12'])
ax.xaxis.set_minor_locator(AutoMinorLocator(4))
ylabel = 'Significant wave\nheight ($\mathregular{cm}$)'
ax.set_ylabel(ylabel, fontsize=11, fontname='Times New Roman', color='black')
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(11) for label in labels]
[label.set_color('black') for label in labels]
ax.tick_params(which='major', direction='in', colors='xkcd:black', length=6, pad=8)
ax.tick_params(which='minor', direction='in', colors='xkcd:black')

ax = axes[2][0]
ax.plot(tt1_model, tau1_1BF, color='black', linestyle='-', linewidth=2, alpha=0.9)
ax.set_xlim(0, nt1_model-1)
ax.set_ylim(0, 0.4)
ax.xaxis.set_ticks(np.arange(0,nt1_model,24))
ax.yaxis.set_ticks(np.linspace(0,0.4,5))
ax.set_xticklabels(['12/10','12/11','12/12'])
ax.xaxis.set_minor_locator(AutoMinorLocator(4))
ylabel = 'Bottom shear\nstress ($\mathregular{Pa}$)'
ax.set_ylabel(ylabel, fontsize=11, fontname='Times New Roman', color='black')
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(11) for label in labels]
[label.set_color('black') for label in labels]
ax.tick_params(which='major', direction='in', colors='xkcd:black', length=6, pad=8)
ax.tick_params(which='minor', direction='in', colors='xkcd:black')

ax = axes[3][0]
ax.plot(tt2_model, h2_1BF, color='black', linestyle='-', linewidth=2, alpha=0.9)
ax.set_xlim(0, nt2_model-1)
ax.set_ylim(60, 180)
ax.xaxis.set_ticks(np.arange(0,nt2_model,24))
ax.yaxis.set_ticks(np.linspace(60,180,5))
ax.set_xticklabels(['4/2','4/3','4/4','4/5'])
ax.xaxis.set_minor_locator(AutoMinorLocator(4))
ylabel = 'Water depth\n($\mathregular{cm}$)'
ax.set_ylabel(ylabel, fontsize=11, fontname='Times New Roman', color='black')
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(11) for label in labels]
[label.set_color('black') for label in labels]
ax.tick_params(which='major', direction='in', colors='xkcd:black', length=6, pad=8)
ax.tick_params(which='minor', direction='in', colors='xkcd:black')

ax = axes[4][0]
ax.plot(tt2_model, Hwav2_1BF, color='black', linestyle='-', linewidth=2, alpha=0.9)
ax.set_xlim(0, nt2_model-1)
ax.set_ylim(0, 60)
ax.xaxis.set_ticks(np.arange(0,nt2_model,24))
ax.yaxis.set_ticks(np.linspace(0,60,5))
ax.set_xticklabels(['4/2','4/3','4/4','4/5'])
ax.xaxis.set_minor_locator(AutoMinorLocator(4))
ylabel = 'Significant wave\nheight ($\mathregular{cm}$)'
ax.set_ylabel(ylabel, fontsize=11, fontname='Times New Roman', color='black')
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(11) for label in labels]
[label.set_color('black') for label in labels]
ax.tick_params(which='major', direction='in', colors='xkcd:black', length=6, pad=8)
ax.tick_params(which='minor', direction='in', colors='xkcd:black')

ax = axes[5][0]
ax.plot(tt2_model, tau2_1BF, color='black', linestyle='-', linewidth=2, alpha=0.9)
ax.plot(tt2_obs, tau_3d_1BF, color='C3', linestyle='-', marker='.', markersize=5)
ax.set_xlim(0, nt2_model-1)
ax.set_ylim(0, 1)
ax.xaxis.set_ticks(np.arange(0,nt2_model,24))
ax.yaxis.set_ticks(np.linspace(0,1,6))
ax.set_xticklabels(['4/2','4/3','4/4','4/5'])
ax.xaxis.set_minor_locator(AutoMinorLocator(4))
ax.set_xlabel('Time', fontsize=11, fontname='Times New Roman', color='black')
ylabel = 'Bottom shear\nstress ($\mathregular{Pa}$)'
ax.set_ylabel(ylabel, fontsize=11, fontname='Times New Roman', color='black')
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(11) for label in labels]
[label.set_color('black') for label in labels]
ax.tick_params(which='major', direction='in', colors='xkcd:black', length=6, pad=8)
ax.tick_params(which='minor', direction='in', colors='xkcd:black')

# 2BF
ax = axes[0][1]
ax.plot(tt1_model, h1_2BF, color='black', linestyle='-', linewidth=2, alpha=0.9)
ax.plot(tt1_obs, h_obs_2BF, color='C3', linestyle='-', marker='.', markersize=5)
ax.set_xlim(0, nt1_model-1)
ax.set_ylim(175, 260)
ax.xaxis.set_ticks(np.arange(0,nt1_model+1,24))
ax.yaxis.set_ticks(np.linspace(180,260,5))
ax.set_xticklabels(['12/10','12/11','12/12'])
ax.xaxis.set_minor_locator(AutoMinorLocator(4))
ax.set_title('2BF', fontsize=14, fontname='Times New Roman', color='black')
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(11) for label in labels]
[label.set_color('black') for label in labels]
ax.tick_params(which='major', direction='in', colors='xkcd:black', length=6, pad=8)
ax.tick_params(which='minor', direction='in', colors='xkcd:black')

ax = axes[1][1]
ax.plot(tt1_model, Hwav1_2BF, color='black', linestyle='-', linewidth=2, alpha=0.9)
ax.plot(tt1_obs, Hwav_obs_2BF, color='C3', linestyle='-', marker='.', markersize=5)
ax.set_xlim(0, nt1_model-1)
ax.set_ylim(0, 50)
ax.xaxis.set_ticks(np.arange(0,nt1_model+1,24))
ax.yaxis.set_ticks(np.linspace(0,50,6))
ax.set_xticklabels(['12/10','12/11','12/12'])
ax.xaxis.set_minor_locator(AutoMinorLocator(4))
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(11) for label in labels]
[label.set_color('black') for label in labels]
ax.tick_params(which='major', direction='in', colors='xkcd:black', length=6, pad=8)
ax.tick_params(which='minor', direction='in', colors='xkcd:black')

ax = axes[2][1]
ax.plot(tt1_model, tau1_2BF, color='black', linestyle='-', linewidth=2, alpha=0.9)
ax.set_xlim(0, nt1_model-1)
ax.set_ylim(0, 0.4)
ax.xaxis.set_ticks(np.arange(0,nt1_model,24))
ax.yaxis.set_ticks(np.linspace(0,0.4,5))
ax.set_xticklabels(['12/10','12/11','12/12'])
ax.xaxis.set_minor_locator(AutoMinorLocator(4))
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(11) for label in labels]
[label.set_color('black') for label in labels]
ax.tick_params(which='major', direction='in', colors='xkcd:black', length=6, pad=8)
ax.tick_params(which='minor', direction='in', colors='xkcd:black')

ax = axes[3][1]
ax.plot(tt2_model, h2_2BF, color='black', linestyle='-', linewidth=2, alpha=0.9)
ax.set_xlim(0, nt2_model-1)
ax.set_ylim(160, 280)
ax.xaxis.set_ticks(np.arange(0,nt2_model,24))
ax.yaxis.set_ticks(np.linspace(160,280,5))
ax.set_xticklabels(['4/2','4/3','4/4','4/5'])
ax.xaxis.set_minor_locator(AutoMinorLocator(4))
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(11) for label in labels]
[label.set_color('black') for label in labels]
ax.tick_params(which='major', direction='in', colors='xkcd:black', length=6, pad=8)
ax.tick_params(which='minor', direction='in', colors='xkcd:black')

ax = axes[4][1]
ax.plot(tt2_model, Hwav2_2BF, color='black', linestyle='-', linewidth=2, alpha=0.9)
ax.set_xlim(0, nt2_model-1)
ax.set_ylim(0, 60)
ax.xaxis.set_ticks(np.arange(0,nt2_model,24))
ax.yaxis.set_ticks(np.linspace(0,60,5))
ax.set_xticklabels(['4/2','4/3','4/4','4/5'])
ax.xaxis.set_minor_locator(AutoMinorLocator(4))
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(11) for label in labels]
[label.set_color('black') for label in labels]
ax.tick_params(which='major', direction='in', colors='xkcd:black', length=6, pad=8)
ax.tick_params(which='minor', direction='in', colors='xkcd:black')

ax = axes[5][1]
ax.plot(tt2_model, tau2_2BF, color='black', linestyle='-', linewidth=2, alpha=0.9)
ax.plot(tt2_obs, tau_3d_2BF, color='C3', linestyle='-', marker='.', markersize=5)
ax.set_xlim(0, nt2_model-1)
ax.set_ylim(0, 1)
ax.xaxis.set_ticks(np.arange(0,nt2_model,24))
ax.yaxis.set_ticks(np.linspace(0,1,6))
ax.set_xticklabels(['4/2','4/3','4/4','4/5'])
ax.xaxis.set_minor_locator(AutoMinorLocator(4))
ax.set_xlabel('Time', fontsize=11, fontname='Times New Roman', color='black')
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(11) for label in labels]
[label.set_color('black') for label in labels]
ax.tick_params(which='major', direction='in', colors='xkcd:black', length=6, pad=8)
ax.tick_params(which='minor', direction='in', colors='xkcd:black')

plt.tight_layout()
fig.savefig('F4.png', dpi=300)
fig.savefig('F4.pdf', dpi=600)
plt.show()