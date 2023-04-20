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
import matplotlib.gridspec as gridspec
from netCDF4 import Dataset
from matplotlib.ticker import AutoMinorLocator

z_1BF = -1.1    # 1BF is -1.1 m and 2BF is -2.1 m
z_2BF = -2.1

# read forcing
t1 = {'d0': 343, 'd1': 345}
t2 = {'d0': 456, 'd1': 459}

rdir = '/Users/tanz151/Documents/Projects/TAI_BGC/Drafts/Outputs/VeniceLagoon/'

filename = rdir + 'force_h.nc'
try:
    nc = Dataset(filename,'r')
    h_var1 = 100 * np.array(nc.variables['h'][t1['d0']*144:t1['d1']*144+1]) # cm
    h_var2 = 100 * np.array(nc.variables['h'][t2['d0']*144:t2['d1']*144+1]) # cm
finally:
    nc.close()
    
h1b_1BF = h_var1 - z_1BF*100
h1b_2BF = h_var1 - z_2BF*100
h2b_1BF = h_var2 - z_1BF*100
h2b_2BF = h_var2 - z_2BF*100
nt1_h = len(h_var1)
tt1_h = np.arange(nt1_h)
nt2_h = len(h_var2)
tt2_h = np.arange(nt2_h)

filename = rdir + 'force_U10.nc'
try:
    nc = Dataset(filename,'r')
    U10_var1 = np.array(nc.variables['U10'][t1['d0']*96:t1['d1']*96+1])
    U10_var2 = np.array(nc.variables['U10'][t2['d0']*96:t2['d1']*96+1])
finally:
    nc.close()

nt1_U10 = np.size(U10_var1)
tt1_U10 = np.arange(nt1_U10)
nt2_U10 = np.size(U10_var2)
tt2_U10 = np.arange(nt2_U10)

# read simulation outputs
day0 = 9
day1 = 11

filename = rdir + 'maces_ecogeom_2002-12-01_2002-12-13_466.nc'
try:
    nc = Dataset(filename,'r')
    x = np.array(nc.variables['x'][:])
    zh = np.array(nc.variables['zh'][0,:])
finally:
    nc.close()
    
# find the site index
index1 = np.argmin(np.abs(zh - z_1BF))
index2 = np.argmin(np.abs(zh - z_2BF))

filename = rdir + 'maces_hydro_2002-12-01_2002-12-13_466.nc'
try:
    nc = Dataset(filename,'r')
    h1_1BF = np.array(nc.variables['h'][day0*24:day1*24+1,index1])
    Hwav1_1BF = np.array(nc.variables['Hwav'][day0*24:day1*24+1,index1])
    tau1_1BF = np.array(nc.variables['tau'][day0*24:day1*24+1,index1])
    h1_2BF = np.array(nc.variables['h'][day0*24:day1*24+1,index2])
    Hwav1_2BF = np.array(nc.variables['Hwav'][day0*24:day1*24+1,index2])
    tau1_2BF = np.array(nc.variables['tau'][day0*24:day1*24+1,index2])
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
filename = rdir + 'maces_ecogeom_2003-03-21_2003-04-06_466.nc'
try:
    nc = Dataset(filename,'r')
    x = np.array(nc.variables['x'][:])
    zh = np.array(nc.variables['zh'][0,:])
finally:
    nc.close()
    
# find the site index
index1 = np.argmin(np.abs(zh - z_1BF))
index2 = np.argmin(np.abs(zh - z_2BF))

filename = rdir + 'maces_hydro_2003-03-21_2003-04-06_466.nc'
try:
    nc = Dataset(filename,'r')
    h2_1BF = np.array(nc.variables['h'][day0*24:day1*24+1,index1])
    Hwav2_1BF = np.array(nc.variables['Hwav'][day0*24:day1*24+1,index1])
    tau2_1BF = np.array(nc.variables['tau'][day0*24:day1*24+1,index1])
    h2_2BF = np.array(nc.variables['h'][day0*24:day1*24+1,index2])
    Hwav2_2BF = np.array(nc.variables['Hwav'][day0*24:day1*24+1,index2])
    tau2_2BF = np.array(nc.variables['tau'][day0*24:day1*24+1,index2])
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
filename = rdir + '1BF_OBS.xls'
df = pd.read_excel(filename, sheet_name='1BF', header=None, skiprows=range(3), 
                   usecols='A,B,F,O,Q')
df.columns = ['Time','Hmo','Hmax','hw','Turbidity']
Hwav_obs_1BF = 100 * np.array(df['Hmo'])[5334:5526] # cm
sed_obs_1BF = np.array(df['Turbidity'])[5334:5526]  # mg/l

filename = rdir + 'WaterLevelClose1BF.xls'
df = pd.read_excel(filename, sheet_name='Valori orari 2002', header=None, 
                   skiprows=range(4), usecols='A:C')
df.columns = ['Date','Hour','hw']
h_obs_1BF = 100 * np.array(df['hw'])[8231:8303]
nt1_obs2 = np.size(h_obs_1BF)
tt1_obs2 = np.arange(nt1_obs2)

filename = rdir + '2BF_OBS.xls'
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

filename = rdir + '2-4Apr03-TauMax1BF&2BF.txt'
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
fig = plt.figure(figsize=(8,10))

gs = gridspec.GridSpec(nrows=6, ncols=2)

plt.style.use('default')

# 1BF
ax = fig.add_subplot(gs[0,0])
ax.plot(tt1_model, h1_1BF, color='black', linestyle='-', linewidth=2, alpha=0.9)
ax.plot(tt1_h, h1b_1BF, color='black', linestyle='--', linewidth=1, alpha=0.9)
ax.plot(tt1_obs2, h_obs_1BF, color='C3', linestyle='-', marker='.', markersize=5)
ax.set_xlim(0, nt1_model-1)
ax.set_ylim(75, 160)
ax.xaxis.set_ticks(np.arange(0,nt1_model,24))
ax.yaxis.set_ticks(np.linspace(80,160,5))
ax.set_xticklabels(['12/10','12/11','12/12'])
ax.xaxis.set_minor_locator(AutoMinorLocator(4))
ax.set_title('1BF', fontsize=14, fontname='Times New Roman', fontweight='bold')
ylabel = 'Water depth\n($\mathregular{cm}$)'
ax.set_ylabel(ylabel, fontsize=11, fontname='Times New Roman', fontweight='bold')
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(11) for label in labels]
[label.set_color('black') for label in labels]
[label.set_fontweight('bold') for label in labels]
ax.tick_params(which='major', direction='in', colors='xkcd:black', length=6, pad=4)
ax.tick_params(which='minor', direction='in', colors='xkcd:black')

rmse_h_1BF = 0.0
count = 0
for ii, tt in enumerate(tt1_obs2):
    if tt in tt1_model and np.isfinite(h_obs_1BF[ii]):
        indx = np.where(tt1_model==tt)
        rmse_h_1BF = rmse_h_1BF + (h1_1BF[indx] - h_obs_1BF[ii])**2
        count = count + 1
rmse_h_1BF = np.sqrt( rmse_h_1BF / count )
nrmse_h_1BF = rmse_h_1BF / np.nanmean(h_obs_1BF)

ax = fig.add_subplot(gs[1,0])
ax.plot(tt1_model, Hwav1_1BF, color='black', linestyle='-', linewidth=2, alpha=0.9)
ax.plot(tt1_obs, Hwav_obs_1BF, color='C3', linestyle='-', marker='.', markersize=5)
ax.set_xlim(0, nt1_model-1)
ax.set_ylim(0, 50)
ax.xaxis.set_ticks(np.arange(0,nt1_model,24))
ax.yaxis.set_ticks(np.linspace(0,50,6))
ax.set_xticklabels(['12/10','12/11','12/12'])
ax.xaxis.set_minor_locator(AutoMinorLocator(4))
ylabel = 'Significant wave\nheight ($\mathregular{cm}$)'
ax.set_ylabel(ylabel, fontsize=11, fontname='Times New Roman', fontweight='bold')
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(11) for label in labels]
[label.set_color('black') for label in labels]
[label.set_fontweight('bold') for label in labels]
ax.tick_params(which='major', direction='in', colors='xkcd:black', length=6, pad=4)
ax.tick_params(which='minor', direction='in', colors='xkcd:black')

axInv = ax.twinx()
axInv.plot(tt1_U10, U10_var1, color='blue', linestyle='--', linewidth=1, alpha=0.9)
axInv.set_xlim(0, nt1_model-1)
axInv.set_ylim(0, 15)
axInv.xaxis.set_ticks([])
axInv.yaxis.set_ticks(np.linspace(0,15,6))
axInv.tick_params(which='major', direction='in', colors='xkcd:black', length=6, pad=4)
axInv.tick_params(which='minor', direction='in', colors='xkcd:black')
axInv.tick_params(axis='y', colors='blue')
axInv.spines['right'].set_color('blue')
axInv.set_ylabel('Wind speed ($\mathregular{m}$ $\mathregular{{s}^{-1}}$)',
                 color='blue', fontsize=11, fontname='Times New Roman', 
                 fontweight='bold')
labels = axInv.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(11) for label in labels]
[label.set_color('blue') for label in labels]
[label.set_fontweight('bold') for label in labels]

rmse_Hwav_1BF = 0.0
count = 0
for ii, tt in enumerate(tt1_obs2):
    if tt in tt1_model and np.isfinite(Hwav_obs_1BF[ii]):
        indx = np.where(tt1_model==tt)
        rmse_Hwav_1BF = rmse_Hwav_1BF + (Hwav1_1BF[indx] - Hwav_obs_1BF[ii])**2
        count = count + 1
rmse_Hwav_1BF = np.sqrt( rmse_Hwav_1BF / count )
nrmse_Hwav_1BF = rmse_Hwav_1BF / np.nanmean(Hwav1_1BF)
print('Max Hwav: ', np.max(Hwav1_1BF), np.max(Hwav1_2BF))

ax = fig.add_subplot(gs[2,0])
ax.plot(tt1_model, tau1_1BF, color='black', linestyle='-', linewidth=2, alpha=0.9)
ax.set_xlim(0, nt1_model-1)
ax.set_ylim(0, 0.4)
ax.xaxis.set_ticks(np.arange(0,nt1_model,24))
ax.yaxis.set_ticks(np.linspace(0,0.4,5))
ax.set_xticklabels(['12/10','12/11','12/12'])
ax.xaxis.set_minor_locator(AutoMinorLocator(4))
ylabel = 'Bottom shear\nstress ($\mathregular{Pa}$)'
ax.set_ylabel(ylabel, fontsize=11, fontname='Times New Roman', fontweight='bold')
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(11) for label in labels]
[label.set_color('black') for label in labels]
[label.set_fontweight('bold') for label in labels]
ax.tick_params(which='major', direction='in', colors='xkcd:black', length=6, pad=4)
ax.tick_params(which='minor', direction='in', colors='xkcd:black')

print('Max tau: ', np.max(tau1_1BF), np.max(tau1_2BF))

ax = fig.add_subplot(gs[3,0])
ax.plot(tt2_model, h2_1BF, color='black', linestyle='-', linewidth=2, alpha=0.9)
ax.plot(tt2_h, h2b_1BF, color='black', linestyle='--', linewidth=1, alpha=0.9)
ax.set_xlim(0, nt2_model-1)
ax.set_ylim(60, 180)
ax.xaxis.set_ticks(np.arange(0,nt2_model,24))
ax.yaxis.set_ticks(np.linspace(60,180,5))
ax.set_xticklabels(['4/2','4/3','4/4','4/5'])
ax.xaxis.set_minor_locator(AutoMinorLocator(4))
ylabel = 'Water depth\n($\mathregular{cm}$)'
ax.set_ylabel(ylabel, fontsize=11, fontname='Times New Roman', fontweight='bold')
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(11) for label in labels]
[label.set_color('black') for label in labels]
[label.set_fontweight('bold') for label in labels]
ax.tick_params(which='major', direction='in', colors='xkcd:black', length=6, pad=4)
ax.tick_params(which='minor', direction='in', colors='xkcd:black')

ax = fig.add_subplot(gs[4,0])
ax.plot(tt2_model, Hwav2_1BF, color='black', linestyle='-', linewidth=2, alpha=0.9)
ax.set_xlim(0, nt2_model-1)
ax.set_ylim(0, 60)
ax.xaxis.set_ticks(np.arange(0,nt2_model,24))
ax.yaxis.set_ticks(np.linspace(0,60,5))
ax.set_xticklabels(['4/2','4/3','4/4','4/5'])
ax.xaxis.set_minor_locator(AutoMinorLocator(4))
ylabel = 'Significant wave\nheight ($\mathregular{cm}$)'
ax.set_ylabel(ylabel, fontsize=11, fontname='Times New Roman', fontweight='bold')
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(11) for label in labels]
[label.set_color('black') for label in labels]
ax.tick_params(which='major', direction='in', colors='xkcd:black', length=6, pad=4)
ax.tick_params(which='minor', direction='in', colors='xkcd:black')

axInv = ax.twinx()
axInv.plot(tt2_U10, U10_var2, color='blue', linestyle='--', linewidth=1, alpha=0.9)
axInv.set_xlim(0, nt1_model-1)
axInv.set_ylim(0, 18)
axInv.xaxis.set_ticks([])
axInv.yaxis.set_ticks(np.linspace(0,18,4))
axInv.tick_params(which='major', direction='in', colors='xkcd:black', length=6, pad=4)
axInv.tick_params(which='minor', direction='in', colors='xkcd:black')
axInv.tick_params(axis='y', colors='blue')
axInv.spines['right'].set_color('blue')
axInv.set_ylabel('Wind speed ($\mathregular{m}$ $\mathregular{{s}^{-1}}$)',
                 color='blue', fontsize=11, fontname='Times New Roman', 
                 fontweight='bold')
labels = axInv.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(11) for label in labels]
[label.set_color('blue') for label in labels]
[label.set_fontweight('bold') for label in labels]

print('Max Hwav: ', np.max(Hwav2_1BF), np.max(Hwav2_2BF))

ax = fig.add_subplot(gs[5,0])
ax.plot(tt2_model, tau2_1BF, color='black', linestyle='-', linewidth=2, alpha=0.9)
ax.plot(tt2_obs, tau_3d_1BF, color='C3', linestyle='-', marker='.', markersize=5)
ax.set_xlim(0, nt2_model-1)
ax.set_ylim(0, 1)
ax.xaxis.set_ticks(np.arange(0,nt2_model,24))
ax.yaxis.set_ticks(np.linspace(0,1,6))
ax.set_xticklabels(['4/2','4/3','4/4','4/5'])
ax.xaxis.set_minor_locator(AutoMinorLocator(4))
ax.set_xlabel('Time', fontsize=11, fontname='Times New Roman', fontweight='bold')
ylabel = 'Bottom shear\nstress ($\mathregular{Pa}$)'
ax.set_ylabel(ylabel, fontsize=11, fontname='Times New Roman', fontweight='bold')
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(11) for label in labels]
[label.set_color('black') for label in labels]
[label.set_fontweight('bold') for label in labels]
ax.tick_params(which='major', direction='in', colors='xkcd:black', length=6, pad=4)
ax.tick_params(which='minor', direction='in', colors='xkcd:black')

rmse_tau_1BF = 0.0
count = 0
for ii, tt in enumerate(tt2_obs):
    if tt in tt2_model and np.isfinite(tau_3d_1BF[ii]):
        indx = np.where(tt2_model==tt)
        rmse_tau_1BF = rmse_tau_1BF + (tau2_1BF[indx] - tau_3d_1BF[ii])**2
        count = count + 1
rmse_tau_1BF = np.sqrt( rmse_tau_1BF / count )
nrmse_tau_1BF = rmse_tau_1BF / np.nanmean(tau_3d_1BF)
print('Max Hwav: ', np.max(Hwav1_1BF), np.max(Hwav1_2BF))
print('Max tau: ', np.max(tau2_1BF), np.max(tau2_2BF))

# 2BF
ax = fig.add_subplot(gs[0,1])
ax.plot(tt1_model, h1_2BF, color='black', linestyle='-', linewidth=2, alpha=0.9)
ax.plot(tt1_h, h1b_2BF, color='black', linestyle='--', linewidth=1, alpha=0.9)
ax.plot(tt1_obs, h_obs_2BF, color='C3', linestyle='-', marker='.', markersize=5)
ax.set_xlim(0, nt1_model-1)
ax.set_ylim(175, 260)
ax.xaxis.set_ticks(np.arange(0,nt1_model+1,24))
ax.yaxis.set_ticks(np.linspace(180,260,5))
ax.set_xticklabels(['12/10','12/11','12/12'])
ax.xaxis.set_minor_locator(AutoMinorLocator(4))
ax.set_title('2BF', fontsize=14, fontname='Times New Roman', fontweight='bold')
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(11) for label in labels]
[label.set_color('black') for label in labels]
[label.set_fontweight('bold') for label in labels]
ax.tick_params(which='major', direction='in', colors='xkcd:black', length=6, pad=4)
ax.tick_params(which='minor', direction='in', colors='xkcd:black')

rmse_h_2BF = 0.0
count = 0
for ii, tt in enumerate(tt1_obs):
    if tt in tt1_model and np.isfinite(h_obs_2BF[ii]):
        indx = np.where(tt1_model==tt)
        rmse_h_2BF = rmse_h_2BF + (h1_2BF[indx] - h_obs_2BF[ii])**2
        count = count + 1
rmse_h_2BF = np.sqrt( rmse_h_2BF / count )
nrmse_h_2BF = rmse_h_2BF / np.nanmean(h_obs_2BF)

print('RMSE_h_1BF', rmse_h_1BF, nrmse_h_1BF)
print('RMSE_h_2BF', rmse_h_2BF, nrmse_h_2BF)

ax = fig.add_subplot(gs[1,1])
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
[label.set_fontweight('bold') for label in labels]
ax.tick_params(which='major', direction='in', colors='xkcd:black', length=6, pad=4)
ax.tick_params(which='minor', direction='in', colors='xkcd:black')

axInv = ax.twinx()
axInv.plot(tt1_U10, U10_var1, color='blue', linestyle='--', linewidth=1, alpha=0.9)
axInv.set_xlim(0, nt1_model-1)
axInv.set_ylim(0, 15)
axInv.xaxis.set_ticks([])
axInv.yaxis.set_ticks(np.linspace(0,15,6))
axInv.tick_params(which='major', direction='in', colors='xkcd:black', length=6, pad=4)
axInv.tick_params(which='minor', direction='in', colors='xkcd:black')
axInv.tick_params(axis='y', colors='blue')
axInv.spines['right'].set_color('blue')
axInv.set_ylabel('Wind speed ($\mathregular{m}$ $\mathregular{{s}^{-1}}$)',
                 color='blue', fontsize=11, fontname='Times New Roman', 
                 fontweight='bold')
labels = axInv.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(11) for label in labels]
[label.set_color('blue') for label in labels]
[label.set_fontweight('bold') for label in labels]

rmse_Hwav_2BF = 0.0
count = 0
for ii, tt in enumerate(tt1_obs):
    if tt in tt1_model and np.isfinite(Hwav_obs_2BF[ii]):
        indx = np.where(tt1_model==tt)
        rmse_Hwav_2BF = rmse_Hwav_2BF + (Hwav1_2BF[indx] - Hwav_obs_2BF[ii])**2
        count = count + 1
rmse_Hwav_2BF = np.sqrt( rmse_Hwav_2BF / count )
nrmse_Hwav_2BF = rmse_Hwav_2BF / np.nanmean(Hwav1_2BF)

print('RMSE_Hwav_1BF', rmse_Hwav_1BF, nrmse_Hwav_1BF)
print('RMSE_Hwav_2BF', rmse_Hwav_2BF, nrmse_Hwav_2BF)

ax = fig.add_subplot(gs[2,1])
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
[label.set_fontweight('bold') for label in labels]
ax.tick_params(which='major', direction='in', colors='xkcd:black', length=6, pad=4)
ax.tick_params(which='minor', direction='in', colors='xkcd:black')

ax = fig.add_subplot(gs[3,1])
ax.plot(tt2_model, h2_2BF, color='black', linestyle='-', linewidth=2, alpha=0.9)
ax.plot(tt2_h, h2b_2BF, color='black', linestyle='--', linewidth=1, alpha=0.9)
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
[label.set_fontweight('bold') for label in labels]
ax.tick_params(which='major', direction='in', colors='xkcd:black', length=6, pad=4)
ax.tick_params(which='minor', direction='in', colors='xkcd:black')

ax = fig.add_subplot(gs[4,1])
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
[label.set_fontweight('bold') for label in labels]
ax.tick_params(which='major', direction='in', colors='xkcd:black', length=6, pad=4)
ax.tick_params(which='minor', direction='in', colors='xkcd:black')

axInv = ax.twinx()
axInv.plot(tt2_U10, U10_var2, color='blue', linestyle='--', linewidth=1, alpha=0.9)
axInv.set_xlim(0, nt1_model-1)
axInv.set_ylim(0, 18)
axInv.xaxis.set_ticks([])
axInv.yaxis.set_ticks(np.linspace(0,18,4))
axInv.tick_params(which='major', direction='in', colors='xkcd:black', length=6, pad=4)
axInv.tick_params(which='minor', direction='in', colors='xkcd:black')
axInv.tick_params(axis='y', colors='blue')
axInv.spines['right'].set_color('blue')
axInv.set_ylabel('Wind speed ($\mathregular{m}$ $\mathregular{{s}^{-1}}$)',
                 color='blue', fontsize=11, fontname='Times New Roman', 
                 fontweight='bold')
labels = axInv.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(11) for label in labels]
[label.set_color('blue') for label in labels]
[label.set_fontweight('bold') for label in labels]

ax = fig.add_subplot(gs[5,1])
ax.plot(tt2_model, tau2_2BF, color='black', linestyle='-', linewidth=2, alpha=0.9)
ax.plot(tt2_obs, tau_3d_2BF, color='C3', linestyle='-', marker='.', markersize=5)
ax.set_xlim(0, nt2_model-1)
ax.set_ylim(0, 1)
ax.xaxis.set_ticks(np.arange(0,nt2_model,24))
ax.yaxis.set_ticks(np.linspace(0,1,6))
ax.set_xticklabels(['4/2','4/3','4/4','4/5'])
ax.xaxis.set_minor_locator(AutoMinorLocator(4))
ax.set_xlabel('Time', fontsize=11, fontname='Times New Roman', fontweight='bold')
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(11) for label in labels]
[label.set_color('black') for label in labels]
[label.set_fontweight('bold') for label in labels]
ax.tick_params(which='major', direction='in', colors='xkcd:black', length=6, pad=4)
ax.tick_params(which='minor', direction='in', colors='xkcd:black')

rmse_tau_2BF = 0.0
count = 0
for ii, tt in enumerate(tt2_obs):
    if tt in tt2_model and np.isfinite(tau_3d_2BF[ii]):
        indx = np.where(tt2_model==tt)
        rmse_tau_2BF = rmse_tau_2BF + (tau2_2BF[indx] - tau_3d_2BF[ii])**2
        count = count + 1
rmse_tau_2BF = np.sqrt( rmse_tau_2BF / count )
nrmse_tau_2BF = rmse_tau_2BF / np.nanmean(tau_3d_2BF)

print('RMSE_tau_1BF', rmse_tau_1BF, nrmse_tau_1BF)
print('RMSE_tau_2BF', rmse_tau_2BF, nrmse_tau_2BF)

plt.tight_layout()
fig.savefig('F4.png', dpi=300)
fig.savefig('F4.jpg', dpi=600)
plt.show()