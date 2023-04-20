#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  9 22:40:28 2020

Validation of MACES simulated hydrodynamics at Plum Island Estuary

@author: Zeli Tan
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import AutoMinorLocator
from netCDF4 import Dataset
from datetime import date

z_nelson_c = -1.57 + 0.840   # channel (Nelson Island)
z_nelson_m = 0.41 + 0.840    # marsh edge (Nelson Island)
z_shad_c = -0.969 + 0.84     # channel (Shad Creek)
z_shad_m = 0.026 + 0.84      # marsh edge (Shad Creek)

rdir = '/Users/tanz151/Documents/Projects/TAI_BGC/Drafts/Outputs/PlumIsland/'

# read forcing data
t1 = {}
t2 = {}
t1['d0'] = (date(2017,7,19) - date(2012,1,1)).days
t1['d1'] = (date(2017,7,23) - date(2012,1,1)).days
t2['d0'] = (date(2017,10,7) - date(2012,1,1)).days
t2['d1'] = (date(2017,10,11) - date(2012,1,1)).days

filename = rdir + 'force_h.nc'
try:
    nc = Dataset(filename,'r')
    h_var1 = 100 * np.array(nc.variables['h'][t1['d0']*96:t1['d1']*96,0]) # cm
    h_var2 = 100 * np.array(nc.variables['h'][t2['d0']*96:t2['d1']*96,0]) # cm
finally:
    nc.close()
nt1_h = np.size(h_var1)
tt1_h = np.arange(nt1_h) / 4.0
nt2_h = np.size(h_var2)
tt2_h = np.arange(nt2_h) / 4.0

h1b_nelson_c = h_var1 - z_nelson_c*100
h1b_nelson_m = h_var1 - z_nelson_m*100
h1b_shad_c = h_var1 - z_shad_c*100
h1b_shad_m = h_var1 - z_shad_m*100

h2b_nelson_c = h_var2 - z_nelson_c*100
h2b_nelson_m = h_var2 - z_nelson_m*100
h2b_shad_c = h_var2 - z_shad_c*100
h2b_shad_m = h_var2 - z_shad_m*100

# read simulated water level
day0 = (date(2017,7,19) - date(2017,7,17)).days
day1 = (date(2017,7,23) - date(2017,7,17)).days

filename = rdir + 'maces_ecogeom_2017-07-17_2017-08-01_4097.nc'
try:
    nc = Dataset(filename, 'r')
    x = np.array(nc.variables['x'][:])
    zh = np.array(nc.variables['zh'][0,:])
finally:
    nc.close()
    
# find the site index
index1 = np.argmin(np.abs(zh - z_nelson_c))
index2 = np.argmin(np.abs(zh - z_nelson_m))
index3 = np.argmin(np.abs(zh - z_shad_c))
index4 = np.argmin(np.abs(zh - z_shad_m))

filename = rdir + 'maces_hydro_2017-07-17_2017-08-01_4097.nc'
try:
    nc = Dataset(filename, 'r')
    h1_nelson_c = 100*np.array(nc.variables['h'][day0*24-1:day1*24-1,index1])
    h1_nelson_m = 100*np.array(nc.variables['h'][day0*24-1:day1*24-1,index2])
    h1_shad_c = 100*np.array(nc.variables['h'][day0*24-1:day1*24-1,index3])
    h1_shad_m = 100*np.array(nc.variables['h'][day0*24-1:day1*24-1,index4])
finally:
    nc.close()
nt1_model = np.size(h1_nelson_c)
tt1_model = np.arange(nt1_model)

day0 = (date(2017,10,7) - date(2017,10,4)).days
day1 = (date(2017,10,11) - date(2017,10,4)).days
filename = rdir + 'maces_hydro_2017-10-04_2017-10-21_4097.nc'
try:
    nc = Dataset(filename, 'r')
    h2_nelson_c = 100*np.array(nc.variables['h'][day0*24-1:day1*24-1,index1])
    h2_nelson_m = 100*np.array(nc.variables['h'][day0*24-1:day1*24-1,index2])
    h2_shad_c = 100*np.array(nc.variables['h'][day0*24-1:day1*24-1,index3])
    h2_shad_m = 100*np.array(nc.variables['h'][day0*24-1:day1*24-1,index4])
finally:
    nc.close()
nt2_model = np.size(h2_nelson_c)
tt2_model = np.arange(nt2_model)

# read measured water level forcings
filename = rdir + 'MAR-RO-Wtable-Nel-2017.xls'
df = pd.read_excel(filename, sheet_name='MAR-RO-Wtable-Nel-2017', header=0, 
                   usecols='A,B,C,K')
df.columns = ['Date', 'Time', 'Gauge', 'N204']
h1_nelson_gauge = 100 * np.array(df['Gauge'])[16895:18047]
h1_nelson_n204 = 100 * np.array(df['N204'])[16895:18047] - 198
h2_nelson_gauge = 100 * np.array(df['Gauge'])[39707:40859] 
h2_nelson_n204 =  100 * np.array(df['N204'])[39707:40859] - 198
h1_nelson_gauge[h1_nelson_gauge<0] = 0
h1_nelson_n204[h1_nelson_n204<0] = 0
h2_nelson_gauge[h2_nelson_gauge<0] = 0
h2_nelson_n204[h2_nelson_n204<0] = 0

nt1_obs = np.size(h1_nelson_gauge)
tt1_obs = np.arange(nt1_obs) / 12

nt2_obs = np.size(h2_nelson_gauge)
tt2_obs = np.arange(nt2_obs) / 12

filename = rdir + 'MAR-RO-Wtable-Shad-2017.xls'
df = pd.read_excel(filename, sheet_name='MAR-RO-Wtable-Shad-2017', header=0, 
                   usecols='A,B,C,G')
df.columns = ['Date', 'Time', 'Gauge', 'S102']
h1_shad_gauge = 100 * np.array(df['Gauge'])[17086:18238]
h1_shad_s102 = 100 * np.array(df['S102'])[17086:18238] - 98.5
h2_shad_gauge = 100 * np.array(df['Gauge'])[39967:41119]
h2_shad_s102 = 100 * np.array(df['S102'])[39967:41119] - 98.5
h1_shad_gauge[h1_shad_gauge<0] = 0
h1_shad_s102[h1_shad_s102<0] = 0
h2_shad_gauge[h2_shad_gauge<0] = 0
h2_shad_s102[h2_shad_s102<0] = 0
    
# plot
plt.clf()
fig = plt.figure(figsize=(8,10))

gs = gridspec.GridSpec(nrows=2, ncols=2)

plt.style.use('default')

# nelson island channel
ax = fig.add_subplot(gs[0,0])
ax.plot(tt1_model, h1_nelson_c, color='black', linestyle='-', linewidth=2, alpha=0.9)
ax.plot(tt1_h, h1b_nelson_c, color='black', linestyle='--', linewidth=1, alpha=0.9)
ax.plot(tt1_obs, h1_nelson_gauge, color='C3', linestyle='--', linewidth=2)
ax.set_xlim(0, nt1_model)
ax.set_ylim(0, 300)
ax.xaxis.set_ticks(np.arange(0,nt1_model+1,24))
ax.yaxis.set_ticks(np.linspace(0,300,6))
ax.set_xticklabels(['7/19','7/20','7/21','7/22','7/23'])
ax.xaxis.set_minor_locator(AutoMinorLocator(4))
#ax.set_xlabel('Time', fontsize=11, fontname='Times New Roman', color='black')
ylabel = 'Water depth ($\mathregular{cm}$)'
ax.set_ylabel(ylabel, fontsize=12, fontname='Times New Roman', color='black')
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(12) for label in labels]
[label.set_color('black') for label in labels]
[label.set_fontweight('bold') for label in labels]
ax.text(0.05, 0.90, 'a', transform=ax.transAxes, fontsize=16,
        fontname='Times New Roman', fontweight='bold')
ax.tick_params(which='major', direction='in', colors='xkcd:black', length=6, pad=8)
ax.tick_params(which='minor', direction='in', colors='xkcd:black')

rmse_h1_nelson_c = 0.0
count = 0
for ii, tt in enumerate(tt1_obs):
    if tt in tt1_model and np.isfinite(h1_nelson_gauge[ii]):
        indx = np.where(tt1_model==tt)
        rmse_h1_nelson_c = rmse_h1_nelson_c + (h1_nelson_c[indx] - h1_nelson_gauge[ii])**2
        count = count + 1
rmse_h1_nelson_c = np.sqrt( rmse_h1_nelson_c / count )

ax = fig.add_subplot(gs[0,1])
ax.plot(tt2_model, h2_nelson_c, color='black', linestyle='-', linewidth=2, alpha=0.9)
ax.plot(tt2_h, h2b_nelson_c, color='black', linestyle='--', linewidth=1, alpha=0.9)
ax.plot(tt2_obs, h2_nelson_gauge, color='C3', linestyle='--', linewidth=2)
ax.set_xlim(0, nt2_model)
ax.set_ylim(0, 300)
ax.xaxis.set_ticks(np.arange(0,nt2_model+1,24))
ax.yaxis.set_ticks(np.linspace(0,300,6))
ax.set_xticklabels(['10/7','10/8','10/9','10/10','10/11'])
ax.xaxis.set_minor_locator(AutoMinorLocator(4))
#ax.set_xlabel('Time', fontsize=11, fontname='Times New Roman', color='black')
#ylabel = 'Water depth ($\mathregular{cm}$)'
#ax.set_ylabel(ylabel, fontsize=11, fontname='Times New Roman', color='black')
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(12) for label in labels]
[label.set_color('black') for label in labels]
[label.set_fontweight('bold') for label in labels]
ax.text(0.05, 0.90, 'b', transform=ax.transAxes, fontsize=16,
        fontname='Times New Roman', fontweight='bold')
ax.tick_params(which='major', direction='in', colors='xkcd:black', length=6, pad=8)
ax.tick_params(which='minor', direction='in', colors='xkcd:black')

rmse_h2_nelson_c = 0.0
count = 0
for ii, tt in enumerate(tt2_obs):
    if tt in tt2_model and np.isfinite(h2_nelson_gauge[ii]):
        indx = np.where(tt2_model==tt)
        rmse_h2_nelson_c = rmse_h2_nelson_c + (h2_nelson_c[indx] - h2_nelson_gauge[ii])**2
        count = count + 1
rmse_h2_nelson_c = np.sqrt( rmse_h2_nelson_c / count )

print('RMSE_h1_nelson_c', rmse_h1_nelson_c)
print('RMSE_h2_nelson_c', rmse_h2_nelson_c)

# nelson island marsh edge
axax = fig.add_subplot(gs[1,0])
h1_nelson_cs = h1_nelson_c - 198
h1_nelson_cs[h1_nelson_cs<0] = 0.0
ax.plot(tt1_model, h1_nelson_m, color='black', linestyle='-', linewidth=2, alpha=0.9)
ax.plot(tt1_h, h1b_nelson_m, color='black', linestyle='--', linewidth=1, alpha=0.9)
ax.plot(tt1_obs, h1_nelson_n204, color='C3', linestyle='--', linewidth=2)
#ax.plot(tt1_model, h1_nelson_cs, color='C0', linestyle='-', linewidth=2, alpha=0.9)
ax.set_xlim(0, nt1_model)
ax.set_ylim(0, 100)
ax.xaxis.set_ticks(np.arange(0,nt1_model+1,24))
ax.yaxis.set_ticks(np.linspace(0,100,6))
ax.set_xticklabels(['7/19','7/20','7/21','7/22','7/23'])
ax.xaxis.set_minor_locator(AutoMinorLocator(4))
#ax.set_xlabel('Time', fontsize=11, fontname='Times New Roman', color='black')
ylabel = 'Water depth ($\mathregular{cm}$)'
ax.set_ylabel(ylabel, fontsize=12, fontname='Times New Roman', color='black')
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(12) for label in labels]
[label.set_color('black') for label in labels]
[label.set_fontweight('bold') for label in labels]
ax.text(0.05, 0.90, 'c', transform=ax.transAxes, fontsize=16,
        fontname='Times New Roman', fontweight='bold')
ax.tick_params(which='major', direction='in', colors='xkcd:black', length=6, pad=8)
ax.tick_params(which='minor', direction='in', colors='xkcd:black')

rmse_h1_nelson_m = 0.0
count = 0
for ii, tt in enumerate(tt1_obs):
    if tt in tt1_model and np.isfinite(h1_nelson_n204[ii]):
        indx = np.where(tt1_model==tt)
        rmse_h1_nelson_m = rmse_h1_nelson_m + (h1_nelson_m[indx] - h1_nelson_n204[ii])**2
        count = count + 1
rmse_h1_nelson_m = np.sqrt( rmse_h1_nelson_m / count )

ax = fig.add_subplot(gs[1,1])
ax.plot(tt2_model, h2_nelson_m, color='black', linestyle='-', linewidth=2, alpha=0.9)
ax.plot(tt2_h, h2b_nelson_m, color='black', linestyle='--', linewidth=1, alpha=0.9)
ax.plot(tt2_obs, h2_nelson_n204, color='C3', linestyle='--', linewidth=2)
ax.set_xlim(0, nt2_model)
ax.set_ylim(0, 100)
ax.xaxis.set_ticks(np.arange(0,nt2_model+1,24))
ax.yaxis.set_ticks(np.linspace(0,100,6))
ax.set_xticklabels(['10/7','10/8','10/9','10/10','10/11'])
ax.xaxis.set_minor_locator(AutoMinorLocator(4))
#ax.set_xlabel('Time', fontsize=11, fontname='Times New Roman', color='black')
#ylabel = 'Water depth ($\mathregular{cm}$)'
#ax.set_ylabel(ylabel, fontsize=11, fontname='Times New Roman', color='black')
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(12) for label in labels]
[label.set_color('black') for label in labels]
[label.set_fontweight('bold') for label in labels]
ax.text(0.05, 0.90, 'd', transform=ax.transAxes, fontsize=16,
        fontname='Times New Roman', fontweight='bold')
ax.tick_params(which='major', direction='in', colors='xkcd:black', length=6, pad=8)
ax.tick_params(which='minor', direction='in', colors='xkcd:black')

rmse_h2_nelson_m = 0.0
count = 0
for ii, tt in enumerate(tt2_obs):
    if tt in tt2_model and np.isfinite(h2_nelson_n204[ii]):
        indx = np.where(tt2_model==tt)
        rmse_h2_nelson_m = rmse_h2_nelson_m + (h2_nelson_m[indx] - h2_nelson_n204[ii])**2
        count = count + 1
rmse_h2_nelson_m = np.sqrt( rmse_h2_nelson_m / count )

print('RMSE_h1_nelson_m', rmse_h1_nelson_m)
print('RMSE_h2_nelson_m', rmse_h2_nelson_m)

## shad creek channel
#ax = axes[2][0]
#ax.plot(tt1_model, h1_shad_c, color='black', linestyle='-', linewidth=2, alpha=0.9)
#ax.plot(tt1_obs, h1_shad_gauge, color='C3', linestyle='--', linewidth=2)
#ax.set_xlim(0, nt1_model)
#ax.set_ylim(0, 300)
#ax.xaxis.set_ticks(np.arange(0,nt1_model+1,24))
#ax.yaxis.set_ticks(np.linspace(0,300,6))
#ax.set_xticklabels(['7/19','7/20','7/21','7/22','7/23'])
#ax.xaxis.set_minor_locator(AutoMinorLocator(4))
##ax.set_xlabel('Time', fontsize=11, fontname='Times New Roman', color='black')
#ylabel = 'Water depth ($\mathregular{cm}$)'
#ax.set_ylabel(ylabel, fontsize=12, fontname='Times New Roman', color='black')
#labels = ax.get_xticklabels() + ax.get_yticklabels()
#[label.set_fontname('Times New Roman') for label in labels]
#[label.set_fontsize(11) for label in labels]
#[label.set_color('black') for label in labels]
#ax.text(0.05, 0.90, 'e', transform=ax.transAxes, fontsize=16,
#        fontname='Times New Roman', fontweight='bold')
#ax.tick_params(which='major', direction='in', colors='xkcd:black', length=6, pad=8)
#ax.tick_params(which='minor', direction='in', colors='xkcd:black')
#
#ax = axes[2][1]
#ax.plot(tt2_model, h2_shad_c, color='black', linestyle='-', linewidth=2, alpha=0.9)
#ax.plot(tt2_obs, h2_shad_gauge, color='C3', linestyle='--', linewidth=2)
#ax.set_xlim(0, nt2_model)
#ax.set_ylim(0, 300)
#ax.xaxis.set_ticks(np.arange(0,nt2_model+1,24))
#ax.yaxis.set_ticks(np.linspace(0,300,6))
#ax.set_xticklabels(['10/7','10/8','10/9','10/10','10/11'])
#ax.xaxis.set_minor_locator(AutoMinorLocator(4))
##ax.set_xlabel('Time', fontsize=11, fontname='Times New Roman', color='black')
##ylabel = 'Water depth ($\mathregular{cm}$)'
##ax.set_ylabel(ylabel, fontsize=11, fontname='Times New Roman', color='black')
#labels = ax.get_xticklabels() + ax.get_yticklabels()
#[label.set_fontname('Times New Roman') for label in labels]
#[label.set_fontsize(11) for label in labels]
#[label.set_color('black') for label in labels]
#ax.text(0.05, 0.90, 'f', transform=ax.transAxes, fontsize=16,
#        fontname='Times New Roman', fontweight='bold')
#ax.tick_params(which='major', direction='in', colors='xkcd:black', length=6, pad=8)
#ax.tick_params(which='minor', direction='in', colors='xkcd:black')
#
## shad creek marsh edge
#ax = axes[3][0]
#ax.plot(tt1_model, h1_shad_m, color='black', linestyle='-', linewidth=2, alpha=0.9)
#ax.plot(tt1_obs, h1_shad_s102, color='C3', linestyle='--', linewidth=2)
#ax.set_xlim(0, nt1_model)
#ax.set_ylim(0, 100)
#ax.xaxis.set_ticks(np.arange(0,nt1_model+1,24))
#ax.yaxis.set_ticks(np.linspace(0,100,6))
#ax.set_xticklabels(['7/19','7/20','7/21','7/22','7/23'])
#ax.xaxis.set_minor_locator(AutoMinorLocator(4))
#ax.set_xlabel('Time', fontsize=12, fontname='Times New Roman', color='black')
#ylabel = 'Water depth ($\mathregular{cm}$)'
#ax.set_ylabel(ylabel, fontsize=12, fontname='Times New Roman', color='black')
#labels = ax.get_xticklabels() + ax.get_yticklabels()
#[label.set_fontname('Times New Roman') for label in labels]
#[label.set_fontsize(11) for label in labels]
#[label.set_color('black') for label in labels]
#ax.text(0.05, 0.90, 'g', transform=ax.transAxes, fontsize=16,
#        fontname='Times New Roman', fontweight='bold')
#ax.tick_params(which='major', direction='in', colors='xkcd:black', length=6, pad=8)
#ax.tick_params(which='minor', direction='in', colors='xkcd:black')
#
#ax = axes[3][1]
#ax.plot(tt2_model, h2_shad_m, color='black', linestyle='-', linewidth=2, alpha=0.9)
#ax.plot(tt2_obs, h2_shad_s102, color='C3', linestyle='--', linewidth=2)
#ax.set_xlim(0, nt2_model)
#ax.set_ylim(0, 100)
#ax.xaxis.set_ticks(np.arange(0,nt2_model+1,24))
#ax.yaxis.set_ticks(np.linspace(0,100,6))
#ax.set_xticklabels(['10/7','10/8','10/9','10/10','10/11'])
#ax.xaxis.set_minor_locator(AutoMinorLocator(4))
#ax.set_xlabel('Time', fontsize=12, fontname='Times New Roman', color='black')
##ylabel = 'Water depth ($\mathregular{cm}$)'
##ax.set_ylabel(ylabel, fontsize=11, fontname='Times New Roman', color='black')
#labels = ax.get_xticklabels() + ax.get_yticklabels()
#[label.set_fontname('Times New Roman') for label in labels]
#[label.set_fontsize(11) for label in labels]
#[label.set_color('black') for label in labels]
#ax.text(0.05, 0.90, 'h', transform=ax.transAxes, fontsize=16,
#        fontname='Times New Roman', fontweight='bold')
#ax.tick_params(which='major', direction='in', colors='xkcd:black', length=6, pad=8)
#ax.tick_params(which='minor', direction='in', colors='xkcd:black')

plt.tight_layout()
fig.savefig('F5.png', dpi=300)
fig.savefig('F5.jpg', dpi=600)
plt.show()