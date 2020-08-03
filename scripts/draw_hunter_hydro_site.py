#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  9 22:40:28 2020

Validation of MACES simulated hydrodynamics at Hunter Estuary

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

day0 = (date(2004,9,28) - date(2004,9,25)).days
day1 = (date(2004,10,1) - date(2004,9,25)).days
# read simulated water depth and suspended sediment
rdir = '/Users/tanz151/Documents/Projects/TAI_BGC/Data/Hydrodynamics_obs/HunterEstuary/Outputs/'
filename = rdir + 'maces_ecogeom_2004-09-25_2004-10-06_11099.nc'
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
    
filename = rdir + 'maces_hydro_2004-09-25_2004-10-06_11099.nc'
try:
    nc = Dataset(filename, 'r')
    h1_sim = 100*np.array(nc.variables['h'][day0*24:day1*24,index1])
    h2_sim = 100*np.array(nc.variables['h'][day0*24:day1*24,index2])
    h3_sim = 100*np.array(nc.variables['h'][day0*24:day1*24,index3])
    h4_sim = 100*np.array(nc.variables['h'][day0*24:day1*24,index4])
finally:
    nc.close()
nt_sim = np.size(h1_sim)
tt_sim = np.arange(nt_sim)

# read water level forcing
day0 = (date(2004,9,28) - date(2004,1,1)).days
day1 = (date(2004,10,1) - date(2004,1,1)).days
filename = '/Users/tanz151/Documents/Projects/TAI_BGC/Data/Hydrodynamics_obs/' + \
    'HunterEstuary/force_h.nc'
try:
    nc = Dataset(filename,'r')
    h_force = 100*np.array(nc.variables['h'][day0*96:day1*96,0])
finally:
    nc.close()
nt_force = np.size(h_force)
tt_force = np.arange(nt_force)

filename = '/Users/tanz151/Documents/Projects/TAI_BGC/Data/Hydrodynamics_obs/' + \
    'HunterEstuary/force_SSC.nc'
try:
    nc = Dataset(filename,'r')
    css_force = np.array(nc.variables['TSM'][day0*96:day1*96,0])
finally:
    nc.close()

# read benchmark estimates
filename = '/Users/tanz151/Documents/Projects/TAI_BGC/Data/Hydrodynamics_obs/' + \
    'HunterEstuary/water_depth_benchmark_2004-09-28_2004-09-30.nc'
try:
    nc = Dataset(filename,'r')
    h1_obs = 100*np.array(nc.variables['h'][:,26,41]) - 5
    h2_obs = 100*np.array(nc.variables['h'][:,26,36])
    h3_obs = 100*np.array(nc.variables['h'][:,26,33]) - 5
    h4_obs = 100*np.array(nc.variables['h'][:,27,27])
finally:
    nc.close()
h1_obs[h1_obs<0] = 0
h2_obs[h2_obs<0] = 0
h3_obs[h3_obs<0] = 0
h4_obs[h4_obs<0] = 0
nt_obs = np.size(h1_obs)
tt_obs = np.arange(nt_obs) / 4.
    
# plot
plt.clf()
fig, axes = plt.subplots(3, 2, figsize=(8,9))

plt.style.use('default')

# forcing
ax = axes[0][0]
ax.plot(tt_force, h_force, color='black', linestyle='-', linewidth=2, alpha=0.9)
ax.set_xlim(0, nt_force)
ax.set_ylim(-80, 80)
ax.xaxis.set_ticks(np.arange(0,nt_force+1,96))
ax.yaxis.set_ticks(np.linspace(-80,80,9))
ax.set_xticklabels(['9/28','9/29','9/30','10/1'])
ax.xaxis.set_minor_locator(AutoMinorLocator(4))
ylabel = 'Water level ($\mathregular{cm}$ $\mathregular{AHD}$)'
ax.set_ylabel(ylabel, fontsize=12, fontname='Times New Roman', color='black')
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(12) for label in labels]
[label.set_color('black') for label in labels]
ax.text(0.05, 0.90, 'a', transform=ax.transAxes, fontsize=16,
        fontname='Times New Roman', fontweight='bold')
ax.tick_params(which='major', direction='in', colors='xkcd:black', length=6, pad=8)
ax.tick_params(which='minor', direction='in', colors='xkcd:black')

# blank subplot
ax = axes[0][1]
ax.plot(tt_force, css_force, color='black', linestyle='-', linewidth=2, alpha=0.9)
ax.set_xlim(0, nt_force)
ax.set_ylim(10, 35)
ax.xaxis.set_ticks(np.arange(0,nt_force+1,96))
ax.yaxis.set_ticks(np.linspace(10,35,6))
ax.set_xticklabels(['9/28','9/29','9/30','10/1'])
ax.xaxis.set_minor_locator(AutoMinorLocator(4))
ylabel = 'Suspended sediment ($\mathregular{mg}$ $\mathregular{l^{-1}}$)'
ax.set_ylabel(ylabel, fontsize=12, fontname='Times New Roman', color='black')
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(12) for label in labels]
[label.set_color('black') for label in labels]
ax.text(0.05, 0.90, 'b', transform=ax.transAxes, fontsize=16,
        fontname='Times New Roman', fontweight='bold')
ax.tick_params(which='major', direction='in', colors='xkcd:black', length=6, pad=8)
ax.tick_params(which='minor', direction='in', colors='xkcd:black')

# channel
ax = axes[1][0]
ax.plot(tt_sim, h1_sim, color='black', linestyle='-', linewidth=2, alpha=0.9)
ax.plot(tt_obs, h1_obs, color='C3', linestyle='--', linewidth=2)
ax.set_xlim(0, nt_sim)
ax.set_ylim(0, 100)
ax.xaxis.set_ticks(np.arange(0,nt_sim+1,24))
ax.yaxis.set_ticks(np.linspace(0,100,6))
ax.set_xticklabels(['9/28','9/29','9/30','10/1'])
ax.xaxis.set_minor_locator(AutoMinorLocator(4))
#ax.set_xlabel('Time', fontsize=11, fontname='Times New Roman', color='black')
ylabel = 'Water depth ($\mathregular{cm}$)'
ax.set_ylabel(ylabel, fontsize=12, fontname='Times New Roman', color='black')
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(12) for label in labels]
[label.set_color('black') for label in labels]
ax.text(0.05, 0.90, 'c', transform=ax.transAxes, fontsize=16,
        fontname='Times New Roman', fontweight='bold')
ax.tick_params(which='major', direction='in', colors='xkcd:black', length=6, pad=8)
ax.tick_params(which='minor', direction='in', colors='xkcd:black')

rmse_h1 = 0.0
count = 0
for ii, tt in enumerate(tt_obs):
    if tt in tt_sim and np.isfinite(h1_obs[ii]):
        indx = np.where(tt_sim==tt)
        rmse_h1 = rmse_h1 + (h1_sim[indx] - h1_obs[ii])**2
        count = count + 1
rmse_h1 = np.sqrt( rmse_h1 / count )

# mangrove edge
ax = axes[1][1]
ax.plot(tt_sim, h2_sim, color='black', linestyle='-', linewidth=2, alpha=0.9)
ax.plot(tt_obs, h2_obs, color='C3', linestyle='--', linewidth=2)
ax.set_xlim(0, nt_sim)
ax.set_ylim(0, 100)
ax.xaxis.set_ticks(np.arange(0,nt_sim+1,24))
ax.yaxis.set_ticks(np.linspace(0,100,6))
ax.set_xticklabels(['9/28','9/29','9/30','10/1'])
ax.xaxis.set_minor_locator(AutoMinorLocator(4))
#ax.set_xlabel('Time', fontsize=11, fontname='Times New Roman', color='black')
ylabel = 'Water depth ($\mathregular{cm}$)'
ax.set_ylabel(ylabel, fontsize=12, fontname='Times New Roman', color='black')
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(12) for label in labels]
[label.set_color('black') for label in labels]
ax.text(0.05, 0.90, 'd', transform=ax.transAxes, fontsize=16,
        fontname='Times New Roman', fontweight='bold')
ax.tick_params(which='major', direction='in', colors='xkcd:black', length=6, pad=8)
ax.tick_params(which='minor', direction='in', colors='xkcd:black')

rmse_h2 = 0.0
count = 0
for ii, tt in enumerate(tt_obs):
    if tt in tt_sim and np.isfinite(h2_obs[ii]):
        indx = np.where(tt_sim==tt)
        rmse_h2 = rmse_h2 + (h2_sim[indx] - h2_obs[ii])**2
        count = count + 1
rmse_h2 = np.sqrt( rmse_h2 / count )

# mangrove interior
ax = axes[2][0]
ax.plot(tt_sim, h3_sim, color='black', linestyle='-', linewidth=2, alpha=0.9)
ax.plot(tt_obs, h3_obs, color='C3', linestyle='--', linewidth=2)
#ax.plot(tt_force, h_force-38, color='C0', linestyle='-', linewidth=2)
ax.set_xlim(0, nt_sim)
ax.set_ylim(0, 50)
ax.xaxis.set_ticks(np.arange(0,nt_sim+1,24))
ax.yaxis.set_ticks(np.linspace(0,50,6))
ax.set_xticklabels(['9/28','9/29','9/30','10/1'])
ax.xaxis.set_minor_locator(AutoMinorLocator(4))
ax.set_xlabel('Time', fontsize=12, fontname='Times New Roman', color='black')
ylabel = 'Water depth ($\mathregular{cm}$)'
ax.set_ylabel(ylabel, fontsize=12, fontname='Times New Roman', color='black')
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(12) for label in labels]
[label.set_color('black') for label in labels]
ax.text(0.05, 0.90, 'e', transform=ax.transAxes, fontsize=16,
        fontname='Times New Roman', fontweight='bold')
ax.tick_params(which='major', direction='in', colors='xkcd:black', length=6, pad=8)
ax.tick_params(which='minor', direction='in', colors='xkcd:black')

rmse_h3 = 0.0
count = 0
for ii, tt in enumerate(tt_obs):
    if tt in tt_sim and np.isfinite(h3_obs[ii]):
        indx = np.where(tt_sim==tt)
        rmse_h3 = rmse_h3 + (h3_sim[indx] - h3_obs[ii])**2
        count = count + 1
rmse_h3 = np.sqrt( rmse_h3 / count )

# marsh edge
ax = axes[2][1]
ax.plot(tt_sim, h4_sim, color='black', linestyle='-', linewidth=2, alpha=0.9)
ax.plot(tt_obs, h4_obs, color='C3', linestyle='--', linewidth=2)
#ax.plot(tt_force, h_force-65, color='C0', linestyle='-', linewidth=1)
ax.set_xlim(0, nt_sim)
ax.set_ylim(0, 10)
ax.xaxis.set_ticks(np.arange(0,nt_sim+1,24))
ax.yaxis.set_ticks(np.linspace(0,10,6))
ax.set_xticklabels(['9/28','9/29','9/30','10/1'])
ax.xaxis.set_minor_locator(AutoMinorLocator(4))
ax.set_xlabel('Time', fontsize=12, fontname='Times New Roman', color='black')
ylabel = 'Water depth ($\mathregular{cm}$)'
ax.set_ylabel(ylabel, fontsize=12, fontname='Times New Roman', color='black')
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(12) for label in labels]
[label.set_color('black') for label in labels]
ax.text(0.05, 0.90, 'f', transform=ax.transAxes, fontsize=16,
        fontname='Times New Roman', fontweight='bold')
ax.tick_params(which='major', direction='in', colors='xkcd:black', length=6, pad=8)
ax.tick_params(which='minor', direction='in', colors='xkcd:black')

rmse_h4 = 0.0
count = 0
for ii, tt in enumerate(tt_obs):
    if tt in tt_sim and np.isfinite(h4_obs[ii]):
        indx = np.where(tt_sim==tt)
        rmse_h4 = rmse_h4 + (h4_sim[indx] - h4_obs[ii])**2
        count = count + 1
rmse_h4 = np.sqrt( rmse_h4 / count )

print('RMSE h: ', rmse_h1, rmse_h2, rmse_h3, rmse_h4)

plt.tight_layout()
fig.savefig('F11.png', dpi=300)
fig.savefig('F11.pdf', dpi=600)
plt.show()