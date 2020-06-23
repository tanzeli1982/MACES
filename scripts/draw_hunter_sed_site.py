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

day0 = (date(2004,9,28) - date(2004,9,25)).days
day1 = (date(2004,10,1) - date(2004,9,25)).days
# read simulated water depth and suspended sediment
filename = '/Users/tanz151/Python_maces/src/maces_ecogeom_2004-09-25_2004-10-06_11099.nc'
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

filename = '/Users/tanz151/Python_maces/src/maces_hydro_2004-09-25_2004-10-06_11099.nc'
try:
    nc = Dataset(filename, 'r')
    sed1_sim = 1e3*np.array(nc.variables['TSM'][day0*24:day1*24,index1])
    sed2_sim = 1e3*np.array(nc.variables['TSM'][day0*24:day1*24,index2])
    sed3_sim = 1e3*np.array(nc.variables['TSM'][day0*24:day1*24,index3])
    sed4_sim = 1e3*np.array(nc.variables['TSM'][day0*24:day1*24,index4])
finally:
    nc.close()
nt_sim = np.size(sed1_sim)
tt_sim = np.arange(nt_sim)

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

# channel
ax = axes[0][0]
ax.plot(tt_sim, sed1_sim, color='black', linestyle='-', linewidth=2, alpha=0.9)
ax.plot(tt_obs, sed1_obs, color='C3', linestyle='--', linewidth=2)
ax.set_xlim(0, nt_sim)
#ax.set_ylim(0, 100)
ax.xaxis.set_ticks(np.arange(0,nt_sim+1,24))
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
ax.plot(tt_sim, sed2_sim, color='black', linestyle='-', linewidth=2, alpha=0.9)
ax.plot(tt_obs, sed2_obs, color='C3', linestyle='--', linewidth=2)
ax.set_xlim(0, nt_sim)
#ax.set_ylim(0, 100)
ax.xaxis.set_ticks(np.arange(0,nt_sim+1,24))
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
ax.plot(tt_sim, sed3_sim, color='black', linestyle='-', linewidth=2, alpha=0.9)
ax.plot(tt_obs, sed3_obs, color='C3', linestyle='--', linewidth=2)
ax.set_xlim(0, nt_sim)
#ax.set_ylim(0, 50)
ax.xaxis.set_ticks(np.arange(0,nt_sim+1,24))
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
ax.plot(tt_sim, sed4_sim, color='black', linestyle='-', linewidth=2, alpha=0.9)
ax.plot(tt_obs, sed4_obs, color='C3', linestyle='--', linewidth=2)
ax.set_xlim(0, nt_sim)
#ax.set_ylim(0, 10)
ax.xaxis.set_ticks(np.arange(0,nt_sim+1,24))
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
#fig.savefig('F12.pdf', dpi=600)
plt.show()