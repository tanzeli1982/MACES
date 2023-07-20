#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 21 06:44:46 2020

Draw the forcing data of Hunter Estuary

@author: Zeli Tan
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import AutoMinorLocator
from netCDF4 import Dataset
from datetime import date

t1 = {}
t1['d0'] = (date(2004,9,28) - date(2004,1,1)).days
t1['d1'] = (date(2004,10,1) - date(2004,1,1)).days

# read data
rdir = '/Users/tanz151/Documents/Projects/TAI_BGC/Drafts/Outputs/HunterEstuary/'
filename = rdir + 'force_h.nc'
try:
    nc = Dataset(filename,'r')
    h_var1 = 100 * np.array(nc.variables['h'][t1['d0']*96:t1['d1']*96,0]) # cm
finally:
    nc.close()
nt1_h = len(h_var1)
tt1_h = np.arange(nt1_h)
    
filename = rdir + 'force_U10.nc'
try:
    nc = Dataset(filename,'r')
    U10_var1 = np.array(nc.variables['U10'][t1['d0']*96:t1['d1']*96,0])
finally:
    nc.close()
nt1_U10 = len(U10_var1)
tt1_U10 = np.arange(nt1_U10)

filename = rdir + 'force_SSC.nc'
try:
    nc = Dataset(filename,'r')
    SSC_var1 = np.array(nc.variables['TSM'][t1['d0']*96:t1['d1']*96,0])
finally:
    nc.close()
nt1_SSC = len(SSC_var1)
tt1_SSC = np.arange(nt1_SSC)
    
# plot
plt.clf()
fig = plt.figure(figsize=(4,6))

gs = gridspec.GridSpec(nrows=2, ncols=1)

plt.style.use('default')

ax = fig.add_subplot(gs[0,0])
ax.plot(tt1_h, h_var1, color='black', linestyle='-', linewidth=2, alpha=0.9)
ax.set_xlim(0, nt1_h)
ax.set_ylim(-80, 80)
ax.xaxis.set_ticks(np.arange(0,nt1_h+1,96))
ax.yaxis.set_ticks(np.linspace(-80,80,9))
ax.set_xticklabels(['9/28','9/29','9/30','10/1'])
ax.xaxis.set_minor_locator(AutoMinorLocator(4))
#ax.set_xlabel('Time', fontsize=12, fontname='Times New Roman', fontweight='bold')
ylabel = 'Water level ($\mathregular{cm}$ $\mathregular{AHD}$)'
ax.set_ylabel(ylabel, fontsize=12, fontname='Times New Roman', fontweight='bold')
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(12) for label in labels]
[label.set_color('black') for label in labels]
[label.set_fontweight('bold') for label in labels]
ax.tick_params(which='major', direction='in', colors='xkcd:black', length=6, pad=8)
ax.tick_params(which='minor', direction='in', colors='xkcd:black')

#ax = fig.add_subplot(gs[0,1])
#ax.plot(tt1_U10, U10_var1, color='black', linestyle='-', linewidth=2, alpha=0.9)
#ax.set_xlim(0, nt1_U10)
#ax.set_ylim(0, 9)
#ax.xaxis.set_ticks(np.arange(0,nt1_U10+1,96))
##ax.yaxis.set_ticks(np.linspace(3,18,6))
#ax.set_xticklabels(['9/28','9/29','9/30','10/1'])
#ax.xaxis.set_minor_locator(AutoMinorLocator(4))
#ax.set_xlabel('Time', fontsize=12, fontname='Times New Roman', fontweight='bold')
#ylabel = 'Wind speed ($\mathregular{m}$ $\mathregular{{s}^{-1}}$)'
#ax.set_ylabel(ylabel, fontsize=12, fontname='Times New Roman', fontweight='bold')
#labels = ax.get_xticklabels() + ax.get_yticklabels()
#[label.set_fontname('Times New Roman') for label in labels]
#[label.set_fontsize(12) for label in labels]
#[label.set_color('black') for label in labels]
#[label.set_fontweight('bold') for label in labels]
#ax.tick_params(which='major', direction='in', colors='xkcd:black', length=6, pad=8)
#ax.tick_params(which='minor', direction='in', colors='xkcd:black')

ax = fig.add_subplot(gs[1,0])
ax.plot(tt1_SSC, SSC_var1, color='black', linestyle='-', linewidth=2, alpha=0.9)
ax.set_xlim(0, nt1_SSC)
ax.set_ylim(10, 35)
ax.xaxis.set_ticks(np.arange(0,nt1_SSC+1,96))
ax.yaxis.set_ticks(np.linspace(10,35,6))
ax.set_xticklabels(['9/28','9/29','9/30','10/1'])
ax.xaxis.set_minor_locator(AutoMinorLocator(4))
ax.set_xlabel('Time', fontsize=12, fontname='Times New Roman', fontweight='bold')
ylabel = 'Suspended sediment ($\mathregular{mg}$ $\mathregular{l^{-1}}$)'
ax.set_ylabel(ylabel, fontsize=12, fontname='Times New Roman', fontweight='bold')
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(12) for label in labels]
[label.set_color('black') for label in labels]
[label.set_fontweight('bold') for label in labels]
ax.tick_params(which='major', direction='in', colors='xkcd:black', length=6, pad=8)
ax.tick_params(which='minor', direction='in', colors='xkcd:black')

#ax = axes[1][1]

plt.tight_layout()
fig.savefig('S4.png', dpi=300)
#fig.savefig('S4.pdf', dpi=600)
plt.show()
