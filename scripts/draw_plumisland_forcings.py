#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 21 06:44:46 2020

Draw the forcing data of Plum Island

@author: Zeli Tan
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import AutoMinorLocator
from netCDF4 import Dataset
from datetime import date

t1 = {}
t2 = {}
t1['d0'] = (date(2017,7,17) - date(2012,1,1)).days
t1['d1'] = (date(2017,7,24) - date(2012,1,1)).days
t2['d0'] = (date(2017,10,6) - date(2012,1,1)).days
t2['d1'] = (date(2017,10,13) - date(2012,1,1)).days

# read data
rdir = '/Users/tanz151/Documents/Projects/TAI_BGC/Drafts/Outputs/PlumIsland/'
filename = rdir + 'force_h.nc'
try:
    nc = Dataset(filename,'r')
    h_var1 = 100 * np.array(nc.variables['h'][t1['d0']*96:t1['d1']*96,0]) # cm
    h_var2 = 100 * np.array(nc.variables['h'][t2['d0']*96:t2['d1']*96,0]) # cm
finally:
    nc.close()
nt1_h = np.size(h_var1)
tt1_h = np.arange(nt1_h)
nt2_h = np.size(h_var2)
tt2_h = np.arange(nt2_h)
    
filename = rdir + 'force_U10.nc'
try:
    nc = Dataset(filename,'r')
    U10_var1 = np.array(nc.variables['U10'][t1['d0']*96:t1['d1']*96,0])
    U10_var2 = np.array(nc.variables['U10'][t2['d0']*96:t2['d1']*96,0])
finally:
    nc.close()
    
nt1_U10 = np.size(U10_var1)
tt1_U10 = np.arange(nt1_U10)
nt2_U10 = np.size(U10_var2)
tt2_U10 = np.arange(nt2_U10)

filename = rdir + 'force_SSC_revised.nc'
try:
    nc = Dataset(filename,'r')
    ssc_var1 = np.array(nc.variables['TSM'][t1['d0']*96:t1['d1']*96,0])
    ssc_var2 = np.array(nc.variables['TSM'][t2['d0']*96:t2['d1']*96,0])
finally:
    nc.close()
    
nt1_ssc = np.size(ssc_var1)
tt1_ssc = np.arange(nt1_ssc)
nt2_ssc = np.size(ssc_var2)
tt2_ssc = np.arange(nt2_ssc)
    
# plot
plt.clf()
fig = plt.figure(figsize=(8,9))
gs = gridspec.GridSpec(nrows=3, ncols=2)

plt.style.use('default')

ax = fig.add_subplot(gs[0,0])
ax.plot(tt1_h, h_var1, color='black', linestyle='-', linewidth=2, alpha=0.9)
ax.set_xlim(0, nt1_h-1)
ax.set_ylim(-200, 200)
ax.xaxis.set_ticks(np.arange(0,nt1_h,96))
ax.yaxis.set_ticks(np.linspace(-200,200,5))
ax.set_xticklabels(['7/17','7/18','7/19','7/20','7/21','7/22','7/23','7/24'])
ax.xaxis.set_minor_locator(AutoMinorLocator(4))
#ax.yaxis.set_minor_locator(AutoMinorLocator(2))
ylabel = 'Water level ($\mathregular{cm}$ $\mathregular{asl}$)'
ax.set_ylabel(ylabel, fontsize=12, fontname='Times New Roman', color='black', 
              fontweight='bold')
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(12) for label in labels]
[label.set_color('black') for label in labels]
[label.set_fontweight('bold') for label in labels]
ax.tick_params(which='major', direction='in', colors='xkcd:black', length=6, pad=8)
ax.tick_params(which='minor', direction='in', colors='xkcd:black')

ax = fig.add_subplot(gs[0,1])
ax.plot(tt2_h, h_var2, color='black', linestyle='-', linewidth=2, alpha=0.9)
ax.set_xlim(0, nt2_h-1)
ax.set_ylim(-200, 200)
ax.xaxis.set_ticks(np.arange(0,nt2_h,96))
ax.yaxis.set_ticks(np.linspace(-200,200,5))
ax.set_xticklabels(['10/6','10/7','10/8','10/9','10/10','10/11','10/12','10/13'])
ax.xaxis.set_minor_locator(AutoMinorLocator(4))
#ax.yaxis.set_minor_locator(AutoMinorLocator(2))
#ax.set_xlabel('Time', fontsize=12, fontname='Times New Roman', color='black')
#ylabel = 'Water level ($\mathregular{cm}$ $\mathregular{AHD}$)'
#ax.set_ylabel(ylabel, fontsize=12, fontname='Times New Roman', color='black')
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(12) for label in labels]
[label.set_color('black') for label in labels]
[label.set_fontweight('bold') for label in labels]
ax.tick_params(which='major', direction='in', colors='xkcd:black', length=6, pad=8)
ax.tick_params(which='minor', direction='in', colors='xkcd:black')

ax = fig.add_subplot(gs[1,0])
ax.plot(tt1_U10, U10_var1, color='black', linestyle='-', linewidth=2, alpha=0.9)
ax.set_xlim(0, nt1_U10-1)
ax.set_ylim(0, 9)
ax.xaxis.set_ticks(np.arange(0,nt1_U10,96))
ax.yaxis.set_ticks(np.linspace(0,8,5))
ax.set_xticklabels(['7/17','7/18','7/19','7/20','7/21','7/22','7/23','7/24'])
ax.xaxis.set_minor_locator(AutoMinorLocator(4))
#ax.yaxis.set_minor_locator(AutoMinorLocator(2))
#ax.set_xlabel('Time', fontsize=12, fontname='Times New Roman', color='black', 
#              fontweight='bold')
ylabel = 'Wind speed ($\mathregular{m}$ $\mathregular{{s}^{-1}}$)'
ax.set_ylabel(ylabel, fontsize=12, fontname='Times New Roman', color='black', 
              fontweight='bold')
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(12) for label in labels]
[label.set_color('black') for label in labels]
[label.set_fontweight('bold') for label in labels]
ax.tick_params(which='major', direction='in', colors='xkcd:black', length=6, pad=8)
ax.tick_params(which='minor', direction='in', colors='xkcd:black')

ax = fig.add_subplot(gs[1,1])
ax.plot(tt2_U10, U10_var2, color='black', linestyle='-', linewidth=2, alpha=0.9)
ax.set_xlim(0, nt2_U10-1)
ax.set_ylim(0, 9)
ax.xaxis.set_ticks(np.arange(0,nt2_U10,96))
ax.yaxis.set_ticks(np.linspace(0,8,5))
ax.set_xticklabels(['10/6','10/7','10/8','10/9','10/10','10/11','10/12','10/13'])
ax.xaxis.set_minor_locator(AutoMinorLocator(4))
#ax.yaxis.set_minor_locator(AutoMinorLocator(2))
#ax.set_xlabel('Time', fontsize=12, fontname='Times New Roman', color='black', 
#              fontweight='bold')
#ylabel = 'Wind speed ($\mathregular{m}$ $\mathregular{{s}^{-1}}$)'
#ax.set_ylabel(ylabel, fontsize=12, fontname='Times New Roman', color='black')
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(12) for label in labels]
[label.set_color('black') for label in labels]
[label.set_fontweight('bold') for label in labels]
ax.tick_params(which='major', direction='in', colors='xkcd:black', length=6, pad=8)
ax.tick_params(which='minor', direction='in', colors='xkcd:black')

ax = fig.add_subplot(gs[2,0])
ax.plot(tt1_ssc, ssc_var1, color='black', linestyle='-', linewidth=2, alpha=0.9)
ax.set_xlim(0, nt1_ssc-1)
ax.set_ylim(10, 50)
ax.xaxis.set_ticks(np.arange(0,nt1_ssc,96))
ax.yaxis.set_ticks(np.linspace(10,50,5))
ax.set_xticklabels(['7/17','7/18','7/19','7/20','7/21','7/22','7/23','7/24'])
ax.xaxis.set_minor_locator(AutoMinorLocator(4))
#ax.yaxis.set_minor_locator(AutoMinorLocator(2))
ax.set_xlabel('Time', fontsize=12, fontname='Times New Roman', color='black', 
              fontweight='bold')
ylabel = 'Suspended sediment ($\mathregular{mg}$ $\mathregular{{l}^{-1}}$)'
ax.set_ylabel(ylabel, fontsize=12, fontname='Times New Roman', color='black', 
              fontweight='bold')
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(12) for label in labels]
[label.set_color('black') for label in labels]
[label.set_fontweight('bold') for label in labels]
ax.tick_params(which='major', direction='in', colors='xkcd:black', length=6, pad=8)
ax.tick_params(which='minor', direction='in', colors='xkcd:black')

ax = fig.add_subplot(gs[2,1])
ax.plot(tt2_ssc, ssc_var2, color='black', linestyle='-', linewidth=2, alpha=0.9)
ax.set_xlim(0, nt2_ssc-1)
ax.set_ylim(10, 50)
ax.xaxis.set_ticks(np.arange(0,nt2_ssc,96))
ax.yaxis.set_ticks(np.linspace(10,50,5))
ax.set_xticklabels(['10/6','10/7','10/8','10/9','10/10','10/11','10/12','10/13'])
ax.xaxis.set_minor_locator(AutoMinorLocator(4))
#ax.yaxis.set_minor_locator(AutoMinorLocator(2))
ax.set_xlabel('Time', fontsize=12, fontname='Times New Roman', color='black', 
              fontweight='bold')
#ylabel = 'Wind speed ($\mathregular{m}$ $\mathregular{{s}^{-1}}$)'
#ax.set_ylabel(ylabel, fontsize=12, fontname='Times New Roman', color='black')
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(12) for label in labels]
[label.set_color('black') for label in labels]
[label.set_fontweight('bold') for label in labels]
ax.tick_params(which='major', direction='in', colors='xkcd:black', length=6, pad=8)
ax.tick_params(which='minor', direction='in', colors='xkcd:black')

plt.tight_layout()
fig.savefig('S3.png', dpi=300)
#fig.savefig('S3.pdf', dpi=600)
plt.show()
