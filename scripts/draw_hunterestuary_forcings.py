#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 21 06:44:46 2020

Draw the forcing data of Hunter Estuary

@author: Zeli Tan
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from netCDF4 import Dataset
from datetime import date

t0 = []
t1 = []
t0.append( (date(2004,9,15) - date(2004,1,1)).days )
t1.append( (date(2004,9,20) - date(2004,1,1)).days )
t0.append( (date(2004,9,27) - date(2004,1,1)).days )
t1.append( (date(2004,10,2) - date(2004,1,1)).days )
t0.append( (date(2004,10,13) - date(2004,1,1)).days )
t1.append( (date(2004,10,18) - date(2004,1,1)).days )

h_vars = []
nt = []
# read data
filename = r'/Users/tanz151/Documents/Projects/TAI_BGC/Data/Hydrodynamics_obs/' + \
    r'HunterEstuary/force_h.nc'
try:
    nc = Dataset(filename,'r')
    for ii, t0v in enumerate(t0):
        t1v = t1[ii]
        h_var = 100 * np.array(nc.variables['h'][96*t0v:96*t1v])
        h_vars.append(h_var)
        nt.append(len(h_var))
finally:
    nc.close()
    
# plot
plt.clf()
fig, axes = plt.subplots(2, 2, figsize=(8,6))

plt.style.use('default')

ax = axes[0][0]
ax.plot(np.arange(nt[0]), h_vars[0], color='black', linestyle='-', linewidth=2, alpha=0.9)
ax.set_xlim(0, nt[0])
ax.set_ylim(-100, 100)
ax.xaxis.set_ticks(np.arange(0,nt[0],96))
#ax.yaxis.set_ticks(np.linspace(-40,60,6))
ax.set_xticklabels(['9/15','9/16','9/17','9/18','9/19'])
ax.xaxis.set_minor_locator(AutoMinorLocator(4))
ylabel = 'Water level ($\mathregular{cm}$ $\mathregular{AHD}$)'
ax.set_ylabel(ylabel, fontsize=12, fontname='Times New Roman', color='black')
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(12) for label in labels]
[label.set_color('black') for label in labels]
ax.tick_params(which='major', direction='in', colors='xkcd:black', length=6, pad=8)
ax.tick_params(which='minor', direction='in', colors='xkcd:black')

ax = axes[1][0]
ax.plot(np.arange(nt[1]), h_vars[1], color='black', linestyle='-', linewidth=2, alpha=0.9)
ax.set_xlim(0, nt[1])
ax.set_ylim(-100, 100)
ax.xaxis.set_ticks(np.arange(0,nt[1],96))
#ax.yaxis.set_ticks(np.linspace(3,18,6))
ax.set_xticklabels(['9/27','9/28','9/29','9/30','10/1'])
ax.xaxis.set_minor_locator(AutoMinorLocator(4))
ax.set_xlabel('Time', fontsize=12, fontname='Times New Roman', color='black')
ylabel = 'Water level ($\mathregular{cm}$ $\mathregular{AHD}$)'
ax.set_ylabel(ylabel, fontsize=12, fontname='Times New Roman', color='black')
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(12) for label in labels]
[label.set_color('black') for label in labels]
ax.tick_params(which='major', direction='in', colors='xkcd:black', length=6, pad=8)
ax.tick_params(which='minor', direction='in', colors='xkcd:black')

ax = axes[0][1]
ax.plot(np.arange(nt[2]), h_vars[2], color='black', linestyle='-', linewidth=2, alpha=0.9)
ax.set_xlim(0, nt[2])
ax.set_ylim(-100, 100)
ax.xaxis.set_ticks(np.arange(0,nt[2],96))
#ax.yaxis.set_ticks(np.linspace(3,18,6))
ax.set_xticklabels(['10/13','10/14','10/15','10/16','10/17'])
ax.xaxis.set_minor_locator(AutoMinorLocator(4))
ylabel = 'Water level ($\mathregular{cm}$ $\mathregular{AHD}$)'
ax.set_ylabel(ylabel, fontsize=12, fontname='Times New Roman', color='black')
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(12) for label in labels]
[label.set_color('black') for label in labels]
ax.tick_params(which='major', direction='in', colors='xkcd:black', length=6, pad=8)
ax.tick_params(which='minor', direction='in', colors='xkcd:black')

#ax = axes[1][1]

plt.tight_layout()
fig.savefig('hunter_forcing.png', dpi=300)
#fig.savefig('F3.pdf', dpi=600)
plt.show()
