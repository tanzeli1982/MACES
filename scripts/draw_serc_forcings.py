#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Mar  9 15:19:45 2019

Compare the simulated CH4 and CO2 emissions at two SERC sites

@author: Zeli Tan
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.io import netcdf
from matplotlib.ticker import AutoMinorLocator

# read ten year simulations
years = np.arange(2001,2011)
nyear = np.size(years)
#days = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

try:
    nc = netcdf.netcdf_file('serc_fco2_fch4.nc','r')
    tair = np.array(nc.variables['tair'][:])
    tgnd = np.array(nc.variables['tgnd'][:])
    tsoi = np.array(nc.variables['tsoi'][:])
    prcp = np.array(nc.variables['prcp'][:])
finally:
    nc.close()

# plot
plt.clf()
fig, axes = plt.subplots(2, 1, figsize=(5,6))

colors = ['#1f77b4', '#d62728']

plt.style.use('default')
tt = np.arange(12*nyear)

ax = axes[0]
#ax.plot(tt, tsoi[:,0]-tsoi[:,1], color=colors[0], linestyle='-', linewidth=2, alpha=0.9)
ax.plot(tt, tsoi[:,0], color=colors[0], linestyle='-', linewidth=2, alpha=0.9)
ax.plot(tt, tsoi[:,1], color=colors[1], linestyle='-', linewidth=2, alpha=0.9)
#ax.legend(['Upland forest habitat','Marsh habitat'], numpoints=1,
#          loc=2, ncol=1, framealpha=0.0, 
#          prop={'family':'Times New Roman', 'size':'large', 'weight':'bold'})
ax.set_xlim(0, nyear*12)
#ax.set_ylim(-7.5, 10)
ax.xaxis.set_ticks(np.arange(0,nyear*12+1,24))
#ax.yaxis.set_ticks(np.linspace(-7.5,10,8))
ax.set_xticklabels(['2001','2003','2005','2007','2009','2011'])
ax.xaxis.set_minor_locator(AutoMinorLocator(2))
ylabel = 'Air temperature (celsius)'
ax.set_ylabel(ylabel, fontsize=14, fontname='Times New Roman', color='black')
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(12) for label in labels]
[label.set_color('black') for label in labels]

ax = axes[1]
ax.plot(tt, prcp[:,0], color=colors[0], linestyle='-', linewidth=2, alpha=0.9)
ax.plot(tt, prcp[:,1], color=colors[1], linestyle='-', linewidth=2, alpha=0.9)
ax.set_xlim(0, nyear*12)
#ax.set_ylim(-100, 100)
ax.xaxis.set_ticks(np.arange(0,nyear*12+1,24))
#ax.yaxis.set_ticks(np.linspace(-100,100,6))
ax.set_xticklabels(['2001','2003','2005','2007','2009','2011'])
ax.xaxis.set_minor_locator(AutoMinorLocator(2))
ylabel = 'Rainfall ($\mathregular{mm}$ $\mathregular{{d}^{-1}}$)'
ax.set_ylabel(ylabel, fontsize=14, fontname='Times New Roman', color='black')
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(12) for label in labels]
[label.set_color('black') for label in labels]

plt.tight_layout()
fig.savefig('serc_forcings.png', dpi=300)
#fig.savefig('serc_fco2_fch4.pdf', dpi=600)
plt.show()