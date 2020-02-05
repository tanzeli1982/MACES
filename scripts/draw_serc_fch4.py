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

try:
    nc = netcdf.netcdf_file('serc_fco2_fch4.nc','r')
    fch4_sat = np.array(nc.variables['fch4_sat'][:])
    fch4_unsat = np.array(nc.variables['fch4_unsat'][:])
    fch4 = np.array(nc.variables['fch4'][:])
finally:
    nc.close()

# plot
plt.clf()
fig, axes = plt.subplots(2, 1, figsize=(4,5))

colors = ['#1f77b4', '#d62728']

plt.style.use('default')
tt = np.arange(12*nyear)

ax = axes[0]
h1, = ax.plot(tt, fch4[:,0], color=colors[0], linestyle='-', linewidth=2, alpha=0.9)
h2, = ax.plot(tt, fch4[:,1], color=colors[1], linestyle='-', linewidth=2, alpha=0.9)
ax.legend([h1,h2], ['Coastal forest','Marshes'],
          numpoints=1, loc=0, ncol=1, framealpha=0.0,
          prop={'family':'Times New Roman', 'size':'medium'})
ax.set_xlim(0, nyear*12)
ax.set_ylim(-100, 2100)
ax.xaxis.set_ticks(np.arange(0,nyear*12+1,24))
ax.yaxis.set_ticks(np.linspace(0,2000,5))
ax.set_xticklabels(['2001','2003','2005','2007','2009','2011'])
ax.xaxis.set_minor_locator(AutoMinorLocator(2))
ylabel = '$\mathregular{{CH}_{4}}$ flux ($\mathregular{mg}$ ' + \
    '$\mathregular{{m}^{-2}}$ $\mathregular{{d}^{-1}}$)'
ax.set_ylabel(ylabel, fontsize=12, fontname='Times New Roman', color='black')
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(12) for label in labels]
[label.set_color('black') for label in labels]
ax.text(0.05, 0.9, 'A', transform=ax.transAxes, fontsize=14,
        fontname='Times New Roman', fontweight='bold')

ax = axes[1]
ax.plot(tt, fch4_unsat[:,0], color=colors[0], linestyle='-', linewidth=2, alpha=0.9)
ax.plot(tt, fch4_sat[:,1], color=colors[1], linestyle='-', linewidth=2, alpha=0.9)
ax.set_xlim(0, nyear*12)
ax.set_ylim(-100, 2100)
ax.xaxis.set_ticks(np.arange(0,nyear*12+1,24))
ax.yaxis.set_ticks(np.linspace(0,2000,5))
ax.set_xticklabels(['2001','2003','2005','2007','2009','2011'])
ax.xaxis.set_minor_locator(AutoMinorLocator(2))
ylabel = '$\mathregular{{CH}_{4}}$ flux ($\mathregular{mg}$ ' + \
    '$\mathregular{{m}^{-2}}$ $\mathregular{{d}^{-1}}$)'
ax.set_ylabel(ylabel, fontsize=12, fontname='Times New Roman', color='black')
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(12) for label in labels]
[label.set_color('black') for label in labels]
ax.text(0.05, 0.9, 'B', transform=ax.transAxes, fontsize=14,
        fontname='Times New Roman', fontweight='bold')

plt.tight_layout()
fig.savefig('serc_fch4.png', dpi=300)
#fig.savefig('serc_fch4.pdf', dpi=600)
plt.show()