#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  3 14:03:30 2020

Draw the dynamics of water flow velocity on the TAI platform

@author: Zeli Tan
"""

import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

# read the simulated water flow velocity
filename = '/Users/tanz151/Python_maces/src/out_hydro_2001-01-01_2001-01-06.nc'
try:
    nc = Dataset(filename,'r')
    xv = np.array(nc.variables['x'][:][0])
    Uw = np.array(nc.variables['U'][:][0])
finally:
    nc.close()
ntime = np.shape(Uw)[0]

tt = -np.arange(ntime)
    
# plot
plt.clf()
fig, ax = plt.subplots(figsize=(8,10))

plt.style.use('default')

cf = ax.contourf(xv, tt, Uw, 10, cmap='RdYlBu_r')
ax.set_xlim([xv[0], xv[-1]])
ax.set_ylim([tt[-1], tt[0]])
cbar = fig.colorbar(cf, ax=ax, orientation='horizontal', pad=0.05)
ylabel = 'Water flow ($\mathregular{m}$ $\mathregular{{s}^{-1}}$)'
cbar.set_label(ylabel, fontsize=12, fontname='Times New Roman', 
               labelpad=-50)
#cbar.ax.xaxis.set_minor_locator(AutoMinorLocator(2))
labels = cbar.ax.get_xticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(12) for label in labels]

labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(12) for label in labels]

plt.tight_layout()
fig.savefig('flow.png', dpi=300)
#fig.savefig('F1.pdf', dpi=600)
plt.show()