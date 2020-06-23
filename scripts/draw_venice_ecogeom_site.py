#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 21 01:18:18 2020

Draw the OM accretion related plot at Venice Lagoon.

@author: Zeli Tan
"""

import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from matplotlib.ticker import AutoMinorLocator
from datetime import date

models = ['M12', 'DA07', 'KM12', 'K16']
om_accr_sim = {}
bg_sim = {}

day0 = (date(2002,7,1) - date(2002,1,1)).days
day1 = (date(2002,8,1) - date(2002,1,1)).days
# read simulation outputs
rdir = '/Users/tanz151/Documents/Projects/TAI_BGC/Data/Hydrodynamics_obs/VeniceLagoon/Outputs/'
filename = rdir + 'maces_ecogeom_2002-01-01_2004-01-01_466.M12M12.nc'
try:
    nc = Dataset(filename,'r')
    x = np.array(nc.variables['x'][:])
    zh = np.array(nc.variables['zh'][0,:])
    pft = np.array(nc.variables['pft'][0,:], dtype=np.int8)
finally:
    nc.close()
index0 = np.argmin(np.abs(zh))
x = x - x[index0]
# extract salt marsh zone
indices_marsh = pft==2
x = 1e-3*x[indices_marsh]   # km

for model in models:
    filename = rdir + 'maces_ecogeom_2002-01-01_2004-01-01_466.M12' + model + '.nc'
    try:
        nc = Dataset(filename,'r')
        Bag = 1e3*np.mean(np.array(nc.variables['Bag'][day0:day1,:]),axis=0)    # g/m2
        om_accr = 0.5*8.64e7*np.sum(np.array(nc.variables['DepOM'][:]),axis=0)  # g/m2/yr
    finally:
        nc.close()
    bg_sim[model] = Bag[indices_marsh]
    om_accr_sim[model] = om_accr[indices_marsh]

# plotting
plt.clf()
fig, axes = plt.subplots(2, 1, figsize=(6,8))

plt.style.use('default')

colors = ["#aee39a", "#643176", "#4be32e", "#e72fc2", "#518413", "#7540fc", 
          "#b3e61c"]
linestyles = ['-', '--', '-.', ':', '-', '--', '-.']

# comparison of aboveground biomass
ax = axes[0]
handles = []
for key in bg_sim:
    indx = len(handles)
    h, = ax.plot(x, bg_sim[key], color=colors[indx], 
                 linestyle=linestyles[indx], linewidth=2, alpha=1)
    handles.append(h)
legend = ax.legend(handles, list(bg_sim.keys()), numpoints=1, loc=1, 
                   prop={'family':'Times New Roman', 'size':'large'}, 
                   framealpha=0.0)
#ax.set_xlim(np.min(x), np.max(x))
#ax.set_ylim(0, 150)
#ax.xaxis.set_ticks(np.arange(0,nt_model+1,24))
#ax.yaxis.set_ticks(np.linspace(0,150,6))
ax.xaxis.set_minor_locator(AutoMinorLocator(5))
ax.set_xlabel('Distance ($\mathregular{km}$)', fontsize=12, 
              fontname='Times New Roman', color='black')
ylabel = 'Aboveground biomass ($\mathregular{g}$ $\mathregular{m^{-2}}$)'
ax.set_ylabel(ylabel, fontsize=12, fontname='Times New Roman', color='black')
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(12) for label in labels]
[label.set_color('black') for label in labels]
ax.text(0.05, 0.93, 'a', transform=ax.transAxes, fontsize=20,
        fontname='Times New Roman', fontweight='bold')
ax.tick_params(which='major', direction='in', colors='xkcd:black', length=6, pad=8)
ax.tick_params(which='minor', direction='in', colors='xkcd:black')

# comparison of simulated OM accretion
ax = axes[1]
handles = []
for key in om_accr_sim:
    indx = len(handles)
    h, = ax.plot(x, om_accr_sim[key], color=colors[indx], 
                 linestyle=linestyles[indx], linewidth=2, alpha=1)
    handles.append(h)
#ax.set_xlim(-1, 2)
#ax.set_ylim(0, 150)
#ax.xaxis.set_ticks(np.linspace(-1,2,7))
#ax.yaxis.set_ticks(np.linspace(0,150,6))
ax.xaxis.set_minor_locator(AutoMinorLocator(5))
ax.set_xlabel('Distance ($\mathregular{km}$)', fontsize=12, 
              fontname='Times New Roman', color='black')
ylabel = 'OM accretion ($\mathregular{g}$ $\mathregular{m^{-2}}$ $\mathregular{yr^{-1}}$)'
ax.set_ylabel(ylabel, fontsize=12, fontname='Times New Roman', color='black')
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(12) for label in labels]
[label.set_color('black') for label in labels]
ax.text(0.05, 0.93, 'b', transform=ax.transAxes, fontsize=20,
        fontname='Times New Roman', fontweight='bold')
ax.tick_params(which='major', direction='in', colors='xkcd:black', length=6, pad=8)
ax.tick_params(which='minor', direction='in', colors='xkcd:black')
    
plt.tight_layout()
fig.savefig('F6.png', dpi=300)
#fig.savefig('F6.pdf', dpi=600)
plt.show()