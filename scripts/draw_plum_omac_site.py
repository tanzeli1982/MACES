#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 29 21:12:25 2020

Plot the comparison of OM accretion at Plum Island 

@author: Zeli Tan
"""

import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from matplotlib.ticker import AutoMinorLocator

models = ['M12', 'DA07', 'KM12', 'K16']
omac_sim_mean = {'Data': 69.9}
omac_sim_std = {'Data': 9.4}

# read grid
rdir = '/Users/tanz151/Documents/Projects/TAI_BGC/Data/Hydrodynamics_obs/PlumIsland/Outputs/'
filename = rdir + 'maces_ecogeom_2017-01-01_2019-01-01_4097.T03DA07.nc'
try:
    nc = Dataset(filename, 'r')
    x = np.array(nc.variables['x'][:])
    zh = np.array(nc.variables['zh'][0,:])
    pft = np.array(nc.variables['pft'][0,:], dtype=np.int8)
finally:
    nc.close()
nx = np.size(x)
dx = np.zeros(nx)
for ii in range(nx):
    if ii==0:
        dx[ii] = 0.5 * (x[ii+1] - x[ii])
    elif ii==nx-1:
        dx[ii] = 0.5 * (x[ii] - x[ii-1])
    else:
        dx[ii] = 0.5 * (x[ii+1] - x[ii-1])
index0 = np.argmin(np.abs(zh))
x = x - x[index0]
indice_obs = np.logical_and(zh>=0, zh<=1.5)
x_obs = x[indice_obs]
indices_marsh = pft==2
x_marsh = x[indices_marsh]
for model in models:
    filename = rdir + 'maces_ecogeom_2017-01-01_2019-01-01_4097.T03' + model + '.nc'
    try:
        nc = Dataset(filename,'r')
        om_accr = 0.5*8.64e7*np.sum(np.array(nc.variables['DepOM'][:]),axis=0)  # g/m2/yr
    finally:
        nc.close()
    omac_sim_mean[model] = np.sum(om_accr[indice_obs]*dx[indice_obs]) / np.sum(dx[indice_obs])
    omac_sim_std[model] = np.std(om_accr[indice_obs])
    
# plotting
plt.clf()
fig, ax = plt.subplots(figsize=(6,4.5))

plt.style.use('default')

colors = ['#d8dcd6', '#7b85d4', '#f37738', '#83c995', '#d7369e', '#c4c9d8', '#859795',
          '#e9d043', '#ad5b50', '#e377c2']
#colors = ["#aee39a", "#643176", "#4be32e", "#e72fc2", "#518413", "#7540fc", 
#          "#b3e61c"]
linestyles = ['-', '--', '-.', ':', '-', '--', '-.']

nbar = len(omac_sim_mean)
xpos = np.arange(1,nbar+1)
hbar = ax.bar(xpos, list(omac_sim_mean.values()), yerr=list(omac_sim_std.values()), 
              align='center', width=0.8, color=colors[:nbar], ecolor='black', 
              capstyle='butt', capsize=2, alpha=1.0)
ax.set_xlim(0, nbar+1)
ax.set_ylim(0, 105)
ax.xaxis.set_ticks(np.arange(1,nbar+1,1))
ax.yaxis.set_ticks(np.linspace(0,100,6))
ax.xaxis.set_minor_locator(AutoMinorLocator(5))
ax.set_xticklabels(omac_sim_mean.keys())
ylabel = 'OM accretion ($\mathregular{g}$ $\mathregular{m^{-2}}$ $\mathregular{yr^{-1}}$)'
ax.set_ylabel(ylabel, fontsize=12, fontname='Times New Roman', color='black')
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(12) for label in labels]
[label.set_color('black') for label in labels]
ax.tick_params(which='major', direction='in', colors='xkcd:black', length=6, pad=8)
ax.tick_params(which='minor', direction='in', colors='xkcd:black')

plt.tight_layout()
fig.savefig('S2.png', dpi=300)
fig.savefig('S2.pdf', dpi=600)
plt.show()