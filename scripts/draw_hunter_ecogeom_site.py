#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 21 02:01:42 2020

Draw the comparison of mineral and OM accretion at Hunter Estuary.

@author: Zeli Tan
"""

import numpy as np
import matplotlib.pyplot as plt
import xml.etree.ElementTree as ET
from matplotlib.ticker import AutoMinorLocator
from netCDF4 import Dataset

# read sediment density and porosity of different mineral accretion models
models = ['F06', 'T03', 'KM12', 'F07', 'VDK05', 'DA07', 'M12']
xmlfile = '/Users/tanz151/Python_maces/src/optpar_minac.xml'
tree = ET.parse(xmlfile)
root = tree.getroot()
rhoSed = {}
porSed = {}
for key in models:
    findstr = "./group/[@id='" + key + 'MOD' + "']/entry"
    for entry in root.findall(findstr):
        if entry.get('id')=='rhoSed':
            rhoSed[key] = float(entry.get('value'))
        elif entry.get('id')=='porSed':
            porSed[key] = float(entry.get('value'))
            
min_accr_sim = {}
om_accr_sim = {}            
# read simulation
filename = '/Users/tanz151/Python_maces/src/maces_ecogeom_2004-01-01_2005-01-01_11099.nc'
try:
    nc = Dataset(filename, 'r')
    x = np.array(nc.variables['x'][:])
    zh = np.array(nc.variables['zh'][0,:])
    pft = np.array(nc.variables['pft'][0,:], dtype=np.int8)
    min_accr = (np.sum(8.64e7*np.array(nc.variables['Dsed'][:]),axis=0) - \
        np.sum(8.64e7*np.array(nc.variables['Esed'][:]),axis=0)) / \
        rhoSed['M12'] / (1.0-porSed['M12']) # mm/yr
    om_accr = 8.64e7*np.sum(np.array(nc.variables['DepOM'][:]),axis=0)  # g/m2/yr
finally:
    nc.close()
index0 = np.argmin(np.abs(zh))
x = x - x[index0]
indices = np.logical_or(pft==2, pft==5)
x = 1e-3 * x[indices]   # km
min_accr_sim['M12'] = min_accr[indices]
om_accr_sim['M12'] = om_accr[indices]

# plotting
plt.clf()
fig, axes = plt.subplots(2, 1, figsize=(6,8))

plt.style.use('default')

colors = ["#aee39a", "#643176", "#4be32e", "#e72fc2", "#518413", "#7540fc", 
          "#b3e61c"]
linestyles = ['-', '--', '-.', ':', '-', '--', '-.']

# comparison of simulted mineral accretion
ax = axes[0]
handles = []
for key in min_accr_sim:
    indx = len(handles)
    h, = ax.plot(x, min_accr_sim[key], color=colors[indx], 
                 linestyle=linestyles[indx], linewidth=3, alpha=1)
    handles.append(h)
legend = ax.legend(handles, list(min_accr_sim.keys()), numpoints=1, loc=1, 
                   prop={'family':'Times New Roman', 'size':'large'}, 
                   framealpha=0.0)
#ax.set_xlim(np.min(x), np.max(x))
#ax.set_ylim(0, 150)
#ax.xaxis.set_ticks(np.arange(0,nt_model+1,24))
#ax.yaxis.set_ticks(np.linspace(0,150,6))
ax.xaxis.set_minor_locator(AutoMinorLocator(5))
ax.set_xlabel('Distance ($\mathregular{km}$)', fontsize=12, 
              fontname='Times New Roman', color='black')
ylabel = 'Mineral accretion ($\mathregular{mm}$ $\mathregular{yr^{-1}}$)'
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
                 linestyle=linestyles[indx], linewidth=3, alpha=1)
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
fig.savefig('F13.png', dpi=300)
#fig.savefig('F13.pdf', dpi=600)
plt.show()