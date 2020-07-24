#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 21 01:18:18 2020

Draw the OM accretion related plot at Venice Lagoon.

@author: Zeli Tan
"""

import numpy as np
import matplotlib.pyplot as plt
import xml.etree.ElementTree as ET
from netCDF4 import Dataset
from matplotlib.ticker import AutoMinorLocator
from datetime import date

om_models = ['M12', 'DA07', 'KM12', 'K16']
min_models = ['F06', 'T03', 'KM12', 'M12', 'F07', 'VDK05', 'DA07']
rhoSed = {}
porSed = {}
xmlfile = '/Users/tanz151/Python_maces/src/optpar_minac.xml'
tree = ET.parse(xmlfile)
root = tree.getroot()
for key in min_models:
    findstr = "./group/[@id='" + key + 'MOD' + "']/entry"
    for entry in root.findall(findstr):
        if entry.get('id')=='rhoSed':
            rhoSed[key] = float(entry.get('value'))
        elif entry.get('id')=='porSed':
            porSed[key] = float(entry.get('value'))

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
x_marsh = 1e-3*x[indices_marsh]   # km

indice_obs = np.logical_and(x>=0, x<=0.2)

nx = np.size(x)
dx = np.zeros(nx)
for ii in range(nx):
    if ii==0:
        dx[ii] = 0.5 * (x[ii+1] - x[ii])
    elif ii==nx-1:
        dx[ii] = 0.5 * (x[ii] - x[ii-1])
    else:
        dx[ii] = 0.5 * (x[ii+1] - x[ii-1])

for model in min_models:
    filename = rdir + 'maces_ecogeom_2002-01-01_2004-01-01_466.' + model + 'DA07.nc'
    try:
        nc = Dataset(filename,'r')
        Esed = 0.5*8.64e7*np.sum(np.array(nc.variables['Esed'][:]),axis=0)
        Dsed = 0.5*8.64e7*np.sum(np.array(nc.variables['Dsed'][:]),axis=0)
    finally:
        nc.close()
    minac_sim = np.sum((Dsed[indice_obs]-Esed[indice_obs])*dx[indice_obs]) / \
        np.sum(dx[indice_obs]) / rhoSed[model] / (1.0-porSed[model]) # mm/yr
    print('MINAC MODEL: ', model, ', ', minac_sim)

for model in om_models:
    filename = rdir + 'maces_ecogeom_2002-01-01_2004-01-01_466.M12' + model + '.nc'
    try:
        nc = Dataset(filename,'r')
        Bag = 1e3*np.mean(np.array(nc.variables['Bag'][day0:day1,:]),axis=0)    # g/m2
        om_accr = 0.5*8.64e7*np.sum(np.array(nc.variables['DepOM'][:]),axis=0)  # g/m2/yr
    finally:
        nc.close()
    omac_sim = np.sum(om_accr[indice_obs]*dx[indice_obs]) / np.sum(dx[indice_obs])
    print('OMAC MODEL: ', model, ', ', omac_sim)
    bg_sim[model] = Bag[indices_marsh]
    om_accr_sim[model] = om_accr[indices_marsh]

# plotting
plt.clf()
fig, axes = plt.subplots(2, 1, figsize=(6,8))

plt.style.use('default')

colors = ['#7b85d4', '#f37738', '#83c995', '#d7369e', '#c4c9d8', '#859795',
          '#e9d043', '#ad5b50', '#e377c2']
#colors = ["#aee39a", "#643176", "#4be32e", "#e72fc2", "#518413", "#7540fc", 
#          "#b3e61c"]
linestyles = ['-', '--', '-.', ':', '-', '--', '-.']

# comparison of aboveground biomass
ax = axes[0]
handles = []
for key in bg_sim:
    indx = len(handles)
    h, = ax.plot(x_marsh, bg_sim[key], color=colors[indx], 
                 linestyle=linestyles[indx], linewidth=2, alpha=1)
    handles.append(h)
legend = ax.legend(handles, list(bg_sim.keys()), numpoints=1, loc=1, 
                   prop={'family':'Times New Roman', 'size':'large'}, 
                   framealpha=0.0)
ax.set_xlim(0, 1.5)
ax.set_ylim(500, 2500)
ax.xaxis.set_ticks(np.linspace(0,1.5,6))
ax.yaxis.set_ticks(np.linspace(500,2500,6))
ax.xaxis.set_minor_locator(AutoMinorLocator(5))
#ax.set_xlabel('Distance ($\mathregular{km}$)', fontsize=12, 
#              fontname='Times New Roman', color='black')
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
    h, = ax.plot(x_marsh, om_accr_sim[key], color=colors[indx], 
                 linestyle=linestyles[indx], linewidth=2, alpha=1)
    handles.append(h)
ax.set_xlim(0, 1.5)
ax.set_ylim(0, 300)
ax.xaxis.set_ticks(np.linspace(0,1.5,6))
ax.yaxis.set_ticks(np.linspace(0,300,7))
ax.xaxis.set_minor_locator(AutoMinorLocator(5))
ax.set_xlabel('Distance ($\mathregular{km}$)', fontsize=12, 
              fontname='Times New Roman', color='black')
ylabel = 'OM accretion ($\mathregular{gC}$ $\mathregular{m^{-2}}$ $\mathregular{yr^{-1}}$)'
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
fig.savefig('F6.pdf', dpi=600)
plt.show()
