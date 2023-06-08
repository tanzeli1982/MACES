#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 21 02:01:42 2020

Draw the comparison of mineral and OM accretion at Hunter Estuary.

@author: Zeli Tan
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import xml.etree.ElementTree as ET
from matplotlib.ticker import AutoMinorLocator
from netCDF4 import Dataset
from datetime import date

# read sediment density and porosity of different mineral accretion models
min_models = ['F06', 'T03', 'KM12', 'M12', 'F07']#, 'VDK05', 'DA07']
om_models = ['M12', 'DA07', 'KM12', 'K16']
case_min = 'F07'
case_om = 'M12'

rdir = '/Users/tanz151/Documents/Projects/TAI_BGC/Drafts/Outputs/HunterEstuary/'
xmlfile = rdir + 'optpar_minac.xml'
tree = ET.parse(xmlfile)
root = tree.getroot()
rhoSed = {}
porSed = {}
rhoOM = {}
for key in min_models:
    findstr = "./group/[@id='" + key + 'MOD' + "']/entry"
    for entry in root.findall(findstr):
        if entry.get('id')=='rhoSed':
            rhoSed[key] = float(entry.get('value'))
        elif entry.get('id')=='porSed':
            porSed[key] = float(entry.get('value'))
            
xmlfile = rdir + 'optpar_omac.xml'
tree = ET.parse(xmlfile)
root = tree.getroot()
for key in om_models:
    findstr = "./group/[@id='" + key + 'MOD' + "']/entry"
    for entry in root.findall(findstr):
        if entry.get('id')=='rhoOM':
            rhoOM[key] = float(entry.get('value'))
            
om_accr_sim = {}
min_accr_sim = {}
tot_accr_sim = {}
bg_sim = {}

day0 = (date(2004,7,1) - date(2004,1,1)).days
day1 = (date(2004,8,1) - date(2004,1,1)).days
site_dem = 0.56

filename = rdir + 'maces_ecogeom_2004-01-01_2005-01-01_' + case_min + '%' + \
    case_om + '_11099.nc'
try:
    nc = Dataset(filename,'r')
    x = np.array(nc.variables['x'][:])
    zh = np.array(nc.variables['zh'][0,:])
    pft = np.array(nc.variables['pft'][0,:], dtype=np.int8)
finally:
    nc.close()
index0 = np.argmin(np.abs(zh))
x = x - x[index0]
nx = len(x)

zh_ref = zh - site_dem
zh_ref[zh_ref<0] = 1e20
index_obs = []
index_obs.append(np.argmin(zh_ref))
zh_ref = zh - site_dem
zh_ref[zh_ref>=0] = -1e20
index_obs.append(np.argmax(zh_ref))
print("index_obs: ", index_obs)

for model in min_models:
    for omodel in om_models:
        filename = rdir + 'maces_ecogeom_2004-01-01_2005-01-01_' + model + \
            '%' + omodel + '_11099.nc'
        try:
            nc = Dataset(filename,'r')
            Esed = 8.64e7*np.sum(np.array(nc.variables['Esed'][:]),axis=0)  # kg/m2/yr
            Dsed = 8.64e7*np.sum(np.array(nc.variables['Dsed'][:]),axis=0)  # kg/m2/yr
            Bag = 1e3*np.mean(np.array(nc.variables['Bag'][day0:day1,:]),axis=0)    # g/m2
            Bmax = 1e3*np.max(np.array(nc.variables['Bag'][day0:day1,:]),axis=0)    # g/m2
            om_accr = 8.64e7*np.sum(np.array(nc.variables['DepOM'][:]),axis=0)  # g/m2/yr
        finally:
            nc.close()
        minac_sim = (np.interp(site_dem, zh, Dsed) - np.interp(site_dem, zh, Esed)) / \
            rhoSed[model] / (1.0-porSed[model]) # mm/yr
        omac_sim = np.interp(site_dem, zh, om_accr)
        key = model + '%' + omodel
        min_accr_sim[key] = (Dsed - Esed) / rhoSed[model] / (1.0-porSed[model]) # mm/yr
        if omodel==case_om:
            print('MINAC MODEL: ', model, ', ', minac_sim)
        if model==case_min:
            print('OMAC MODEL: ', omodel, ', ', omac_sim)
            print('Bmax: ', np.max(Bag))
        bg_sim[key] = Bag
        om_accr_sim[key] = om_accr
        tot_accr_sim[key] = min_accr_sim[key] + om_accr_sim[key] / 0.44 / \
            rhoOM[omodel] / (1.0-porSed[model]) # mm/yr

# plotting
plt.clf()
fig = plt.figure(figsize=(8.5,7))

gs = gridspec.GridSpec(nrows=2, ncols=2)

plt.style.use('default')

colors = ['#7b85d4', '#f37738', '#83c995', '#d7369e', '#c4c9d8', '#859795',
          '#e9d043', '#ad5b50', '#e377c2']
#colors = ["#aee39a", "#643176", "#4be32e", "#e72fc2", "#518413", "#7540fc", 
#          "#b3e61c"]
linestyles = ['-', '--', '-.', ':', '-', '--', '-.']

# aboveground biomass
ax = fig.add_subplot(gs[0,0])
handles = []
for indx, omodel in enumerate(om_models):
    key = case_min + '%' + omodel
    h, = ax.plot(x, bg_sim[key], color=colors[indx], 
                 linestyle=linestyles[indx], linewidth=2, alpha=1)
    handles.append(h)
legend = ax.legend(handles, om_models, numpoints=1, loc='upper center', 
                   prop={'family':'Times New Roman', 'size':'large', 'weight': 'bold'}, 
                   framealpha=0.0)
ax.set_xlim(0, 80)
ax.set_ylim(600, 1600)
ax.xaxis.set_ticks(np.linspace(0,80,5))
ax.yaxis.set_ticks(np.linspace(600,1600,6))
ax.xaxis.set_minor_locator(AutoMinorLocator(4))
ax.set_xlabel('Distance ($\mathregular{m}$)', fontsize=12, 
              fontname='Times New Roman', color='black', fontweight='bold')
ylabel = 'Aboveground biomass ($\mathregular{g}$ $\mathregular{m^{-2}}$)'
ax.set_ylabel(ylabel, fontsize=12, fontname='Times New Roman', color='black', 
              fontweight='bold')
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(12) for label in labels]
[label.set_color('black') for label in labels]
[label.set_fontweight('bold') for label in labels]
ax.text(0.05, 0.9, 'a', transform=ax.transAxes, fontsize=16,
        fontname='Times New Roman', fontweight='bold')
ax.tick_params(which='major', direction='in', colors='xkcd:black', length=6, pad=8)
ax.tick_params(which='minor', direction='in', colors='xkcd:black')

# comparison of simulated OM accretion
ax = fig.add_subplot(gs[0,1])
handles = []
for indx, omodel in enumerate(om_models):
    key = case_min + '%' + omodel
    h, = ax.plot(x, om_accr_sim[key], color=colors[indx], 
                 linestyle=linestyles[indx], linewidth=2, alpha=1)
    handles.append(h)
legend = ax.legend(handles, om_models, numpoints=1, loc='upper center', 
                   prop={'family':'Times New Roman', 'size':'large', 'weight': 'bold'}, 
                   framealpha=0.0)
ax.plot(np.interp(site_dem, zh, x), 105, color='black', marker='*', mec='black', 
        mfc='black', ms=15, alpha=1.0)
ax.set_xlim(0, 80)
ax.set_ylim(50, 250)
ax.xaxis.set_ticks(np.linspace(0,80,5))
ax.yaxis.set_ticks(np.linspace(50,250,5))
ax.xaxis.set_minor_locator(AutoMinorLocator(4))
ax.set_xlabel('Distance ($\mathregular{m}$)', fontsize=12, 
              fontname='Times New Roman', color='black', fontweight='bold')
ylabel = 'OM accretion ($\mathregular{gC}$ $\mathregular{m^{-2}}$ $\mathregular{yr^{-1}}$)'
ax.set_ylabel(ylabel, fontsize=12, fontname='Times New Roman', color='black', 
              fontweight='bold')
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(12) for label in labels]
[label.set_color('black') for label in labels]
[label.set_fontweight('bold') for label in labels]
ax.text(0.05, 0.9, 'b', transform=ax.transAxes, fontsize=16,
        fontname='Times New Roman', fontweight='bold')
ax.tick_params(which='major', direction='in', colors='xkcd:black', length=6, pad=8)
ax.tick_params(which='minor', direction='in', colors='xkcd:black')

# comparison of simulted mineral accretion
ax = fig.add_subplot(gs[1,0])
handles = []
for indx, model in enumerate(min_models):
    key = model + '%' + case_om
    h, = ax.plot(x, min_accr_sim[key], color=colors[indx], 
                 linestyle=linestyles[indx], linewidth=2, alpha=1)
    handles.append(h)
ax.plot(np.interp(site_dem, zh, x), 3.66, color='black', marker='*', mec='black', 
        mfc='black', ms=15, alpha=1.0)
legend = ax.legend(handles, min_models, numpoints=1, loc="upper right", 
                   prop={'family':'Times New Roman', 'size':'large', 'weight': 'bold'}, 
                   framealpha=0.0, ncol=2)
ax.set_xlim(0, 80)
ax.set_ylim(0, 50)
ax.xaxis.set_ticks(np.linspace(0,80,5))
ax.yaxis.set_ticks(np.linspace(0,50,6))
ax.xaxis.set_minor_locator(AutoMinorLocator(4))
ax.set_xlabel('Distance ($\mathregular{m}$)', fontsize=12, 
              fontname='Times New Roman', color='black', fontweight='bold')
ylabel = 'Mineral accretion ($\mathregular{mm}$ $\mathregular{{yr}^{-1}}$)'
ax.set_ylabel(ylabel, fontsize=12, fontname='Times New Roman', color='black', 
              fontweight='bold')
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(12) for label in labels]
[label.set_color('black') for label in labels]
[label.set_fontweight('bold') for label in labels]
ax.text(0.05, 0.9, 'c', transform=ax.transAxes, fontsize=16,
        fontname='Times New Roman', fontweight='bold')
ax.tick_params(which='major', direction='in', colors='xkcd:black', length=6, pad=8)
ax.tick_params(which='minor', direction='in', colors='xkcd:black')

# ensemble estimates of total accretion
nmodel = len(min_models) * len(om_models)
tot_accr_ensemble = np.zeros((nmodel,nx))
fom_accr_ensemble = np.zeros((nmodel,nx))
indx = 0
for model in min_models:
    for omodel in om_models:
        key = model + '%' + omodel
        tot_accr_ensemble[indx] = tot_accr_sim[key]
        fom_accr_ensemble[indx] = 100 * (1.0 - min_accr_sim[key]/tot_accr_sim[key])
        indx = indx + 1

ax = fig.add_subplot(gs[1,1])
ax.plot(x, np.percentile(tot_accr_ensemble,50,axis=0), color='black', 
        linestyle='-', linewidth=2, alpha=0.9)
ax.fill_between(x, np.percentile(tot_accr_ensemble,25,axis=0), 
                np.percentile(tot_accr_ensemble,75,axis=0), alpha=0.3, 
                facecolor='black')
ax.set_xlim(0, 80)
ax.set_ylim(0, 50)
ax.xaxis.set_ticks(np.linspace(0,80,5))
ax.yaxis.set_ticks(np.linspace(0,50,6))
ax.xaxis.set_minor_locator(AutoMinorLocator(4))
ax.set_xlabel('Distance ($\mathregular{m}$)', fontsize=12, 
              fontname='Times New Roman', color='black', fontweight='bold')
ylabel = 'Total accretion ($\mathregular{mm}$ $\mathregular{{yr}^{-1}}$)'
ax.set_ylabel(ylabel, fontsize=12, fontname='Times New Roman', color='black', 
              fontweight='bold')
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(12) for label in labels]
[label.set_color('black') for label in labels]
[label.set_fontweight('bold') for label in labels]
ax.text(0.05, 0.9, 'd', transform=ax.transAxes, fontsize=16,
        fontname='Times New Roman', fontweight='bold')
ax.tick_params(which='major', direction='in', colors='xkcd:black', length=6, pad=8)
ax.tick_params(which='minor', direction='in', colors='xkcd:black')

axInv = ax.twinx()
axInv.plot(x, np.percentile(fom_accr_ensemble,50,axis=0), color='C3', 
           linestyle='-', linewidth=2, alpha=0.9)
axInv.fill_between(x, np.percentile(fom_accr_ensemble,25,axis=0), 
                   np.percentile(fom_accr_ensemble,75,axis=0), alpha=0.3, 
                   facecolor='C3')
axInv.set_xlim(0, 80)
axInv.set_ylim(0, 100)
axInv.xaxis.set_ticks(np.arange(0,80,5))
axInv.yaxis.set_ticks(np.linspace(0,100,6))
axInv.tick_params(which='major', direction='in', colors='xkcd:black', length=6, pad=4)
axInv.tick_params(which='minor', direction='in', colors='xkcd:black')
axInv.tick_params(axis='y', colors='C3')
axInv.spines['right'].set_color('C3')
axInv.set_ylabel('Percentage of OM accretion', color='C3', fontsize=11, 
                 fontname='Times New Roman', fontweight='bold')
labels = axInv.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(12) for label in labels]
[label.set_color('C3') for label in labels]
[label.set_fontweight('bold') for label in labels]

ax.xaxis.set_ticks(np.linspace(0,80,5))
ax.yaxis.set_ticks(np.linspace(0,50,6))
ax.xaxis.set_minor_locator(AutoMinorLocator(4))
    
plt.tight_layout()
fig.savefig('F12.png', dpi=300)
#fig.savefig('F12.jpg', dpi=600)
plt.show()
