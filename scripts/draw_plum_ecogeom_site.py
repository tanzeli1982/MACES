#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 09:48:40 2020

Compare the simulated and observed sediment accretion and biomass.

@author: Zeli Tan
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import xml.etree.ElementTree as ET
from matplotlib.ticker import AutoMinorLocator
from netCDF4 import Dataset
from datetime import date

# read sediment density and porosity of different mineral accretion models
min_models = ['F06', 'T03', 'KM12', 'M12', 'F07', 'VDK05', 'DA07']
om_models = ['M12', 'DA07', 'KM12', 'K16']
case_min = 'F07'
case_om = 'M12'

rhoSed = {}
porSed = {}
rdir = '/Users/tanz151/Documents/Projects/TAI_BGC/Drafts/Outputs/PlumIsland/'
xmlfile = rdir + 'optpar_minac.xml'
tree = ET.parse(xmlfile)
root = tree.getroot()        
for key in min_models:
    findstr = "./group/[@id='" + key + 'MOD' + "']/entry"
    for entry in root.findall(findstr):
        if entry.get('id')=='rhoSed':
            rhoSed[key] = float(entry.get('value'))
        elif entry.get('id')=='porSed':
            porSed[key] = float(entry.get('value'))
 
rhoOM = {}           
xmlfile = rdir + 'optpar_omac.xml'
tree = ET.parse(xmlfile)
root = tree.getroot()
for key in om_models:
    findstr = "./group/[@id='" + key + 'MOD' + "']/entry"
    for entry in root.findall(findstr):
        if entry.get('id')=='rhoOM':
            rhoOM[key] = float(entry.get('value'))

# site elevation
z_LAC = 1.1     # Spartina alterniflora-dominated high salt marsh
z_LPC = 1.4     # Spartina patens-dominated high salt marsh
z_MRS = 0.89    # tall Spartina alterniflora low salt marsh
z_MAR = {'A': [0,0.5], 'B': [0.5,1.0], 'C': [1.0,1.5]}

x_LAC = 8.8
x_LPC = 40.9
x_MRS = 2.0

# read simulation
filename = rdir + 'maces_ecogeom_2017-01-01_2019-01-01_' + case_min + '%' + \
    case_om + '_4097.nc'
try:
    nc = Dataset(filename, 'r')
    x = np.array(nc.variables['x'][:])
    zh = np.array(nc.variables['zh'][0,:])
    pft = np.array(nc.variables['pft'][0,:], dtype=np.int8)
finally:
    nc.close()
indx0 = np.argmin(np.abs(zh))
x = x - x[indx0]
nx = len(x)

om_accr_sim = {}
min_accr_sim = {}
tot_accr_sim = {}
for model in min_models:
    for omodel in om_models:
        filename = rdir + 'maces_ecogeom_2017-01-01_2019-01-01_' + model + \
            '%' + omodel + '_4097.nc'
        try:
            nc = Dataset(filename,'r')
            Esed = 8.64e7*np.sum(np.array(nc.variables['Esed'][:]),axis=0)  # kg/m2/yr
            Dsed = 8.64e7*np.sum(np.array(nc.variables['Dsed'][:]),axis=0)  # kg/m2/yr
            om_accr = 8.64e7*np.sum(np.array(nc.variables['DepOM'][:]),axis=0)  # g/m2/yr
        finally:
            nc.close()
        minac_sim_LAC = (np.interp(x_LAC, x, Dsed) - np.interp(x_LAC, x, Esed)) / \
            rhoSed[model] / (1.0-porSed[model]) # mm/yr
        minac_sim_LPC = (np.interp(x_LPC, x, Dsed) - np.interp(x_LPC, x, Esed)) / \
            rhoSed[model] / (1.0-porSed[model]) # mm/yr
        minac_sim_MRS = (np.interp(x_MRS, x, Dsed) - np.interp(x_LAC, x, Esed)) / \
            rhoSed[model] / (1.0-porSed[model]) # mm/yr
        key = model + '%' + omodel
        min_accr_sim[key] = (Dsed - Esed) / rhoSed[model] / (1.0-porSed[model]) # mm/yr
        if omodel==case_om:
            print('MINAC MODEL: ', model, ', ', minac_sim_MRS, minac_sim_LAC, minac_sim_LPC)
        om_accr_sim[key] = om_accr
        tot_accr_sim[key] = min_accr_sim[key] + om_accr_sim[key] / 0.44 / \
            rhoOM[omodel] / (1.0-porSed[model]) # mm/yr
            
# from LTE-MP-LAC-elevationmeans_1.xls
minac_mean_obs_LAC = 5.3    # mm/yr
minac_std_obs_LAC = 0.1     # mm/yr

# from LTE-MP_LPC-elevationmeans_1.xls
minac_mean_obs_LPC = 2.3    # mm/yr
minac_std_obs_LPC = 0.1     # mm/yr

# from Wilson et al. (2014)
minac_mean_obs_MRS = 6.9    # mm/yr
minac_std_obs_MRS = 0.9     # mm/yr

minac_mean_obs = np.array([minac_mean_obs_MRS, minac_mean_obs_LAC, minac_mean_obs_LPC])
minac_std_obs = np.array([minac_std_obs_MRS, minac_std_obs_LAC, minac_std_obs_LPC])
print('MINAC obs: ', minac_mean_obs)
 
Bag_sim_mean_LAC = {}          
# read Law's Point marsh biomass
filename = rdir + 'LTE-MP-LAC-biomassmeans_3.xls'
df = pd.read_excel(filename, sheet_name='LTE-MP-LAC-biomassmeans', header=0, 
                   usecols='A:F')
df.columns = ['Site', 'Year', 'Month', 'Day', 'Bag', 'Bag_std']
Bag_mean_obs_LAC = np.array(df['Bag'])[96:106]      # g/m2
Bag_std_obs_LAC = np.array(df['Bag_std'])[96:106]   # g/m2
Year_LAC = np.array(df['Year'],dtype=np.int32)[96:106]
Month_LAC = np.array(df['Month'],dtype=np.int32)[96:106]
for model in min_models:
    for omodel in om_models:
        filename = rdir + 'maces_ecogeom_2017-01-01_2019-01-01_' + model + \
            '%' + omodel + '_4097.nc'
        try:
            nc = Dataset(filename,'r')
            ntime = len(nc.dimensions['time'])
            Bag = 1e3 * np.array(nc.variables['Bag'][:])    # gC/m2
        finally:
            nc.close()
        Bag_LAC = np.zeros(ntime)
        for jj in range(ntime):
            Bag_LAC[jj] = np.interp(z_LAC, zh, Bag[jj,:])
        Bag_sim_LAC = np.NaN * np.ones_like(Bag_mean_obs_LAC)
        for ii, year in enumerate(Year_LAC):
            month = Month_LAC[ii]
            day0 = (date(year,month,1) - date(2017,1,1)).days
            day1 = (date(year,month+1,1) - date(2017,1,1)).days
            Bag_sim_LAC[ii] = np.mean(Bag_LAC[day0:day1])
        key = model + '%' + omodel
        Bag_sim_mean_LAC[key] = Bag_sim_LAC
        if model==case_min:
            print('OMAC MODEL: ', omodel, ', LAC Bag ', Bag_sim_mean_LAC[key])
        
Bag_sim_mean_MAR = {}
# read 2018/07 biomass in three elevation bands: 0-0.5, 0.5-1.0 and 1.0-1.5
filename = rdir + 'MAR-SO-Biomass-S-alt-2018.xls'
df = pd.read_excel(filename, sheet_name='MAR-SO-Biomass-S-alt-2018', header=0, 
                   usecols='E,F')
df.columns = ['Bag', 'Elevation']
Bag_tmp = np.array(df['Bag'])           # g/m2
Elev_tmp = np.array(df['Elevation'])    # m
Bag_mean_obs_MAR = []
Bag_std_obs_MAR = []
z_MAR_corr = {'A': [-0.1,0.5], 'B': [0.5,1.0], 'C': [1.0,1.5]}
for zstr in ['A','B','C']:
    elev0 = z_MAR[zstr][0]
    elev1 = z_MAR[zstr][1]
    indices = np.logical_and(Elev_tmp>=elev0, Elev_tmp<elev1)
    Bag_mean_obs_MAR.append( np.mean(Bag_tmp[indices]) )
    Bag_std_obs_MAR.append( np.std(Bag_tmp[indices]) )
day0 = (date(2018,7,1) - date(2017,1,1)).days
day1 = (date(2018,8,1) - date(2017,1,1)).days
for model in min_models:
    for omodel in om_models:
        filename = rdir + 'maces_ecogeom_2017-01-01_2019-01-01_' + model + \
            '%' + omodel + '_4097.nc'
        try:
            nc = Dataset(filename,'r')
            ntime = len(nc.dimensions['time'])
            Bag = 1e3 * np.array(nc.variables['Bag'][:])    # gC/m2
        finally:
            nc.close()
        Bag_avg = np.mean(Bag[day0:day1,:],axis=0)
        zh_interp = np.arange(-0.1,1.6,0.1)
        Bag_avg_interp = np.interp(zh_interp, zh, Bag_avg)
        x_interp = np.interp(zh_interp, zh, x)
        nx_interp = len(x_interp)
        dx = np.zeros(nx_interp)
        for jj in range(nx_interp):
            if jj==0:
                dx[jj] = x_interp[jj+1] - x_interp[jj]
            elif jj==nx_interp-1:
                dx[jj] = x_interp[jj] - x_interp[jj-1]
            else:
                dx[jj] = 0.5 * (x_interp[jj+1] - x_interp[jj-1])
        key = model + '%' + omodel
        Bag_sim_mean_MAR[key] = []
        for zstr in ['A','B','C']:
            elev0 = z_MAR[zstr][0]
            elev1 = z_MAR[zstr][1]
            indices = np.logical_and(zh_interp>elev0, zh_interp<=elev1)
            Bag_mean = np.sum(Bag_avg_interp[indices] * dx[indices]) / np.sum(dx[indices])
            Bag_sim_mean_MAR[key].append(Bag_mean)
        if model==case_min:
            print('OMAC MODEL: ', omodel, ', MAR Bag ', Bag_sim_mean_MAR[key])
    
# plot
plt.clf()
fig = plt.figure(figsize=(8,10.5))

gs = gridspec.GridSpec(nrows=3, ncols=2, height_ratios=[1,1,1])

plt.style.use('default')

colors = ['#7b85d4', '#f37738', '#83c995', '#d7369e', '#c4c9d8', '#859795',
          '#e9d043', '#ad5b50', '#e377c2']
#colors = ["#aee39a", "#643176", "#4be32e", "#e72fc2", "#518413", "#7540fc", 
#          "#b3e61c"]
linestyles = ['-', '--', '-.', ':', '-', '--', '-.']

# biomass vs elevation 
ax = fig.add_subplot(gs[0,0])
xpos = np.arange(1,4)
hbar = ax.bar(xpos, Bag_mean_obs_MAR, yerr=Bag_std_obs_MAR, align='center', 
              width=0.8, color='#d8dcd6', ecolor='black', capstyle='butt', 
              capsize=2, alpha=1.0)
handles = []
for indx, omodel in enumerate(om_models):
    key = case_min + '%' + omodel
    h, = ax.plot(xpos, Bag_sim_mean_MAR[key], color=colors[indx], marker='o',
                 ms=5, mfc=colors[indx], mec=colors[indx], 
                 linestyle=linestyles[indx], linewidth=2, alpha=1)
    handles.append(h)
legend = ax.legend(handles, om_models, numpoints=1, loc='upper right', 
                   prop={'family':'Times New Roman', 'size':'large', 'weight': 'bold'}, 
                   framealpha=0.0)
ax.set_xlim(0, 4)
ax.set_ylim(0, 1800)
ax.xaxis.set_ticks(np.arange(1,4,1))
ax.yaxis.set_ticks(np.linspace(0,1800,7))
ax.set_xticklabels(['0–0.5 m','0.5–1 m','1–1.5 m'])
ax.set_xlabel('Elevation', fontsize=12, fontname='Times New Roman', 
              color='black', fontweight='bold')
ylabel = 'Aboveground biomass ($\mathregular{g}$ $\mathregular{{m}^{-2}}$)'
ax.set_ylabel(ylabel, fontsize=12, fontname='Times New Roman', color='black', 
              fontweight='bold')
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(12) for label in labels]
[label.set_color('black') for label in labels]
[label.set_fontweight('bold') for label in labels]
ax.text(0.06, 0.9, 'a', transform=ax.transAxes, fontsize=16,
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
ax.set_xlim(0, 200)
ax.set_ylim(50, 250)
ax.xaxis.set_ticks(np.linspace(0,200,5))
ax.yaxis.set_ticks(np.linspace(50,250,5))
ax.xaxis.set_minor_locator(AutoMinorLocator(5))
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

# biomass seasonality
ax = fig.add_subplot(gs[1,:])
xpos = np.arange(1,11)
hbar = ax.bar(xpos, Bag_mean_obs_LAC, yerr=Bag_std_obs_LAC, align='center', 
              width=0.8, color='#d8dcd6', ecolor='black', capstyle='butt', 
              capsize=2, alpha=1.0)
handles = []
for indx, omodel in enumerate(om_models):
    key = case_min + '%' + omodel
    h, = ax.plot(xpos, Bag_sim_mean_LAC[key], color=colors[indx], marker='o', 
                 ms=5, mfc=colors[indx], mec=colors[indx], 
                 linestyle=linestyles[indx], linewidth=2, alpha=1)
    handles.append(h)
legend = ax.legend(handles, om_models, numpoints=1, loc='upper center',
                   prop={'family':'Times New Roman', 'size':'large', 'weight': 'bold'},
                   framealpha=0.0)
ax.set_xlim(0, 11)
ax.set_ylim(0, 1000)
ax.xaxis.set_ticks(np.arange(1,11,1))
ax.yaxis.set_ticks(np.linspace(0,1000,6))
ax.set_xticklabels(['5/2017','6/2017','7/2017','8/2017','10/2017','5/2018',
                    '6/2018','7/2018','8/2018','10/2018'])
ax.set_xlabel('Time', fontsize=12, fontname='Times New Roman', color='black', 
              fontweight='bold')
ylabel = 'Aboveground biomass ($\mathregular{g}$ $\mathregular{{m}^{-2}}$)'
ax.set_ylabel(ylabel, fontsize=12, fontname='Times New Roman', color='black', 
              fontweight='bold')
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(12) for label in labels]
[label.set_color('black') for label in labels]
[label.set_fontweight('bold') for label in labels]
ax.text(0.03, 0.9, 'c', transform=ax.transAxes, fontsize=16,
        fontname='Times New Roman', fontweight='bold')
ax.tick_params(which='major', direction='in', colors='xkcd:black', length=6, pad=8)
ax.tick_params(which='minor', direction='in', colors='xkcd:black')

# comparison of simulated mineral accretion
ax = fig.add_subplot(gs[2,0])
handles = []
for indx, model in enumerate(min_models):
    key = model + '%' + case_om
    h, = ax.plot(x, min_accr_sim[key], color=colors[indx],
                 linestyle=linestyles[indx], linewidth=2, alpha=1)
    handles.append(h)
hbar = ax.bar([x_MRS,x_LAC,x_LPC], minac_mean_obs, yerr=minac_std_obs, align='center', 
              width=5, color='#d8dcd6', ecolor='black', linewidth=0.5, alpha=1.0)   
#ax.errorbar([x_MRS,x_LAC,x_LPC], minac_mean_obs, xerr=None, yerr=minac_std_obs, 
#            linestyle='None', marker='o', mec='black', mfc='black', ms=5, 
#            ecolor='black', elinewidth=1, capsize=5, alpha=1.0)
legend = ax.legend(handles, min_models, numpoints=1, loc="upper right",
                   prop={'family':'Times New Roman', 'size':'large', 'weight': 'bold'},
                   framealpha=0.0, ncol=1)
ax.set_xlim(-0.8, 200)
ax.set_ylim(0, 12)
ax.xaxis.set_ticks(np.linspace(0,200,5))
ax.yaxis.set_ticks(np.linspace(0,12,5))
ax.xaxis.set_minor_locator(AutoMinorLocator(5))
ax.text(x_MRS+2, minac_mean_obs[0]+0.5, 'MRS', fontsize=12, 
        fontname='Times New Roman', fontweight='bold', alpha=0.8)
ax.text(x_LAC, minac_mean_obs[1]+0.5, 'LAC', fontsize=12, 
        fontname='Times New Roman', fontweight='bold', alpha=0.8)
ax.text(x_LPC, minac_mean_obs[2]+0.5, 'LPC', fontsize=12, 
        fontname='Times New Roman', fontweight='bold', alpha=0.8)
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
ax.text(0.05, 0.9, 'd', transform=ax.transAxes, fontsize=16,
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

ax = fig.add_subplot(gs[2,1])
ax.plot(x, np.percentile(tot_accr_ensemble,50,axis=0), color='black',
        linestyle='-', linewidth=2, alpha=0.9)
ax.fill_between(x, np.percentile(tot_accr_ensemble,25,axis=0),
                np.percentile(tot_accr_ensemble,75,axis=0), alpha=0.3,
                facecolor='black')
ax.set_xlim(0, 200)
ax.set_ylim(0, 12)
ax.xaxis.set_ticks(np.linspace(0,200,5))
ax.yaxis.set_ticks(np.linspace(0,12,5))
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
ax.text(0.05, 0.9, 'e', transform=ax.transAxes, fontsize=16,
        fontname='Times New Roman', fontweight='bold')
ax.tick_params(which='major', direction='in', colors='xkcd:black', length=6, pad=8)
ax.tick_params(which='minor', direction='in', colors='xkcd:black')

axInv = ax.twinx()
axInv.plot(x, np.percentile(fom_accr_ensemble,50,axis=0), color='C3',
           linestyle='-', linewidth=2, alpha=0.9)
axInv.fill_between(x, np.percentile(fom_accr_ensemble,25,axis=0),
                   np.percentile(fom_accr_ensemble,75,axis=0), alpha=0.3,
                   facecolor='C3')
axInv.set_xlim(0, 200)
axInv.set_ylim(0, 100)
axInv.xaxis.set_ticks(np.arange(0,200,5))
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

ax.xaxis.set_ticks(np.linspace(0,200,5))
ax.yaxis.set_ticks(np.linspace(0,12,5))
ax.xaxis.set_minor_locator(AutoMinorLocator(5))

plt.tight_layout()
fig.savefig('F11.png', dpi=300)
#fig.savefig('F11.jpg', dpi=600)
plt.show()
