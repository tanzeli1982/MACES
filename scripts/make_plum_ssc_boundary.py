#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  3 13:48:31 2023

Produce SSC boundary for Plum Island site

@author: Zeli Tan
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from netCDF4 import Dataset
from datetime import date
from statsmodels.formula.api import ols

rdir = '/Users/tanz151/Documents/Projects/TAI_BGC/Drafts/Outputs/PlumIsland/'

# read forcing data
filename = rdir + 'force_h.nc'
try:
    nc = Dataset(filename,'r')
    h_var = np.array(nc.variables['h'][:,0]) # m
finally:
    nc.close()
nt_h = np.size(h_var)
tt_h = range(nt_h)

filename = rdir + 'force_U10.nc'
try:
    nc = Dataset(filename,'r')
    U10_var = np.array(nc.variables['U10'][:,0]) # m/s
finally:
    nc.close()

# read SSC data
filename = rdir + 'EST-RO-TC-WE-TurbiditySuspSed.xlsx'
df = pd.read_excel(filename, sheet_name='EST-RO-TC-WE-TurbiditySuspSed', 
                   header=0, usecols='A,B,C,E,F')
df.columns = ['Date','Time','Site','SSC','Depth']
sed_obs = np.array(df['SSC'])
sed_site = list(df['Site'])
sed_date = list(df['Date'])
sed_time = list(df['Time'])

date0 = date(2012, 1, 1)
sed_channel_obs = []
h_var_sel = []
U10_var_sel = []
for ii, name in enumerate(sed_site):
    if name=='EST-TC-West-LP-channel' and np.isfinite(sed_obs[ii]):
        tt = int( 4*( 24*(sed_date[ii].date() - date0).days + \
                 sed_time[ii].hour + sed_time[ii].minute/60.0) )
        indx = tt_h.index(tt)
        sed_channel_obs.append(sed_obs[ii])
        h_var_sel.append(h_var[indx])
        U10_var_sel.append(U10_var[indx])
        
sed_channel_obs = np.array(sed_channel_obs)
h_var_sel = np.array(h_var_sel)
U10_var_sel = np.array(U10_var_sel)

# linear regression
plt.clf()
fig = plt.figure(figsize=(4,8))

plt.style.use('default')
gs = gridspec.GridSpec(nrows=2, ncols=1)
   
# h 
ax = fig.add_subplot(gs[0,0])

data = {'h': h_var_sel, 'ssc': sed_channel_obs}
df = pd.DataFrame(data = data)
fitted = ols('ssc ~ h', data = df).fit()
print( fitted.params[0], fitted.params[1], fitted.rsquared, fitted.f_pvalue )
h_fitted = np.linspace(np.min(h_var_sel), np.max(h_var_sel), 100)
ssc_fitted = fitted.params[0] + fitted.params[1]*h_fitted

ssc_ts = np.maximum( 2.0*(fitted.params[0] + fitted.params[1]*h_var), 0.0)

ax.scatter(h_var_sel, sed_channel_obs, s=10, c='black', marker='.')
ax.plot(h_fitted, ssc_fitted, color='darkgreen', ls='-', lw=2, alpha=0.8)
ax.plot([h_fitted[0],h_fitted[-1]], [ssc_fitted[0],ssc_fitted[-1]], 
        color='grey', ls='--', lw=2, alpha=0.8)
ax.set_xlim(h_fitted[0], h_fitted[-1])
ax.set_ylim(np.min(sed_channel_obs), np.max(sed_channel_obs))
xlabel = 'h ($\mathregular{m}$)'
ylabel = 'SSC ($\mathregular{mg}$ $\mathregular{{l}^{-1}}$)'
ax.set_xlabel(xlabel, color='black', fontsize=14, fontname='Times New Roman', 
              fontweight='bold')
ax.set_ylabel(ylabel, color='black', fontsize=14, fontname='Times New Roman', 
              fontweight='bold')
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(12) for label in labels]
[label.set_color('black') for label in labels]
[label.set_fontweight('bold') for label in labels]

# U10
ax = fig.add_subplot(gs[1,0])

data = {'U10': U10_var_sel, 'ssc': sed_channel_obs}
df = pd.DataFrame(data = data)
fitted = ols('ssc ~ U10', data = df).fit()
print( fitted.params[1], fitted.rsquared, fitted.f_pvalue )
U10_fitted = np.linspace(np.min(U10_var_sel), np.max(U10_var_sel), 100)
ssc_fitted = fitted.params[0] + fitted.params[1]*U10_fitted

ax.scatter(U10_var_sel, sed_channel_obs, s=10, c='black', marker='.')
ax.plot(U10_fitted, ssc_fitted, color='darkgreen', ls='-', lw=2, alpha=0.8)
ax.plot([U10_fitted[0],U10_fitted[-1]], [ssc_fitted[0],ssc_fitted[-1]], 
        color='grey', ls='--', lw=2, alpha=0.8)
ax.set_xlim(U10_fitted[0], U10_fitted[-1])
ax.set_ylim(np.min(sed_channel_obs), np.max(sed_channel_obs))
xlabel = 'U10 ($\mathregular{m}$ $\mathregular{{s}^{-1}}$)'
ylabel = 'SSC ($\mathregular{mg}$ $\mathregular{{l}^{-1}}$)'
ax.set_xlabel(xlabel, color='black', fontsize=14, fontname='Times New Roman', 
              fontweight='bold')
ax.set_ylabel(ylabel, color='black', fontsize=14, fontname='Times New Roman', 
              fontweight='bold')
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(12) for label in labels]
[label.set_color('black') for label in labels]
[label.set_fontweight('bold') for label in labels]

plt.tight_layout()
plt.show()

# produce sediment boundary file
print('SSC: ', np.min(ssc_ts), np.max(ssc_ts))
filename = rdir + 'force_SSC_revised.nc'
try:
    nc = Dataset(filename, 'w', format='NETCDF3_64BIT_OFFSET')
    nc.createDimension('time', None)
    nc.createDimension('site', 1)
    # global variables
    nc.history = 'Plum Island 15-min suspended sediment forcing data'
    # variable
    date_var = nc.createVariable('date', 'i4')
    date_var[:] = 20120101
    ssc_var = nc.createVariable('TSM', 'f8', ('time','site',))
    ssc_var.long_name = 'suspended sediment concentration'
    ssc_var.units = 'mg/l'
    ssc_var[:] = ssc_ts
finally:
    nc.close()
    
#filename = rdir + 'force_h.nc'
#try:
#    nc = Dataset(filename, 'r+')
#    h_ovar = nc.variables['h']
#    h_ovar[:] = h_var
#finally:
#    nc.close()