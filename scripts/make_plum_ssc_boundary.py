#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  3 13:48:31 2023

Produce SSC boundary for Plum Island site

@author: Zeli Tan
"""

import numpy as np
import pandas as pd
from netCDF4 import Dataset
from datetime import date

rdir = '/Users/tanz151/Documents/Projects/TAI_BGC/Drafts/Outputs/PlumIsland/'

# read forcing data
filename = rdir + 'force_h.nc'
try:
    nc = Dataset(filename,'r')
    h_var = 100 * np.array(nc.variables['h'][:,0]) # cm
finally:
    nc.close()
nt_h = np.size(h_var)
tt_h = np.arange(nt_h) / 4.0

filename = rdir + 'force_U10.nc'
try:
    nc = Dataset(filename,'r')
    U10_var = 100 * np.array(nc.variables['U10'][:,0]) # m/s
finally:
    nc.close()

# read SSC data
filename = rdir + 'EST-RO-TC-WE-TurbiditySuspSed.xlsx'
df = pd.read_excel(filename, sheet_name='EST-RO-TC-WE-TurbiditySuspSed', 
                   header=0, usecols='A,B,C,E,F')
df.columns = ['Date','Time','Site','SSC','Depth']
sed_c = np.array(df['SSC'])[20350:20734]
sed_m = np.array(df['SSC'])[87632:88016]
sed_obs_c = np.NaN * np.ones(24*(day1-day0))
sed_obs_m = np.NaN * np.ones(24*(day1-day0))
nt_obs = np.size(sed_obs_c)
for ii in range(nt_obs):
    sed_obs_c[ii] = np.nanmean(sed_c[4*ii:4*(ii+1)])
    sed_obs_m[ii] = np.nanmean(sed_m[4*ii:4*(ii+1)])

nt_obs = np.size(sed_obs_c)
tt_obs = np.arange(nt_obs)