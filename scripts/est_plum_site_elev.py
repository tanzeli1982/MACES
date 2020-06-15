#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 15 11:13:10 2020

Estimate Plum Island site elevation

@author: Zeli Tan
"""

import numpy as np
import pandas as pd
from netCDF4 import Dataset
from datetime import date

day0 = (date(2017,7,19) - date(2012,1,1)).days
day1 = (date(2017,7,23) - date(2012,1,1)).days
# read forcing water level
filename = '/Users/tanz151/Documents/Projects/TAI_BGC/Data/Hydrodynamics_obs/' + \
    'PlumIsland/force_h.nc'
try:
    nc = Dataset(filename, 'r')
    h_var = 100 * np.array(nc.variables['h'][96*day0:96*day1,0])
finally:
    nc.close()
    
# read site water depth at gauge
filename = r'/Users/tanz151/Documents/Projects/TAI_BGC/Data/Hydrodynamics_obs/' + \
    r'PlumIsland/NelsonIsland/MAR-RO-Wtable-Nel-2017.xls'
df = pd.read_excel(filename, sheet_name='MAR-RO-Wtable-Nel-2017', header=0, 
                   usecols='A,B,C,K')
df.columns = ['Date', 'Time', 'Gauge', 'N204']
h_nelson_gauge = 100 * np.array(df['Gauge'])[16895:18047]
h_nelson_n204 = 100 * np.array(df['N204'])[16895:18047] - 198
h_nelson_gauge[h_nelson_gauge<0] = 0
h_nelson_n204[h_nelson_n204<0] = 0

filename = r'/Users/tanz151/Documents/Projects/TAI_BGC/Data/Hydrodynamics_obs/' + \
    r'PlumIsland/ShadCreek/MAR-RO-Wtable-Shad-2017.xls'
df = pd.read_excel(filename, sheet_name='MAR-RO-Wtable-Shad-2017', header=0, 
                   usecols='A,B,C,G')
df.columns = ['Date', 'Time', 'Gauge', 'S102']
h_shad_gauge = 100 * np.array(df['Gauge'])[17086:18238]
h_shad_s102 = 100 * np.array(df['S102'])[17086:18238] - 98.5
h_shad_gauge[h_shad_gauge<0] = 0
h_shad_s102[h_shad_s102<0] = 0