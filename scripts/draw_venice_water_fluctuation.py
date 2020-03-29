#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 27 13:56:13 2020

Compare the water level fluctuation between 1BF and 2BF

@author: Zeli Tan
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# read data for 1BF
filename = '/Users/tanz151/Documents/Projects/TAI_BGC/Data/Hydrodynamics_obs/' + \
    'VeniceLagoon/WaterLevelClose1BF.xls'
df = pd.read_excel(filename, sheet_name='Valori orari 2002', header=None, 
                   skiprows=range(4), usecols='A:C')
df.columns = ['Date','Hour','hw']
h_obs_1BF = 100 * np.array(df['hw'])[8015:8759]

nt_obs_1BF = np.size(h_obs_1BF)
tt_obs_1BF = np.arange(nt_obs_1BF)

# read data for 2BF
filename = '/Users/tanz151/Documents/Projects/TAI_BGC/Data/Hydrodynamics_obs/' + \
    'VeniceLagoon/2BF_OBS.xls'
df = pd.read_excel(filename, sheet_name='2BF', header=None, skiprows=range(3), 
                   usecols='A,O')
df.columns = ['Time','hw']
h_obs_2BF = 100 * np.array(df['hw'])[4455:7431]

nt_obs_2BF = np.size(h_obs_2BF)
tt_obs_2BF = np.arange(nt_obs_2BF)/4

# plot
plt.clf()
fig, ax = plt.subplots(figsize=(8,6))

plt.style.use('default')

ax.plot(tt_obs_1BF, h_obs_1BF, color='black', linestyle='-', linewidth=2, alpha=0.9)
ax.plot(tt_obs_2BF, h_obs_2BF, color='C3', linestyle='-', linewidth=2, alpha=0.9)
ax.set_xlim(0, nt_obs_1BF-1)
#ax.set_ylim(-100, 100)
ax.xaxis.set_ticks(np.arange(0,nt_obs_1BF,24))
#ax.yaxis.set_ticks(np.linspace(-100,100,5))
labels = ['{:d}'.format(ii+1) for ii in range(31)]
ax.set_xticklabels(labels)
#ax.xaxis.set_minor_locator(AutoMinorLocator(24))
ylabel = 'Water depth ($\mathregular{cm}$)'
ax.set_ylabel(ylabel, fontsize=12, fontname='Times New Roman', color='black')
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(12) for label in labels]
[label.set_color('black') for label in labels]

plt.tight_layout()
fig.savefig('venice_water_fluctuation.png', dpi=300)
#fig.savefig('venice_water_fluctuation.pdf', dpi=600)
plt.show()
