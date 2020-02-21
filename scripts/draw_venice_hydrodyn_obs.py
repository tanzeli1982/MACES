#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 21 07:21:54 2020

Draw hydrodynamics measurements at Venice Lagoon

@author: Zeli Tan
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

# read data
filename = r'/Users/tanz151/Documents/Projects/TAI_BGC/Data/Hydrodynamics_obs/' + \
    'VeniceLagoon/1BF_OBS.xls'
df = pd.read_excel(filename, sheet_name='1BF', header=None, skiprows=range(3), 
                   usecols='A,B,O,Q')
df.columns = ['Time','Hmo','hw','Turbidity']
Hmo_1BF = 100 * np.array(df['Hmo'])[5334:5526]
hw_1BF = 100 * np.array(df['hw'])[5334:5526]
turb_1BF = np.array(df['Turbidity'])[5334:5526]

filename = r'/Users/tanz151/Documents/Projects/TAI_BGC/Data/Hydrodynamics_obs/' + \
    'VeniceLagoon/2BF_OBS.xls'
df = pd.read_excel(filename, sheet_name='2BF', header=None, skiprows=range(3), 
                   usecols='A,B,O,Q')
df.columns = ['Time','Hmo','hw','Turbidity']
Hmo_2BF = 100 * np.array(df['Hmo'])[5319:5511]
hw_2BF = 100 * np.array(df['hw'])[5319:5511]
turb_2BF = np.array(df['Turbidity'])[5319:5511]

tt = np.arange(192)

# plot
plt.clf()
fig, axes = plt.subplots(3, 1, figsize=(8,9))

plt.style.use('default')

ax = axes[0]
ax.plot(tt, Hmo_2BF, color='C3', marker='.', markersize=5, linestyle=None, 
        alpha=0.9)
ax.set_xlim(0, 193)
ax.set_ylim(0, 100)
ax.xaxis.set_ticks(np.arange(0,193,96))
ax.yaxis.set_ticks(np.linspace(0,100,6))
ax.set_xticklabels(['12/10','12/11','12/12'])
ax.xaxis.set_minor_locator(AutoMinorLocator(24))
ylabel = 'Significant wave height ($\mathregular{cm}$)'
ax.set_ylabel(ylabel, fontsize=14, fontname='Times New Roman', color='black')
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(14) for label in labels]
[label.set_color('black') for label in labels]

ax = axes[1]
ax.plot(tt, hw_2BF, color='C3', marker='.', markersize=5, linestyle=None, 
        linewidth=None, alpha=0.9)
ax.set_xlim(0, 193)
ax.set_ylim(-100, 100)
ax.xaxis.set_ticks(np.arange(0,193,96))
ax.yaxis.set_ticks(np.linspace(-100,100,5))
ax.set_xticklabels(['12/10','12/11','12/12'])
ax.xaxis.set_minor_locator(AutoMinorLocator(24))
ylabel = 'Wave depth ($\mathregular{cm}$)'
ax.set_ylabel(ylabel, fontsize=14, fontname='Times New Roman', color='black')
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(14) for label in labels]
[label.set_color('black') for label in labels]

ax = axes[2]
ax.plot(tt, turb_2BF, color='C3', marker='.', markersize=5, linestyle=None, 
        alpha=0.9)
ax.set_xlim(0, 193)
ax.set_ylim(0, 150)
ax.xaxis.set_ticks(np.arange(0,193,96))
ax.yaxis.set_ticks(np.linspace(0,150,6))
ax.set_xticklabels(['12/10','12/11','12/12'])
ax.xaxis.set_minor_locator(AutoMinorLocator(24))
ylabel = 'Suspended sediment ($\mathregular{mg}$ $\mathregular{{l}^{-1}}$)'
ax.set_ylabel(ylabel, fontsize=14, fontname='Times New Roman', color='black')
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(14) for label in labels]
[label.set_color('black') for label in labels]

plt.tight_layout()
fig.savefig('2BF_hydrodynamics.png', dpi=300)
#fig.savefig('1BF_hydrodynamics.pdf', dpi=600)
plt.show()
