#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 21 08:59:21 2020

Draw simulated bottom shear stress at 1BF and 2BF.

@author: Zeli Tan
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

# read data
filename = r'/Users/tanz151/Documents/Projects/TAI_BGC/Data/Hydrodynamics_obs/' + \
    r'VeniceLagoon/2-4Apr03-TauMax1BF&2BF.txt'
tau_1BF = []
tau_2BF = []
try:
    f = open(filename, 'r')
    f.readline()    # skip header
    for line in f:
        line = line.strip()
        columns = line.split()
        tau_1BF.append(float(columns[1]))
        tau_2BF.append(float(columns[2]))
finally:
    f.close()
    
tau_1BF_np = np.array(tau_1BF)
tau_2BF_np = np.array(tau_2BF)
tt = np.arange(145)

# plot
plt.clf()
fig, axes = plt.subplots(2, 1, figsize=(5,7))

plt.style.use('default')

ax = axes[0]
ax.plot(tt, tau_1BF_np, color='black', linestyle='-', linewidth=2, alpha=0.9)
ax.set_xlim(0, 146)
ax.set_ylim(0, 1)
ax.xaxis.set_ticks(np.arange(0,146,48))
ax.yaxis.set_ticks(np.linspace(0,1,6))
ax.set_xticklabels(['4/2','4/3','4/4','4/5'])
ax.xaxis.set_minor_locator(AutoMinorLocator(24))
ylabel = 'Bottom shear stress ($\mathregular{Pa}$)'
ax.set_ylabel(ylabel, fontsize=12, fontname='Times New Roman', color='black')
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(12) for label in labels]
[label.set_color('black') for label in labels]

ax = axes[1]
ax.plot(tt, tau_2BF_np, color='black', linestyle='-', linewidth=2, alpha=0.9)
ax.set_xlim(0, 146)
ax.set_ylim(0, 1)
ax.xaxis.set_ticks(np.arange(0,146,48))
ax.yaxis.set_ticks(np.linspace(0,1,6))
ax.set_xticklabels(['4/2','4/3','4/4','4/5'])
ax.xaxis.set_minor_locator(AutoMinorLocator(24))
ylabel = 'Bottom shear stress ($\mathregular{Pa}$)'
ax.set_ylabel(ylabel, fontsize=12, fontname='Times New Roman', color='black')
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(12) for label in labels]
[label.set_color('black') for label in labels]

plt.tight_layout()
fig.savefig('venice_tau.png', dpi=300)
#fig.savefig('venice_tau.pdf', dpi=600)
plt.show()