#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Mar  9 15:19:45 2019

Compare the simulated CH4 and CO2 emissions at two SERC sites

@author: Zeli Tan
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.io import netcdf
from matplotlib.ticker import AutoMinorLocator

# read ten year simulations
years = np.arange(2001,2011)
nyear = np.size(years)
#days = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
#rdir = '/pic/scratch/tanz151/csmruns/output/serc_1pt.1.CLM_USRDAT.' + \
#    'I20TRGSWCNPECACNTBC.constance.transient/run/serc_1pt.1.CLM_USRDAT.' + \
#    'I20TRGSWCNPECACNTBC.constance.transient.clm2.h0.' 
#fch4_sat = np.zeros((nyear*12,2))       # saturated CH4 flux (mg/m2/d)
#fch4_unsat = np.zeros((nyear*12,2))     # unsaturated CH4 flux (mg/m2/d)
#fch4 = np.zeros((nyear*12,2))           # total CH4 flux (mg/m2/d)
#fco2 = np.zeros((nyear*12,2))           # total CO2 flux (g/m2/d)
#tair = np.zeros((nyear*12,2))           # air temperature (celsius)
#tgnd = np.zeros((nyear*12,2))           # ground temperature (celsius)
#prcp = np.zeros((nyear*12,2))           # precipitation (mm/d)
#for ii, year in enumerate(years):
#    for jj in range(12):
#        for kk in range(days[jj]):
#            filename = rdir + '{:d}'.format(year) + '-' + \
#                '{:02d}'.format(jj+1) + '-' + '{:02d}'.format(kk+1) + '-00000.nc'
#            try:
#                nc = netcdf.netcdf_file(filename,'r')
#                fch4_plt_sat = np.array(nc.variables['CH4_SURF_AERE_SAT'][:])   # mol/m2/s
#                fch4_plt_unsat = np.array(nc.variables['CH4_SURF_AERE_UNSAT'][:])
#                fch4_diff_sat = np.array(nc.variables['CH4_SURF_DIFF_SAT'][:])
#                fch4_diff_unsat = np.array(nc.variables['CH4_SURF_DIFF_UNSAT'][:])
#                fch4_ebul_sat = np.array(nc.variables['CH4_SURF_EBUL_SAT'][:])
#                fch4_ebul_unsat = np.array(nc.variables['CH4_SURF_EBUL_UNSAT'][:])
#                fsat = np.array(nc.variables['FINUNDATED'][:])  # fraction
#                nee = np.array(nc.variables['NEE'][:])  # gC/m2/s
#                tair_var = np.array(nc.variables['TBOT'][:])   # K
#                tgnd_var = np.array(nc.variables['TG'][:])     # K
#                rain_var = np.array(nc.variables['RAIN'][:])   # mm/s
#            finally:
#                nc.close()
#            fch4_unsat_day = 16e3 * 8.64e4 * (fch4_plt_unsat + fch4_diff_unsat + fch4_ebul_unsat)
#            fch4_sat_day = 16e3 * 8.64e4 * (fch4_plt_sat + fch4_diff_sat + fch4_ebul_sat)
#            fch4_unsat[12*ii+jj,:] = fch4_unsat[12*ii+jj,:] + fch4_unsat_day / float(days[jj])
#            fch4_sat[12*ii+jj,:] = fch4_sat[12*ii+jj,1] + fch4_sat_day / float(days[jj])
#            fco2[12*ii+jj,:] = fco2[12*ii+jj,:] + 8.64e4 * nee / float(days[jj])
#            fch4[12*ii+jj,:] = fch4[12*ii+jj,:] + (fch4_sat_day*fsat + \
#                fch4_unsat_day*(1-fsat)) / float(days[jj])
#            tair[12*ii+jj,:] = tair[12*ii+jj,:] + (tair_var - 273.15) / float(days[jj])
#            tgnd[12*ii+jj,:] = tgnd[12*ii+jj,:] + (tgnd_var - 273.15) / float(days[jj])
#            prcp[12*ii+jj,:] = prcp[12*ii+jj,:] + 8.64e4 * rain_var / float(days[jj])
#            
## save as a netcdf file
#filename = '/people/tanz151/scripts/TAI/serc_fco2_fch4.nc'
#f = netcdf.netcdf_file(filename, 'w')
#try:
#    f.history = r'monthly averaged daily CO2 and CH4 fluxes'
#    f.createDimension('time', None)     # unlimited dimension
#    f.createDimension('site', 2)
#    fco2_var = f.createVariable('fco2', 'f4', ('time','site',))
#    fco2_var.long_name = r'mean daily CO2 flux'
#    fco2_var.units = 'gC/m^2/d'
#    fco2_var._FillValue = np.float32(1e36)
#    fco2_var[:] = fco2
#    fch4_var = f.createVariable('fch4', 'f4', ('time','site',))
#    fch4_var.long_name = r'mean daily CH4 flux'
#    fch4_var.units = 'mg/m^2/d'
#    fch4_var._FillValue = np.float32(1e36)
#    fch4_var[:] = fch4
#    fch4_sat_var = f.createVariable('fch4_sat', 'f4', ('time','site',))
#    fch4_sat_var.long_name = r'mean daily CH4 flux from saturated domain'
#    fch4_sat_var.units = 'mg/m^2/d'
#    fch4_sat_var._FillValue = np.float32(1e36)
#    fch4_sat_var[:] = fch4_sat
#    fch4_unsat_var = f.createVariable('fch4_unsat', 'f4', ('time','site',))
#    fch4_unsat_var.long_name = r'mean daily CH4 flux from unsaturated domain'
#    fch4_unsat_var.units = 'mg/m^2/d'
#    fch4_unsat_var._FillValue = np.float32(1e36)
#    fch4_unsat_var[:] = fch4_unsat  
#    tair_var = f.createVariable('tair', 'f4', ('time','site',))
#    tair_var.long_name = r'mean daily air temperature'
#    tair_var.units = 'celsius'
#    tair_var._FillValue = np.float32(1e36)
#    tair_var[:] = tair
#    tgnd_var = f.createVariable('tgnd', 'f4', ('time','site',))
#    tgnd_var.long_name = r'mean daily ground temperature'
#    tgnd_var.units = 'celsius'
#    tgnd_var._FillValue = np.float32(1e36)
#    tgnd_var[:] = tgnd
#    prcp_var = f.createVariable('prcp', 'f4', ('time','site',))
#    prcp_var.long_name = r'mean daily rainfall'
#    prcp_var.units = 'mm/d'
#    prcp_var._FillValue = np.float32(1e36)
#    prcp_var[:] = prcp
#finally:
#    f.close()

try:
    nc = netcdf.netcdf_file('serc_fco2_fch4.nc','r')
    fco2 = np.array(nc.variables['fco2'][:])
    fch4_sat = np.array(nc.variables['fch4_sat'][:])
    fch4_unsat = np.array(nc.variables['fch4_unsat'][:])
    fch4 = np.array(nc.variables['fch4'][:])
finally:
    nc.close()

# plot
plt.clf()
fig = plt.figure(figsize=(8.5,6))

ax1 = plt.subplot2grid((2, 4), (0, 1), colspan=2)
ax2 = plt.subplot2grid((2, 4), (1, 0), colspan=2)
ax3 = plt.subplot2grid((2, 4), (1, 2), colspan=2)

axes = [ax1, ax2, ax3]

colors = ['#1f77b4', '#d62728']

plt.style.use('default')
tt = np.arange(12*nyear)

ax = axes[0]
ax.plot(tt, fco2[:,0], color=colors[0], linestyle='-', linewidth=2, alpha=0.9)
ax.plot(tt, fco2[:,1], color=colors[1], linestyle='-', linewidth=2, alpha=0.9)
ax.legend(['Upland forest habitat','Marsh habitat'], numpoints=1,
          loc=2, ncol=1, framealpha=0.0, 
          prop={'family':'Times New Roman', 'size':'large', 'weight':'bold'})
ax.set_xlim(0, nyear*12)
ax.set_ylim(-7.5, 12.5)
ax.xaxis.set_ticks(np.arange(0,nyear*12+1,24))
ax.yaxis.set_ticks(np.linspace(-5,10,4))
ax.set_xticklabels(['2001','2003','2005','2007','2009','2011'])
ax.xaxis.set_minor_locator(AutoMinorLocator(2))
ylabel = '$\mathregular{{CO}_{2}}$ flux ($\mathregular{gC}$ ' + \
    '$\mathregular{{m}^{-2}}$ $\mathregular{{d}^{-1}}$)'
ax.set_ylabel(ylabel, fontsize=14, fontname='Times New Roman', color='black')
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(12) for label in labels]
[label.set_color('black') for label in labels]
ax.text(0.9, 0.9, 'A', transform=ax.transAxes, fontsize=14,
        fontname='Times New Roman', fontweight='bold')

ax = axes[1]
ax.plot(tt, fch4[:,0], color=colors[0], linestyle='-', linewidth=2, alpha=0.9)
ax.plot(tt, fch4[:,1], color=colors[1], linestyle='-', linewidth=2, alpha=0.9)
ax.set_xlim(0, nyear*12)
ax.set_ylim(-100, 2100)
ax.xaxis.set_ticks(np.arange(0,nyear*12+1,24))
ax.yaxis.set_ticks(np.linspace(0,2000,5))
ax.set_xticklabels(['2001','2003','2005','2007','2009','2011'])
ax.xaxis.set_minor_locator(AutoMinorLocator(2))
ylabel = '$\mathregular{{CH}_{4}}$ flux ($\mathregular{mg}$ ' + \
    '$\mathregular{{m}^{-2}}$ $\mathregular{{d}^{-1}}$)'
ax.set_ylabel(ylabel, fontsize=14, fontname='Times New Roman', color='black')
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(12) for label in labels]
[label.set_color('black') for label in labels]
ax.text(0.05, 0.9, 'B', transform=ax.transAxes, fontsize=14,
        fontname='Times New Roman', fontweight='bold')

ax = axes[2]
ax.plot(tt, fch4_unsat[:,0], color=colors[0], linestyle='-', linewidth=2, alpha=0.9)
ax.plot(tt, fch4_sat[:,1], color=colors[1], linestyle='-', linewidth=2, alpha=0.9)
ax.set_xlim(0, nyear*12)
ax.set_ylim(-100, 2100)
ax.xaxis.set_ticks(np.arange(0,nyear*12+1,24))
ax.yaxis.set_ticks(np.linspace(0,2000,5))
ax.set_xticklabels(['2001','2003','2005','2007','2009','2011'])
ax.xaxis.set_minor_locator(AutoMinorLocator(2))
ylabel = '$\mathregular{{CH}_{4}}$ flux ($\mathregular{mg}$ ' + \
    '$\mathregular{{m}^{-2}}$ $\mathregular{{d}^{-1}}$)'
ax.set_ylabel(ylabel, fontsize=14, fontname='Times New Roman', color='black')
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
[label.set_fontsize(12) for label in labels]
[label.set_color('black') for label in labels]
ax.text(0.05, 0.9, 'C', transform=ax.transAxes, fontsize=14,
        fontname='Times New Roman', fontweight='bold')

plt.tight_layout()
fig.savefig('serc_fco2_fch4.png', dpi=300)
#fig.savefig('serc_fco2_fch4.pdf', dpi=600)
plt.show()