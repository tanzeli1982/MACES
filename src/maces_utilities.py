#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 30 22:03:14 2019

Utility functions for the MACES

@author: Zeli Tan
"""

import numpy as np
from scipy import constants
from scipy.io import netcdf

G = constants.g
Karman = 0.42
Roul = 1028.0
visc = 1e-6     # kinematic viscosity of seawater (m2/s)
        
def construct_tai_platform(diva_segments):
    """Construct the MACES TAI platform.
       DIVA elevations are fixed at -12.5, -8.5, -5.5, -4.5, -3.5, -2.5, -1.5,
       -0.5, 0, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 8.5, 12.5, 16.5 msl.
    Arguments:
        diva_segments['length'] : DIVA segment length (km)
        diva_segments['zhs'] : DIVA segment elevation (m)
        diva_segments['pop'] : DIVA segment population density (people/km^2)
    Returns :
        x_tai   : platform grid coordinate (m)
        zh_tai  : platform grid elevation (m)
        pop_tai : platform grid population density (people/km^2)
    """
    assert len(diva_segments['zhs'])-1 == len(diva_segments['length']), \
        "DIVA segments do not match with elevation nodes"
    Nx = 0
    xRes = 50.0
    nmax = 50
    for length in diva_segments['length']:
        nnode = int( 1e3 * length / xRes )
        Nx = Nx + min( max(nnode,2), nmax )
    Nx = Nx + 1     # the end node
    x_tai = np.zeros(Nx, dtype=np.float64)
    zh_tai = np.zeros(Nx, dtype=np.float64)
    pop_tai = np.zeros(Nx, dtype=np.float64)
    indx = 0
    x0 = 0.0
    zhs = diva_segments['zhs']
    for ii, length in enumerate(diva_segments['length']):
        nnode = min( max( int(1e3*length/xRes), 2 ), nmax )
        x_tai[indx:indx+nnode] = x0 + 1e3*length*np.arange(0,nnode)/nnode
        zh_tai[indx:indx+nnode] = zhs[ii] + (zhs[ii+1]-zhs[ii])*np.arange(0,nnode)/nnode
        pop_tai[indx:indx+nnode] = diva_segments['pop'][ii]
        indx = indx + nnode
        x0 = x0 + 1e3*length
    # the end node
    x_tai[-1] = x0 
    zh_tai[-1] = zhs[-1]
    pop_tai[-1] = diva_segments['pop'][-1]
    return x_tai, zh_tai, pop_tai

def construct_platform_pft(pft_grids, x_tai):
    """Construct the pft on the MACES platform.
    Arguments:
        pft_grids['x'] : pft grid cell coordinate (m)
        pft_grids['pft'] : grid cell pft
        x_tai : coordinate of MACES platform nodes (m)
    Returns : The pft on the MACES platform
    """
    x_arr = pft_grids['x']
    pft_arr = pft_grids['pft']
    pft_tai = np.zeros_like(x_tai, dtype=np.int32)
    for ii, x in enumerate(x_tai):
        if x<x_arr[0]:
            pft_tai[ii] = 1     # tidal flats
        elif x>=x_arr[-1]:
            pft_tai[ii] = pft_arr[-1]
        else:
            indx = np.nonzero(x_arr<=x)[0][-1]
            if x-x_arr[indx]<=x_arr[indx+1]-x:
                pft_tai[ii] = pft_arr[indx]
            else:
                pft_tai[ii] = pft_arr[indx+1]
    return pft_tai
            
def get_refshore_coordinate(diva_segments):
    """Get the coordinate of the shore at msl.
    Arguments:
        diva_segments['length'] : DIVA segment length (km)
        diva_segments['zhs'] : DIVA segment elevation (m)
    Returns : the coordinate of shore at msl (m)
    """
    xlens = diva_segments['length']
    zhs = diva_segments['zhs']
    assert len(zhs)-1 == len(xlens), "DIVA segments do not match with elevation nodes"
    indx = np.nonzero(zhs==0)[0]
    assert len(indx)==1, "Shore node does not exist in DIVA segments"
    return np.sum(xlens[:indx[0]])

def write_outputs(odir, sid, uhydro_out, ecogeom_out):
    """Write model outputs into a nc file.
    Arguments:
        odir : output directory
        sid : site id
        uhydro_out  : hourly hydrodynamic model outputs
        ecogeom_out : daily eco-geomorphology outputs 
    Returns : 
    """
    zh = ecogeom_out['zh']
    nday = np.shape(zh)[0]
    nx = np.shape(zh)[1]
    filename = odir + '{:d}'.format(sid) + '.nc'
    try:
        nc = netcdf.netcdf_file(filename, version=2)
        nc.history = r'Simulated hourly hydrodynamics and daily eco-geomorphology by MACES'
        nc.contact = r'Zeli.Tan@pnnl.gov'
        nc.createDimension('hour', None)
        nc.createDimension('day', nday)
        nc.createDimension('x', nx)
        # create and write hourly variables
        h_var = nc.createVariable('h', 'f4', ('hour','x',))
        h_var.long_name = r'water depth'
        h_var.units = 'm'
        h_var._FillValue = np.float32(1e20)
        h_var[:] = uhydro_out['h']
        U_var = nc.createVariable('U', 'f4', ('hour','x',))
        U_var.long_name = r'tide flow velocity (direction + speed)'
        U_var.units = 'm/s'
        U_var._FillValue = np.float32(1e20)
        U_var[:] = uhydro_out['U']
        Hwav_var = nc.createVariable('Hwav', 'f4', ('hour','x',))
        Hwav_var.long_name = r'significant wave height'
        Hwav_var.units = 'm'
        Hwav_var._FillValue = np.float32(1e20)
        Hwav_var[:] = uhydro_out['Hwav']
        tau_var = nc.createVariable('tau', 'f4', ('hour','x',))
        tau_var.long_name = r'bottom shear stress'
        tau_var.units = 'Pa'
        tau_var._FillValue = np.float32(1e20)
        tau_var[:] = uhydro_out['tau']
        Css_var = nc.createVariable('Css', 'f4', ('hour','x',))
        Css_var.long_name = r'suspended sediment concentration'
        Css_var.units = 'kg/m3'
        Css_var._FillValue = np.float32(1e20)
        Css_var[:] = uhydro_out['Css']
        Cj_var = nc.createVariable('Cj', 'f4', ('hour','x',))
        Cj_var.long_name = r'water salinity'
        Cj_var.units = 'PSU'
        Cj_var._FillValue = np.float32(1e20)
        Cj_var[:] = uhydro_out['Cj']   
        # create and write daily variables
        zh_var = nc.createVariable('zh', 'f4', ('day','x',))
        zh_var.long_name = r'platform surface elevation'
        zh_var.units = 'msl'
        zh_var._FillValue = np.float32(1e20)
        zh_var[:] = ecogeom_out['zh']
        Esed_var = nc.createVariable('Esed', 'f4', ('day','x',))
        Esed_var.long_name = r'sediment erosion rate'
        Esed_var.units = 'kg/m2/s'
        Esed_var._FillValue = np.float32(1e20)
        Esed_var[:] = ecogeom_out['Esed']
        Dsed_var = nc.createVariable('Dsed', 'f4', ('day','x',))
        Dsed_var.long_name = r'sediment deposition rate'
        Dsed_var.units = 'kg/m2/s'
        Dsed_var._FillValue = np.float32(1e20)
        Dsed_var[:] = ecogeom_out['Dsed']
        DepOM_var = nc.createVariable('DepOM', 'f4', ('day','x',))
        DepOM_var.long_name = r'Organic matter accretion rate'
        DepOM_var.units = 'kg/m2/s'
        DepOM_var._FillValue = np.float32(1e20)
        DepOM_var[:] = ecogeom_out['DepOM']
        Bag_var = nc.createVariable('Bag', 'f4', ('day','x',))
        Bag_var.long_name = r'platform aboveground biomass'
        Bag_var.units = 'kg/m2'
        Bag_var._FillValue = np.float32(1e20)
        Bag_var[:] = ecogeom_out['Bag']
    finally:
        nc.close()
    
    