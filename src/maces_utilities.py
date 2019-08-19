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

MAX_OF_STEP = 1800  # maximum simulation time step (s)

def get_date_from_julian(julian):
    """Get Julian day number from date
    Returns : year, month, day
    """
    nf = julian + 1401 + int((int((4*julian+274277)/146097)*3)/4) - 38
    ne = 4*nf + 3
    ng = int(np.mod(ne,1461)/4)
    nh = 5*ng + 2
    day = int(np.mod(nh,153)/5) + 1
    month = np.mod(int(nh/153)+2,12) + 1
    year = int(ne/1461) - 4716 + int((14-month)/12)
    return year, month, day
    
def get_julian_from_date(year, month, day):
    """Get Julian day number from date
    Returns : Julian day number
    """
    a = int((14-month)/12)
    y = year + 4800 - a
    m = month + 12*a - 3
    date = 1e4 * year + 1e2 * month + day
    if date>=15821015:
        julian = day + int((153*m+2)/5) + 365*y + int(y/4) - int(y/100) + \
            int(y/400) - 32045
    else:
        julian = day + int((153*m+2)/5) + 365*y + int(y/4) - 32083
    return julian
        
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
    x_tai = np.zeros(Nx, dtype=np.float64, order='F')
    zh_tai = np.zeros(Nx, dtype=np.float64, order='F')
    pop_tai = np.zeros(Nx, dtype=np.float64, order='F')
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
        pft_grids['x'] : pft grid cell coordinate (km)
        pft_grids['pft'] : grid cell pft
        x_tai : coordinate of MACES platform nodes (m)
    Returns : The pft on the MACES platform
    """
    x_arr = 1e3 * pft_grids['x']    # km -> m
    pft_arr = pft_grids['pft']
    pft_tai = np.zeros_like(x_tai, dtype=np.int8, order='F')
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
            
def get_refshore_coordinate(x, zh):
    """Get the coordinate of the shore at msl.
    Arguments:
        x : longitudinal coordinate (m)
        zh : platform surface elevation (msl)
    Returns : the coordinate of shore at msl (km)
    """
    indx = np.nonzero(zh==0)[0]
    assert len(indx)==1, "Shore node does not exist in DIVA segments"
    return 1e-3 * x[indx[0]]

def get_platform_slope(x, zh):
    """Get the platform slope.
    Arguments:
        x : longitudinal coordinate (m)
        zh : platform surface elevation (msl)
    Returns : the platform cell slope (m/m)
    """
    Nx = np.size(x)
    dzh = np.ones(Nx, dtype=np.float64, order='F')
    for ii in range(Nx):
        if ii>0 and ii<Nx-1:
            dzh[ii] = (zh[ii+1] - zh[ii-1]) / (x[ii+1] - x[ii-1])
        elif ii==0:
            dzh[ii] = (zh[ii+1] - zh[ii]) / (x[ii+1] - x[ii])
        else:
            dzh[ii] = (zh[ii] - zh[ii-1]) / (x[ii] - x[ii-1])
    return dzh

def estimate_Hwav_seaward(U10, h0):
    """Estimate the significant wave height at the seaward side.
    Arguments:
        U10 : 10-m wind speed (m/s)
        h0 : water depth at the seaward side (m)
    Returns : significant wave height (m)
    """
    Hwav_ocean = 0.27 * U10**2 / G
    Hwav0 = np.sinh(Karman*h0) / np.sinh(Karman*30.0) * Hwav_ocean
    return Hwav0

def write_hydro_outputs(odir, sid, uhydro_out):
    """Write model hydrodynamics outputs into a nc file.
    Arguments:
        odir : output directory
        sid : site id
        uhydro_out  : hourly hydrodynamic model outputs
    Returns : 
    """
    nx = np.shape(uhydro_out['h'])[1]
    filename = odir + '{:d}'.format(sid) + '_hydro.nc'
    try:
        nc = netcdf.netcdf_file(filename, 'w', version=2)
        nc.history = r'MACES simulated hourly hydrodynamics'
        nc.contact = r'Zeli.Tan@pnnl.gov'
        nc.createDimension('time', None)
        nc.createDimension('x', nx)
        # create and write variables
        h_var = nc.createVariable('h', 'f4', ('time','x',))
        h_var.long_name = r'water depth'
        h_var.units = 'm'
        h_var._FillValue = np.float32(1e20)
        h_var[:] = uhydro_out['h']
        U_var = nc.createVariable('U', 'f4', ('time','x',))
        U_var.long_name = r'tide flow velocity (direction * speed)'
        U_var.units = 'm/s'
        U_var._FillValue = np.float32(1e20)
        U_var[:] = uhydro_out['U']
        Hwav_var = nc.createVariable('Hwav', 'f4', ('time','x',))
        Hwav_var.long_name = r'significant wave height'
        Hwav_var.units = 'm'
        Hwav_var._FillValue = np.float32(1e20)
        Hwav_var[:] = uhydro_out['Hwav']
        tau_var = nc.createVariable('tau', 'f4', ('time','x',))
        tau_var.long_name = r'bottom shear stress'
        tau_var.units = 'Pa'
        tau_var._FillValue = np.float32(1e20)
        tau_var[:] = uhydro_out['tau']
        Css_var = nc.createVariable('Css', 'f4', ('time','x',))
        Css_var.long_name = r'suspended sediment concentration'
        Css_var.units = 'kg/m3'
        Css_var._FillValue = np.float32(1e20)
        Css_var[:] = uhydro_out['Css']
        Cj_var = nc.createVariable('Cj', 'f4', ('time','x',))
        Cj_var.long_name = r'water salinity'
        Cj_var.units = 'PSU'
        Cj_var._FillValue = np.float32(1e20)
        Cj_var[:] = uhydro_out['Cj']
    finally:
        nc.close()
        
def write_ecogeom_outputs(odir, sid, ecogeom_out):
    """Write model outputs into a nc file.
    Arguments:
        odir : output directory
        sid : site id
        ecogeom_out : eco-geomorphology outputs at user-defined time interval 
    Returns : 
    """
    nx = np.shape(ecogeom_out['zh'])[1]
    filename = odir + '{:d}'.format(sid) + '_ecogeom.nc'
    try:
        nc = netcdf.netcdf_file(filename, 'w', version=2)
        nc.history = r'MACES simulated eco-geomorphology'
        nc.contact = r'Zeli.Tan@pnnl.gov'
        nc.createDimension('time', None)
        nc.createDimension('pool', 2)
        nc.createDimension('x', nx)  
        # create and write variables
        pft_var = nc.createVariable('pft', 'i1', ('time','x',))
        pft_var.long_name = r'platform plant function type'
        pft_var.units = '0 to 8'
        pft_var._FillValue = np.int8(-1)
        pft_var[:] = ecogeom_out['pft'] 
        zh_var = nc.createVariable('zh', 'f4', ('time','x',))
        zh_var.long_name = r'platform surface elevation'
        zh_var.units = 'msl'
        zh_var._FillValue = np.float32(1e20)
        zh_var[:] = ecogeom_out['zh']
        Esed_var = nc.createVariable('Esed', 'f4', ('time','x',))
        Esed_var.long_name = r'sediment erosion rate'
        Esed_var.units = 'kg/m2/s'
        Esed_var._FillValue = np.float32(1e20)
        Esed_var[:] = ecogeom_out['Esed']
        Dsed_var = nc.createVariable('Dsed', 'f4', ('time','x',))
        Dsed_var.long_name = r'suspended sediment deposition rate'
        Dsed_var.units = 'kg/m2/s'
        Dsed_var._FillValue = np.float32(1e20)
        Dsed_var[:] = ecogeom_out['Dsed']
        Lbed_var = nc.createVariable('Lbed', 'f4', ('time','x',))
        Lbed_var.long_name = r'sand bed load rate'
        Lbed_var.units = 'kg/m2/s'
        Lbed_var._FillValue = np.float32(1e20)
        Lbed_var[:] = ecogeom_out['Lbed']
        DepOM_var = nc.createVariable('DepOM', 'f4', ('time','x',))
        DepOM_var.long_name = r'Organic matter deposition rate'
        DepOM_var.units = 'kg/m2/s'
        DepOM_var._FillValue = np.float32(1e20)
        DepOM_var[:] = ecogeom_out['DepOM']
        Bag_var = nc.createVariable('Bag', 'f4', ('time','x',))
        Bag_var.long_name = r'platform aboveground biomass'
        Bag_var.units = 'kg/m2'
        Bag_var._FillValue = np.float32(1e20)
        Bag_var[:] = ecogeom_out['Bag']
        Bbg_var = nc.createVariable('Bbg', 'f4', ('time','x',))
        Bbg_var.long_name = r'platform belowground biomass'
        Bbg_var.units = 'kg/m2'
        Bbg_var._FillValue = np.float32(1e20)
        Bbg_var[:] = ecogeom_out['Bbg']
        SOM_var = nc.createVariable('SOM', 'f4', ('time','x','pool',))
        SOM_var.long_name = r'platform column-integrated soil organic matter'
        SOM_var.units = 'kg/m2'
        SOM_var._FillValue = np.float32(1e20)
        SOM_var[:] = ecogeom_out['SOM']
    finally:
        nc.close()
    
def run_tai_maces(input_data, models, spinup, verbose):
    """Write model outputs into a nc file.
    Arguments:
        odir : output directory
        sid : site id
        spinup : True = spinup, otherwise regular
    Returns : 
        tai_state : model state variables
        uhydro_out : hydrodynamic archives
        ecogeom_out : eco-geomorphology archives
    """
    # input settings
    rk4_mode = input_data['rk4_mode']
    date0 = input_data['date0']
    date1 = input_data['date1']
    site_x = input_data['x']
    site_pft = input_data['pft']
    site_zh = input_data['zh']
    site_Bag = input_data['Bag']
    site_Bbg = input_data['Bbg']
    site_OM = input_data['SOM']
    site_Esed = input_data['Esed']
    site_Dsed = input_data['Dsed']
    site_MHT = input_data['MHT']
    U10 = input_data['U10']
    Twav = input_data['Twav']
    h0 = input_data['h0']
    U0 = input_data['U0']
    Css0 = input_data['Css0']
    Cj0 = input_data['Cj0']
    rslr = input_data['rslr']
    uhydro_tol = input_data['tol']
    dyncheck = input_data['dyncheck']
    # hydrodynamic and eco-geomorphology model instances
    nx = np.size(site_x)
    dx = np.zeros(nx, dtype=np.float64)
    for ii in range(nx):
        if ii==0:
            dx[ii] = 0.5*(site_x[ii+1]-site_x[ii])
        elif ii==nx-1:
            dx[ii] = 0.5*(site_x[ii]-site_x[ii-1])
        else:
            dx[ii] = 0.5*(site_x[ii+1]-site_x[ii-1])
    taihydro = models['taihydro']
    mac_mod = models['mac_mod']
    omac_mod = models['omac_mod']
    rhoSed = mac_mod.m_params['rhoSed']
    porSed = mac_mod.m_params['porSed']
    rhoOM = omac_mod.m_params['rhoOM']
    # output variables
    jdn = get_julian_from_date(date0.year, date0.month, date0.day)
    if input_data['tstep']=='day':
        ntime = (date1-date0).days
    elif input_data['tstep']=='month':
        ntime = 12*(date1.year-date0.year) + (date1.month-date0.month)
    elif input_data['tstep']=='year':
        ntime = date1.year - date0.year
    nday = (date1-date0).days
    nhour = 24 * nday
    npool = np.shape(site_OM)[1]
    uhydro_out = {}
    if not spinup:
        uhydro_out['h'] = 1e20 * np.ones((nhour,nx), dtype=np.float32)
        uhydro_out['U'] = 1e20 * np.ones((nhour,nx), dtype=np.float32)
        uhydro_out['Hwav'] = 1e20 * np.ones((nhour,nx), dtype=np.float32)
        uhydro_out['tau'] = 1e20 * np.ones((nhour,nx), dtype=np.float32)
        uhydro_out['Css'] = 1e20 * np.ones((nhour,nx), dtype=np.float32)
        uhydro_out['Cj'] = 1e20 * np.ones((nhour,nx), dtype=np.float32)
    ecogeom_out = {}
    if not spinup:
        ecogeom_out['pft'] = -1 * np.ones((ntime,nx), dtype=np.int8)
        ecogeom_out['zh'] = 1e20 * np.ones((ntime,nx), dtype=np.float32)
        ecogeom_out['Esed'] = 1e20 * np.ones((ntime,nx), dtype=np.float32)
        ecogeom_out['Dsed'] = 1e20 * np.ones((ntime,nx), dtype=np.float32)
        ecogeom_out['Lbed'] = 1e20 * np.ones((ntime,nx), dtype=np.float32)
        ecogeom_out['DepOM'] = 1e20 * np.ones((ntime,nx), dtype=np.float32)
        ecogeom_out['Bag'] = 1e20 * np.ones((ntime,nx), dtype=np.float32)
        ecogeom_out['Bbg'] = 1e20 * np.ones((ntime,nx), dtype=np.float32)
        ecogeom_out['SOM'] = 1e20 * np.ones((ntime,nx,npool), dtype=np.float32)
    t = 0
    tf = 8.64e4 * nday
    hindx = -1
    dindx = -1
    tindx = -1
    ncount = 0
    isHourNode = False
    isTimeNode = False
    curstep = 50.0
    nextstep = MAX_OF_STEP
    # start simulation
    tstep = 0.0
    Esed_tot = np.zeros(nx, dtype=np.float64, order='F')
    Dsed_tot = np.zeros(nx, dtype=np.float64, order='F')
    Lbed_tot = np.zeros(nx, dtype=np.float64, order='F')
    DepOM_tot = np.zeros(nx, dtype=np.float64, order='F')
    Bag_tot = np.zeros(nx, dtype=np.float64, order='F')
    Bbg_tot = np.zeros(nx, dtype=np.float64, order='F')
    while t < tf:
        if t>=3.6e3*(hindx+1) and hindx+1<nhour:
            isHourNode = True
            hindx = hindx + 1
            if verbose and spinup:
                print('spinup time step', int(hindx))
            elif verbose:
                print('regular time step', int(hindx))
            if np.mod(hindx,24)==0:
                dindx = dindx + 1
                if dindx>0:
                    year, month, day = get_date_from_julian(jdn+dindx)
                    if input_data['tstep']=='day':
                        isTimeNode = True
                        tindx = tindx + 1
                    elif input_data['tstep']=='month' and day==1:
                        isTimeNode = True
                        tindx = tindx + 1
                    elif input_data['tstep']=='year' and month==1 and day==1:
                        isTimeNode = True
                        tindx = tindx + 1
        # simulate hydrodynamics
        if isHourNode:
            h0_abs = -site_zh[0] + h0[hindx]
            Hwav0 = estimate_Hwav_seaward(U10[dindx], h0_abs)
            #np.set_printoptions(precision=3, suppress=True)
            np.set_printoptions(precision=3, suppress=False)
            #print(np.array([U10[dindx],h0_abs,U0[hindx],Hwav0,Css0[dindx],Cj0[dindx]]))
        taihydro.modelsetup(site_zh, site_pft, site_Bag, site_Esed, site_Dsed, 
                            Twav, U10[dindx], h0_abs, U0[hindx], Hwav0, 
                            Css0[dindx], Cj0[dindx])
        curstep, nextstep, error = taihydro.modelrun(rk4_mode, uhydro_tol, 
                                                     dyncheck, curstep)
        assert error==0, "runge-Kutta iteration is more than MAXITER"
        taihydro.modelcallback()
        sim_h, sim_U, sim_Hwav, sim_tau, sim_Css, sim_Cj = taihydro.getmodelsims(nx)
        assert np.all(np.isfinite(sim_h)), "NaN h found"
        assert np.all(np.isfinite(sim_U)), "NaN U found"
        assert np.all(np.isfinite(sim_Hwav)), "NaN Hwav found"
        assert np.all(np.isfinite(sim_tau)), "NaN tau found"
        assert np.all(np.isfinite(sim_Css)), "NaN Css found"
        assert np.all(np.isfinite(sim_Cj)), "NaN Cj found"
        # get wet area length
        if isTimeNode:
            tmp_zh = site_zh - h0[hindx]
            wetL_potential = site_x[tmp_zh<0][-1]
            vol_water = np.sum(dx[sim_h>0]*sim_h[sim_h>0])
            #print(np.array([h0_abs, U0[hindx], sim_h[1], sim_U[1]]))
            #print(np.array([wetL_potential, vol_water, h0_abs, sim_h[0]]))
            print(site_zh)
        # simulate eco-geomorphology
        mac_inputs = {'x': site_x, 'pft': site_pft, 'Css': sim_Css, 'tau': sim_tau}
        site_Esed = mac_mod.mineral_suspension(mac_inputs)
        site_Dsed = mac_mod.mineral_deposition(mac_inputs)
        site_Lbed = mac_mod.bed_loading(mac_inputs)
        omac_inputs = {'x': site_x, 'zh': site_zh, 'MHT': site_MHT, 
                       'pft': site_pft, 'SOM': site_OM}
        site_Bag = omac_mod.aboveground_biomass(omac_inputs)
        omac_inputs['Bag'] = site_Bag
        site_Bbg = omac_mod.belowground_biomass(omac_inputs)
        site_DepOM = omac_mod.organic_deposition(omac_inputs)
        site_DecayOM = omac_mod.soilcarbon_decay(omac_inputs)
        # update soil OM pool
        DepOM_pools = np.zeros((nx,npool), dtype=np.float64, order='F')
        DepOM_pools[:,0] = 0.158 * site_DepOM
        DepOM_pools[:,1] = 0.842 * site_DepOM
        site_OM = site_OM + (DepOM_pools - site_DecayOM) * curstep
        # update platform elevation
        if not spinup:
            site_zh = site_zh + ((site_Dsed/rhoSed + site_Lbed/rhoSed + \
                site_DepOM/rhoOM - site_Esed/rhoSed)/(1.0-porSed) - \
                rslr) * curstep
        Esed_tot = Esed_tot + site_Esed * curstep
        Dsed_tot = Dsed_tot + site_Dsed * curstep
        Lbed_tot = Lbed_tot + site_Lbed * curstep
        DepOM_tot = DepOM_tot + site_DepOM * curstep
        Bag_tot = Bag_tot + site_Bag * curstep
        Bbg_tot = Bbg_tot + site_Bbg * curstep
        # archive hydrodynamic state variables
        if isHourNode and (not spinup):
            uhydro_out['h'][hindx] = sim_h
            uhydro_out['U'][hindx] = sim_U
            uhydro_out['Hwav'][hindx] = sim_Hwav
            uhydro_out['tau'][hindx] = sim_tau
            uhydro_out['Css'][hindx] = sim_Css
            uhydro_out['Cj'][hindx] = sim_Cj
        if isTimeNode and (not spinup):
            # archive daily mean eco-geomorphology variables
            ecogeom_out['zh'][tindx] = site_zh
            ecogeom_out['SOM'][tindx] = site_OM
            ecogeom_out['pft'][tindx] = site_pft
            ecogeom_out['Esed'][tindx] = Esed_tot / tstep
            ecogeom_out['Dsed'][tindx] = Dsed_tot / tstep
            ecogeom_out['Lbed'][tindx] = Lbed_tot / tstep
            ecogeom_out['DepOM'][tindx] = DepOM_tot / tstep
            ecogeom_out['Bag'][tindx] = Bag_tot / tstep
            ecogeom_out['Bbg'][tindx] = Bbg_tot / tstep
            # reset temporary memory
            tstep = 0.0
            Esed_tot[:] = 0.0 
            Dsed_tot[:] = 0.0
            Lbed_tot[:] = 0.0
            DepOM_tot[:] = 0.0
            Bag_tot[:] = 0.0
            Bbg_tot[:] = 0.0
        # check small time step
        isHourNode = False
        isTimeNode = False
        if curstep<0.1:
            ncount = ncount + 1
            err_msg = 'run diverge at step ' + '{:d}'.format(hindx)
            assert ncount<=100, err_msg
            nextstep = 50.0
        else:
            ncount = 0
        t = t + curstep
        tstep = tstep + curstep
        curstep = nextstep
        nextstep = MAX_OF_STEP
    # returns
    tai_state = {'pft': site_pft, 'zh': site_zh, 'Bag': site_Bag, 
                 'Bbg': site_Bbg, 'SOM': site_OM, 'Esed': site_Esed, 
                 'Dsed': site_Dsed}
    return tai_state, uhydro_out, ecogeom_out