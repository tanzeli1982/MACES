#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 21 00:13:18 2019

Simulation coupler

@author: Zeli Tan
"""

import sys
import numpy as np
import maces_utilities as utils
from datetime import date

MAX_OF_STEP = 1800  # maximum simulation time step (s)
rk4_mode = 101      # Runge-kutta iteration mode

def run_tai_maces(input_data, models, spinup):
    """Write model outputs into a nc file.
    Arguments:
        input_data : various input data
        models : model objects
        spinup : True = spinup, otherwise regular
    Returns : 
        tai_state : model state variables
        uhydro_out : hydrodynamic archives
        ecogeom_out : eco-geomorphology archives
    """
    # namelist settings
    namelist = input_data['namelist']
    verbose = namelist['Verbose']
    uhydro_tol = namelist['HYDRO_TOL']
    dyncheck = namelist['DYN_CHECK']
    date0_str = namelist['RUN_STARTDATE'].split('-')
    date1_str = namelist['RUN_STOPDATE'].split('-')
    date0 = date(int(date0_str[0]), int(date0_str[1]), int(date0_str[2]))
    if not spinup:
        date1 = date(int(date1_str[0]), int(date1_str[1]), int(date1_str[2]))
    else:
        date1 = utils.get_spinup_stop_date(date0, namelist['SPINUP_OPTION'], 
                                           namelist['SPINUP_N'])
    
    # input settings
    x = input_data['coord']['x']
    #site_dx = input_data['coord']['dx']
    pft = input_data['state']['pft']
    zh = input_data['state']['zh']
    Bag = input_data['state']['Bag']
    Bbg = input_data['state']['Bbg']
    OM = input_data['state']['OM']
    trng = input_data['forcings']['trng']
    mhws = input_data['forcings']['mhws']
    U10 = input_data['forcings']['U10']
    Twav = input_data['forcings']['Twav']
    Tair = input_data['forcings']['Tair']
    h0 = input_data['forcings']['h0']
    Cs0 = input_data['forcings']['Cs0']
    sal = input_data['forcings']['sal']
    rslr = input_data['forcings']['rslr']
    nx = len(x)
    npool = np.shape(OM)[1]
    xref = utils.get_refshore_coordinate(x, zh)
    
    # hydrodynamic and eco-geomorphology model objects
    taihydro = models['taihydro']
    mac_mod = models['mac_mod']
    omac_mod = models['omac_mod']
    wavero_mod = models['wavero_mod']
    lndmgr_mod = models['lndmgr_mod']
    rhoSed = mac_mod.m_params['rhoSed']
    porSed = mac_mod.m_params['porSed']
    rhoOM = omac_mod.m_params['rhoOM']
    wave_mod = namelist['WAVE_TYPE']
    
    # output variables
    jdn = utils.get_julian_from_date(date0.year, date0.month, date0.day)
    nday = (date1 - date0).days
    nhour = 24 * nday
    nt_hydro = utils.get_shr_output_num(date0, date1, namelist['HYDRO_TSTEP'])
    nt_ecogeom = utils.get_lng_output_num(date0, date1, namelist['ECOGEOM_TSTEP'])
    uhydro_out = {}
    if (not spinup) and (nt_hydro>0):
        uhydro_out['x'] = np.reshape(np.float32(x),(1,nx))
        uhydro_out['h'] = 1e20 * np.ones((1,nt_hydro,nx), dtype=np.float32)
        uhydro_out['U'] = 1e20 * np.ones((1,nt_hydro,nx), dtype=np.float32)
        uhydro_out['Hwav'] = 1e20 * np.ones((1,nt_hydro,nx), dtype=np.float32)
        uhydro_out['Uwav'] = 1e20 * np.ones((1,nt_hydro,nx), dtype=np.float32)
        uhydro_out['tau'] = 1e20 * np.ones((1,nt_hydro,nx), dtype=np.float32)
        uhydro_out['Css'] = 1e20 * np.ones((1,nt_hydro,nx), dtype=np.float32)
        uhydro_out['Esed'] = 1e20 * np.ones((1,nt_hydro,nx), dtype=np.float32)
        uhydro_out['Dsed'] = 1e20 * np.ones((1,nt_hydro,nx), dtype=np.float32)
    ecogeom_out = {}
    ecogeom_tot = np.zeros(nt_ecogeom, dtype=np.float32)
    if (not spinup) and (nt_ecogeom>0):
        ecogeom_out['x'] = np.reshape(np.float32(x),(1,nx))
        ecogeom_out['pft'] = -1 * np.ones((1,nt_ecogeom,nx), dtype=np.int8)
        ecogeom_out['zh'] = 1e20 * np.ones((1,nt_ecogeom,nx), dtype=np.float32)
        ecogeom_out['Esed'] = np.zeros((1,nt_ecogeom,nx), dtype=np.float32)
        ecogeom_out['Dsed'] = np.zeros((1,nt_ecogeom,nx), dtype=np.float32)
        ecogeom_out['Lbed'] = np.zeros((1,nt_ecogeom,nx), dtype=np.float32)
        ecogeom_out['DepOM'] = np.zeros((1,nt_ecogeom,nx), dtype=np.float32)
        ecogeom_out['Bag'] = np.zeros((1,nt_ecogeom,nx), dtype=np.float32)
        ecogeom_out['Bbg'] = np.zeros((1,nt_ecogeom,nx), dtype=np.float32)
        ecogeom_out['OM'] = 1e20 * np.ones((1,nt_ecogeom,nx,npool), dtype=np.float32)
        
    # temporal variables
    Esed = np.zeros(nx, dtype=np.float64, order='F')
    Dsed = np.zeros(nx, dtype=np.float64, order='F')
    Lbed = np.zeros(nx, dtype=np.float64, order='F')
    DepOM = np.zeros(nx, dtype=np.float64, order='F')
    DecayOM = np.zeros((nx,npool), dtype=np.float64, order='F')
    sources = np.zeros(nx, dtype=np.float64, order='F')
    sinks = np.zeros(nx, dtype=np.float64, order='F')
    tau_old = np.zeros(nx, dtype=np.float64, order='F')
    
    # start simulation
    t = 0.0
    tf = 8.64e4 * nday
    hindx = -1
    dindx = -1
    ncount = 0
    curstep = 50.0
    nextstep = MAX_OF_STEP
    hydro_indx = -1
    ecogeom_indx = -1
    while t <= tf:
        if t>=3.6e3*(hindx+1) and hindx+1<=nhour:
            hindx = hindx + 1
            if verbose and spinup:
                print('spinup time step', int(hindx))
                sys.stdout.flush()
            elif verbose:
                print('regular time step', int(hindx))
                sys.stdout.flush()
            if np.mod(hindx,24)==0:
                dindx = dindx + 1
                year, month, day = utils.get_date_from_julian(jdn+dindx)
                date_cur = date(year, month, day)
                doy = min(date_cur.timetuple().tm_yday, 365)
                slope = utils.get_platform_slope(x, zh)
                        
        # get instant boundary conditions
        indx = utils.get_forcing_index(t, 'minute', namelist['U10_TSTEP'])
        U10_inst = U10[indx]
        indx = utils.get_forcing_index(t, 'hour', namelist['Tair_TSTEP'])
        Tair_inst = Tair[indx]
        indx = utils.get_forcing_index(t, 'minute', namelist['h_TSTEP'])
        h0_inst = h0[indx] - zh[0]
        indx = utils.get_forcing_index(t, 'minute', namelist['Wave_TSTEP'])
        Twav_inst = Twav[indx]
        indx = int( (year-date0.year)/namelist['SLR_TSTEP'] )
        rslr_inst = rslr[indx]
        
        # simulate hydrodynamics
        if mac_mod.m_update_Css:
            sources[:] = Esed
            sinks[:] = Dsed
        else:
            sources[:] = 0.0
            sinks[:] = 0.0
        taihydro.modelsetup(sources, sinks, zh, pft, Bag, xref, Twav_inst,
                            h0_inst, U10_inst, Cs0)
        curstep, nextstep, error = taihydro.modelrun(rk4_mode, uhydro_tol, 
                                                     dyncheck, curstep)
        assert error==0, "runge-Kutta iteration is more than MAXITER"
        taihydro.modelcallback(wave_mod)
        assert np.all(np.isfinite(taihydro.sim_h)), "NaN h found"
        assert np.all(np.isfinite(taihydro.sim_u)), "NaN U found"
        assert np.all(np.isfinite(taihydro.sim_uwav)), "NaN Uwav found"
        assert np.all(np.isfinite(taihydro.sim_hwav)), "NaN Hwav found"
        assert np.all(np.isfinite(taihydro.sim_tau)), "NaN tau found"
        assert np.all(np.isfinite(taihydro.sim_css)), "NaN Css found"
        dtau = taihydro.sim_tau - tau_old
        tau_old[:] = taihydro.sim_tau
        
        # simulate mineral accretion
        mac_inputs = {'x': x, 'xref': xref, 'pft': pft, 'zh': zh, 
                      'Css': taihydro.sim_css, 'tau': taihydro.sim_tau, 
                      'U': taihydro.sim_u, 'h': taihydro.sim_h, 
                      'Uwav': taihydro.sim_uwav, 'Bag': Bag, 'Bbg': Bbg, 
                      'Esed': Esed, 'Dsed': Dsed, 'Lbed': Lbed, 'S': slope, 
                      'dtau': dtau, 'TR': trng, 'dt': curstep, 'refCss': Cs0}
        Esed = mac_mod.mineral_suspension(mac_inputs)
        Dsed = mac_mod.mineral_deposition(mac_inputs)
        Lbed = mac_mod.bed_loading(mac_inputs)
        
        # simulate organic matter accretion
        omac_inputs = {'x': x, 'zh': zh, 'S': slope, 'pft': pft, 'OM': OM, 
                       'Bag': Bag, 'Bbg': Bbg, 'DepOM': DepOM, 
                       'DecayOM': DecayOM, 'TR': trng, 'MHHW': mhws, 
                       'month': month, 'doy': doy, 'Tair': Tair_inst, 
                       'dt': curstep}
        Bag = omac_mod.aboveground_biomass(omac_inputs)
        Bbg = omac_mod.belowground_biomass(omac_inputs)
        DepOM = omac_mod.organic_deposition(omac_inputs)
        DecayOM = omac_mod.soilcarbon_decay(omac_inputs)
        # update soil OM pool
        DepOM_pools = np.zeros((nx,npool), dtype=np.float64, order='F')
        DepOM_pools[:,0] = 0.158 * DepOM
        DepOM_pools[:,1] = 0.842 * DepOM
        OM += (DepOM_pools - DecayOM) * curstep
        
        # simulate wave-driven lateral erosion
        wavero_inputs = {'x': x}
        x = wavero_mod.wave_erosion(wavero_inputs)
        
        # simulate landward migration
        lndmgr_inputs = {'pft': pft, 'sal': sal}
        pft = lndmgr_mod.landward_migration(lndmgr_inputs)
        
        # update platform elevation
        if not spinup:
            #zh = utils.update_platform_elev(zh, Esed, Dsed, Lbed, DepOM, \
            #          rhoSed, rhoOM, porSed, rslr_inst, curstep)
            xref = utils.get_refshore_coordinate(x, zh)
             
        # archive short-term hydrodynamic state variables
        if (not spinup) and (nt_hydro>0):
            indx = utils.get_shr_output_index(t, namelist['HYDRO_TSTEP'])
            if indx>hydro_indx:
                hydro_indx = indx
                uhydro_out['h'][0,indx] = taihydro.sim_h
                uhydro_out['U'][0,indx] = taihydro.sim_u
                uhydro_out['Hwav'][0,indx] = taihydro.sim_hwav
                uhydro_out['Uwav'][0,indx] = taihydro.sim_uwav
                uhydro_out['tau'][0,indx] = taihydro.sim_tau
                uhydro_out['Css'][0,indx] = taihydro.sim_css
                uhydro_out['Esed'][0,indx] = Esed
                uhydro_out['Dsed'][0,indx] = Dsed
        
        # archive long-term mean eco-geomorphology variables
        if (not spinup) and (nt_ecogeom>0):
            indx = utils.get_lng_output_index(date0, date_cur, \
                namelist['ECOGEOM_TSTEP'])
            ecogeom_tot[indx] = ecogeom_tot[indx] + curstep
            ecogeom_out['Esed'][0,indx] = ecogeom_out['Esed'][0,indx] + \
                Esed * curstep
            ecogeom_out['Dsed'][0,indx] = ecogeom_out['Dsed'][0,indx] + \
                Dsed * curstep
            ecogeom_out['Lbed'][0,indx] = ecogeom_out['Lbed'][0,indx] + \
                Lbed * curstep
            ecogeom_out['DepOM'][0,indx] = ecogeom_out['DepOM'][0,indx] + \
                DepOM * curstep
            ecogeom_out['Bag'][0,indx] = ecogeom_out['Bag'][0,indx] + \
                Bag * curstep
            ecogeom_out['Bbg'][0,indx] = ecogeom_out['Bbg'][0,indx] + \
                Bbg * curstep
            if indx>ecogeom_indx:
                ecogeom_indx = indx
                ecogeom_out['zh'][0,indx] = zh
                ecogeom_out['OM'][0,indx] = OM
                ecogeom_out['pft'][0,indx] = pft
            
        # check small time step
        if curstep<0.1:
            ncount = ncount + 1
            err_msg = 'run diverge at step ' + '{:d}'.format(hindx)
            assert ncount<=100, err_msg
            nextstep = 50.0
        else:
            ncount = 0
        t = t + curstep
        curstep = nextstep
        nextstep = MAX_OF_STEP
    
    # returns
    if (not spinup) and (nt_ecogeom>0):
        for jj in range(nx):
            ecogeom_out['Esed'][:,:,jj] = ecogeom_out['Esed'][:,:,jj]/ecogeom_tot
            ecogeom_out['Dsed'][:,:,jj] = ecogeom_out['Dsed'][:,:,jj]/ecogeom_tot
            ecogeom_out['Lbed'][:,:,jj] = ecogeom_out['Lbed'][:,:,jj]/ecogeom_tot
            ecogeom_out['DepOM'][:,:,jj] = ecogeom_out['DepOM'][:,:,jj]/ecogeom_tot
            ecogeom_out['Bag'][:,:,jj] = ecogeom_out['Bag'][:,:,jj]/ecogeom_tot
            ecogeom_out['Bbg'][:,:,jj] = ecogeom_out['Bbg'][:,:,jj]/ecogeom_tot
    tai_state = {'pft': pft, 'zh': zh, 'Bag': Bag, 'Bbg': Bbg, 'OM': OM}
    return tai_state, uhydro_out, ecogeom_out
