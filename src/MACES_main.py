#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 12 22:08:10 2019

Main program to run the MACES model

@author: Zeli Tan
"""

import sys
import importlib
import TAIMODSuper
import pandas as pd
import numpy as np
import maces_utilities as utils
import maces_coupler as cpl
from datetime import date
from mpi4py import MPI
from optparse import OptionParser
from TAIHydroMOD import tai_hydro_mod as taihydro
        
if __name__=='__main__':
    
    # use the OptionParser to parse the command line options
    parser = OptionParser()
    parser.add_option("-f", "--file", type="string", dest="filename", 
                      default='namelist.maces.xml')
    (options, args) = parser.parse_args()
    
    # MPI commands
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    numprocs = comm.Get_size()
    
    if rank == 0:
        master_process = True
    else:
        master_process = False
        
    # set how the run status is shown
    #np.set_printoptions(precision=3, suppress=True)
    np.set_printoptions(precision=3, suppress=False)
        
    # read MACES namelist xml file and parameter xml files
    if master_process:
        # MACES namelist xml file
        namelist = utils.parseXML_namelist(options.filename)
        # MINAC model xml file
        xmlfile = namelist['MINAC_FILE']
        mac_name = namelist['MINAC_TYPE']
        mac_params = utils.parseXML_params(xmlfile, mac_name)
        # OMAC model xml file
        xmlfile = namelist['OMAC_FILE']
        omac_name = namelist['OMAC_TYPE']
        omac_params = utils.parseXML_params(xmlfile, omac_name)
        # WAVERO model xml file
        xmlfile = namelist['WAVERO_FILE']
        wavero_name = namelist['WAVERO_TYPE']
        wavero_params = utils.parseXML_params(xmlfile, wavero_name)
        # LNDMGR model xml file
        xmlfile = namelist['LNDMGR_FILE']
        lndmgr_name = namelist['LNDMGR_TYPE']
        lndmgr_params = utils.parseXML_params(xmlfile, lndmgr_name)
        # hydrodynamic model xml file
        xmlfile = namelist['HYDRO_FILE']
        hydro_params = utils.parseXML_hydro_params(xmlfile)
    else:
        namelist = None
        mac_params = None
        omac_params = None
        wavero_params = None
        lndmgr_params = None
        hydro_params = None
    namelist = comm.bcast(namelist, root=0)
    mac_params = comm.bcast(mac_params, root=0)
    omac_params = comm.bcast(omac_params, root=0)
    wavero_params = comm.bcast(wavero_params, root=0)
    lndmgr_params = comm.bcast(lndmgr_params, root=0)
    hydro_params = comm.bcast(hydro_params, root=0)
    # make parameters of different components consistent
    if 'alphaA' in mac_params:
        hydro_params['alphaA'] = mac_params['alphaA']
    if 'betaA' in mac_params:
        hydro_params['betaA'] = mac_params['betaA']
    if 'alphaD' in mac_params:
        hydro_params['alphaD'] = mac_params['alphaD']
    if 'betaD' in mac_params:
        hydro_params['betaD'] = mac_params['betaD']
    if 'Bmax' in omac_params:
        mac_params['Bmax'] = omac_params['Bmax']
    
    # read site database excel file
    if master_process:
        df = pd.read_excel(namelist['SITE_FILE'], sheet_name="diva", header=0, 
                           usecols="A:AV")
    else:
        df = None
    df = comm.bcast(df, root=0)
    site_ids = np.array(df['DIVA_ID'], dtype=np.int32)
    site_coastline = np.array(df['coastline'], dtype=np.float64)    # km
    site_fetchagl = np.array(df['fetchagl'], dtype=np.float64)      # degree
    site_trng = np.array(df['mtidalrng'], dtype=np.float64)         # m
    site_mhws = np.array(df['mhws'], dtype=np.float64)              # m
    site_mhwn = np.array(df['mhwn'], dtype=np.float64)              # m
    site_uplift = np.array(df['uplift'], dtype=np.float64)          # mm/yr
    site_TSM = 1e-3 * np.array(df['TSM'], dtype=np.float64)         # kg/m3
    site_sal = np.array(df['salinity'], dtype=np.float64)           # PSU
    nsite = np.size(site_ids)
    site_n = min(namelist['LAST_ID'], nsite)
    site_1 = namelist['FIRST_ID'] - 1
    nrun = site_n - site_1
    site_ids = site_ids[site_1:site_n]
    site_coastline = site_coastline[site_1:site_n]
    site_fetchagl = site_fetchagl[site_1:site_n]
    site_trng = site_trng[site_1:site_n]
    site_mhws = site_mhws[site_1:site_n]
    site_mhwn = site_mhwn[site_1:site_n]
    site_uplift = site_uplift[site_1:site_n]
    site_TSM = site_TSM[site_1:site_n]
    site_sal = site_sal[site_1:site_n]
    npft = TAIMODSuper.npft
    
    diva_segments = []
    pft_segments = []
    pft_orders = []
    for ii, site in enumerate(np.arange(site_1,site_n)):
        area = np.zeros(utils.NTOPSEG, dtype=np.float64)
        for jj in range(utils.NTOPSEG):
            code = 'area' + '{:02d}'.format(jj+1)
            area[jj] = float(df[code][site])
        segments = area / site_coastline[ii]
        diva_segments.append(segments)
        segments = np.zeros(npft, dtype=np.float64)
        orders = np.zeros(npft, dtype=np.int32)
        for jj in range(npft):
            code = 'pft' + '{:d}'.format(jj)
            segments[jj] = float(df[code][site])
            code = 'pft' + '{:d}'.format(jj) + '_order'
            orders[jj] = int(df[code][site])
            assert (segments[jj]>0)==(orders[jj]>=0), "inconsistent pft segment found"
        pft_segments.append(segments)
        pft_orders.append(orders)
        
    # read driving data
    # units: SLR (mm/yr), Tair (K), U10 (m/s), h0 (m), U0 (m/s), 
    #        Hwav0 (m), Twav (s)
    sid_range = [site_1, site_n]
    date0_str = namelist['RUN_STARTDATE'].split('-')
    date1_str = namelist['RUN_STOPDATE'].split('-')
    run_date0 = date(int(date0_str[0]), int(date0_str[1]), int(date0_str[2]))
    run_date1 = date(int(date1_str[0]), int(date1_str[1]), int(date1_str[2]))
    if master_process:
        SLR = utils.read_force_data(namelist['FILE_SLR'], 'SLR', \
            run_date0, run_date1, namelist['SLR_TSTEP'], 'year', sid_range)
        Tair_avg = utils.read_force_data(namelist['FILE_Tair'], 'Tmean', \
            run_date0, run_date1, 1, 'year', sid_range)
        Tair_summer = utils.read_force_data(namelist['FILE_Tair'], 'Tsummer', \
            run_date0, run_date1, 1, 'year', sid_range)
        U10 = utils.read_force_data(namelist['FILE_U10'], 'U10', \
            run_date0, run_date1, namelist['U10_TSTEP'], 'minute', sid_range)
        h0 = utils.read_force_data(namelist['FILE_h'], 'h', \
            run_date0, run_date1, namelist['h_TSTEP'], 'minute', sid_range)
        Twav = utils.read_force_data(namelist['FILE_Wave'], 'Twav', \
            run_date0, run_date1, namelist['Wave_TSTEP'], 'minute', sid_range)
        if len(namelist['FILE_SSC'])>0:
            SSC = 1e-3 * utils.read_force_data(namelist['FILE_SSC'], 'TSM', \
                run_date0, run_date1, namelist['SSC_TSTEP'], 'minute', sid_range)
        else:
            nt_ssc = np.shape(h0)[0] * int(namelist['SSC_TSTEP']/namelist['h_TSTEP'])
            nsite_ssc = np.shape(h0)[1]
            SSC = np.zeros((nt_ssc,nsite_ssc))
            for ii, site in enumerate(np.arange(site_1,site_n)):
                SSC[:,ii] = site_TSM[ii]
    else:
        SLR = None
        Tair_avg = None
        Tair_summer = None
        U10 = None
        h0 = None
        Twav = None
        SSC = None
    SLR = comm.bcast(SLR, root=0)
    Tair_avg = comm.bcast(Tair_avg, root=0)
    Tair_summer = comm.bcast(Tair_summer, root=0)
    U10 = comm.bcast(U10, root=0)
    h0 = comm.bcast(h0, root=0)
    Twav = comm.bcast(Twav, root=0)
    SSC = comm.bcast(SSC, root=0)
        
    # load ecogeomorphology modules
    mac_module = importlib.import_module('minac_mod')
    mac_class = getattr(mac_module, namelist['MINAC_TYPE'])
    
    omac_module = importlib.import_module('omac_mod')
    omac_class = getattr(omac_module, namelist['OMAC_TYPE'])
    
    wavero_module = importlib.import_module('wavero_mod')
    wavero_class = getattr(wavero_module, namelist['WAVERO_TYPE'])
    
    lndmgr_module = importlib.import_module('lndmgr_mod')
    lndmgr_class = getattr(lndmgr_module, namelist['LNDMGR_TYPE'])
        
    # run simulations (in each iteration, a processor tests the sensitivity of 
    # one site for all parameters)
    niter = int( np.ceil( float(nrun) / float(numprocs) ) )
    for ii in range(niter):
        iid = np.mod( ii*numprocs + rank, nrun )
        site_id = site_ids[iid]
        print( "Simulate site ", site_id )
        sys.stdout.flush()
        
        try:
            # construct site platform
            xres = namelist['CELL_RES']
            xnum = namelist['CELL_NUM']
            site_x, site_zh, site_fetch = \
                utils.construct_tai_platform(diva_segments[iid], site_coastline[iid],
                                             site_fetchagl[iid], xres, xnum)
            nx = len(site_x)
            site_dx = np.zeros(nx, dtype=np.float64, order='F')
            for jj in range(nx):
                if jj==0:
                    site_dx[jj] = 0.5*(site_x[jj+1]-site_x[jj])
                elif jj==nx-1:
                    site_dx[jj] = 0.5*(site_x[jj]-site_x[jj-1])
                else:
                    site_dx[jj] = 0.5*(site_x[jj+1]-site_x[jj-1])
            coords = {'x': site_x, 'dx': site_dx}
            
            # construct pft distribution
            orders = pft_orders[iid]
            segments = pft_segments[iid]
            pfts = np.arange(npft)
            indice = orders>=0
            orders = orders[indice]
            segments = segments[indice]
            pfts = pfts[indice]
            indices = sorted(range(len(orders)), key=lambda k: orders[k])
            segments = segments[indices]
            pfts = pfts[indices]
            site_pft = utils.construct_platform_pft(segments, pfts, site_x)
            
            # instantiate hydrodynamics model
            nvar = len(namelist['HYDRO_TOL'])
            taihydro.inithydromod(site_x, site_zh, site_fetch, site_TSM[iid], 
                                  nvar, npft)
            taihydro.setmodelparams(hydro_params['d50'], hydro_params['Cz0'], 
                                    hydro_params['Kdf'], hydro_params['cbc'], 
                                    hydro_params['cwc'], hydro_params['fr'], 
                                    hydro_params['alphaA'], hydro_params['betaA'], 
                                    hydro_params['alphaD'], hydro_params['betaD'], 
                                    hydro_params['cD0'], hydro_params['ScD'])
        
            # instantiate ecogeomorphology models
            mac_mod = mac_class(mac_params)
            omac_mod = omac_class(omac_params)
            wavero_mod = wavero_class(wavero_params)
            lndmgr_mod = lndmgr_class(lndmgr_params)
            
            models = {'taihydro': taihydro, 'mac_mod': mac_mod, 
                      'omac_mod': omac_mod, 'wavero_mod': wavero_mod, 
                      'lndmgr_mod': lndmgr_mod}
        
            # first run the spin-up
            site_Bag = np.zeros(nx, dtype=np.float64, order='F')
            site_Bbg = np.zeros(nx, dtype=np.float64, order='F')
            tai_state = {'pft': site_pft, 'zh': site_zh, 'Bag': site_Bag, 
                         'Bbg': site_Bbg}
            
            rslr = SLR[:,iid] - site_uplift[iid]
            sal = site_sal[iid]
            forcings = {'U10': U10[:,iid], 'Tmean': Tair_avg[:,iid], 
                        'Tsummer': Tair_summer[:,iid], 'h0': h0[:,iid],
                        'Twav': Twav[:,iid], 'Cs0': SSC[:,iid], 'sal': sal, 
                        'rslr': rslr, 'trng': site_trng[iid], 'mhws': site_mhws[iid], 
                        'mhwn': site_mhwn[iid], 'refCss': site_TSM[iid]}
            
            input_data = {'coord': coords, 'state': tai_state, 
                          'forcings': forcings, 'namelist': namelist}
            tai_state, __, __ = cpl.run_tai_maces(input_data, models, True)
            
            # then do the formal run
            input_data = {'coord': coords, 'state': tai_state, 
                          'forcings': forcings, 'namelist': namelist}
            __, uhydro_out, ecogeom_out = cpl.run_tai_maces(input_data, \
                models, False)
            
        except AssertionError as errstr:
            # print error message and exit the program
            print("Model stops due to that", errstr)
            taihydro.finalizehydromod()
            sys.exit()
            
        # deallocate
        taihydro.finalizehydromod()
            
        # archive outputs
        if namelist['OUTPUT_HYDRO']:
            filename_hydro = namelist['DOUT_ROOT'] + '/maces_hydro_' + \
                namelist['RUN_STARTDATE'] + '_' + namelist['RUN_STOPDATE'] + \
                '_' + namelist['CASE'] + '_' + '{:d}'.format(site_id) + '.nc'
            utils.write_hydro_outputs(filename_hydro, namelist['HYDRO_TSTEP'], 
                                      uhydro_out)
        filename_ecogeom = namelist['DOUT_ROOT'] + '/maces_ecogeom_' + \
            namelist['RUN_STARTDATE'] + '_' + namelist['RUN_STOPDATE'] + \
            '_' + namelist['CASE'] + '_' + '{:d}'.format(site_id) + '.nc'
        utils.write_ecogeom_outputs(filename_ecogeom, namelist['ECOGEOM_TSTEP'], 
                                    ecogeom_out)
