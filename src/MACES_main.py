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
    
    # read site database excel file
    if master_process:
        df = pd.read_excel(namelist['SITE_FILE'], sheet_name="diva", header=0, 
                           usecols="A:AU")
    else:
        df = None
    df = comm.bcast(df, root=0)
    site_ids = np.array(df['DIVA_ID'], dtype=np.int32)
    site_length = np.array(df['coastline'], dtype=np.float64)   # km
    site_trng = np.array(df['mtidalrng'], dtype=np.float64)     # m
    site_mhws = np.array(df['mhws'], dtype=np.float64)          # m
    site_mhwn = np.array(df['mhwn'], dtype=np.float64)          # m
    site_uplift = np.array(df['uplift'], dtype=np.float64)      # mm/yr
    site_TSM = 1e-3 * np.array(df['TSM'], dtype=np.float64)     # kg/m3
    site_sal = np.array(df['salinity'], dtype=np.float64)       # PSU
    nsite = np.size(site_ids)
    site_n = min(namelist['LAST_ID'], nsite)
    site_1 = namelist['FIRST_ID'] - 1
    nrun = site_n - site_1
    site_ids = site_ids[site_1:site_n]
    site_length = site_length[site_1:site_n]
    site_trng = site_trng[site_1:site_n]
    site_mhws = site_mhws[site_1:site_n]
    site_mhwn = site_mhwn[site_1:site_n]
    site_uplift = site_uplift[site_1:site_n]
    site_TSM = site_TSM[site_1:site_n]
    site_sal = site_sal[site_1:site_n]
    npft = TAIMODSuper.npft
    npool = TAIMODSuper.npool
    
    top_segments = []
    for jj in range(utils.NTOPSEG):
        code = 'area' + '{:02d}'.format(jj+1)
        segment = np.array(df[code],dtype=np.float64)[site_1:site_n] / site_length
        top_segments.append(segment)
    pft_segments = []
    pft_orders = []
    for jj in np.arange(npft):
        code = 'pft' + '{:d}'.format(jj)
        segment = np.array(df[code],dtype=np.float64)[site_1:site_n]
        pft_segments.append(segment)
        code = 'pft' + '{:d}'.format(jj) + '_order'
        order = np.array(df[code],dtype=np.int32)[site_1:site_n]
        pft_orders.append(order)
        
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
        Tair = utils.read_force_data(namelist['FILE_Tair'], 'Tair', \
            run_date0, run_date1, namelist['Tair_TSTEP'], 'hour', sid_range)
        U10 = utils.read_force_data(namelist['FILE_U10'], 'U10', \
            run_date0, run_date1, namelist['U10_TSTEP'], 'hour', sid_range)
        h0 = utils.read_force_data(namelist['FILE_h'], 'h', \
            run_date0, run_date1, namelist['h_TSTEP'], 'minute', sid_range)
        U0 = utils.read_force_data(namelist['FILE_U'], 'U', \
            run_date0, run_date1, namelist['U_TSTEP'], 'minute', sid_range)
        Hwav0 = utils.read_force_data(namelist['FILE_Wave'], 'Hwav', \
            run_date0, run_date1, namelist['Wave_TSTEP'], 'minute', sid_range)
        Twav = utils.read_force_data(namelist['FILE_Wave'], 'Twav', \
            run_date0, run_date1, namelist['Wave_TSTEP'], 'minute', sid_range)
    else:
        SLR = None
        Tair = None
        U10 = None
        h0 = None
        U0 = None
        Hwav0 = None
        Twav = None
    SLR = comm.bcast(SLR, root=0)
    Tair = comm.bcast(Tair, root=0)
    U10 = comm.bcast(U10, root=0)
    h0 = comm.bcast(h0, root=0)
    U0 = comm.bcast(U0, root=0)
    Hwav0 = comm.bcast(Hwav0, root=0)
    Twav = comm.bcast(Twav, root=0)
        
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
            diva_segments = np.zeros(utils.NTOPSEG, dtype=np.float64)
            for jj in range(utils.NTOPSEG):
                diva_segments[jj] = top_segments[jj][iid]
            site_x, site_zh = utils.construct_tai_platform(diva_segments, 
                                                           xres, xnum)
            nx = len(site_x)
            site_dx = np.zeros(nx, dtype=np.float64, order='F')
            for jj in range(nx):
                if jj==0:
                    site_dx[jj] = 0.5*(site_x[jj+1]-site_x[jj])
                elif jj==nx-1:
                    site_dx[jj] = 0.5*(site_x[jj]-site_x[jj-1])
                else:
                    site_dx[jj] = 0.5*(site_x[jj+1]-site_x[jj-1])
            xref = utils.get_refshore_coordinate(site_x, site_zh)
            coords = {'x': site_x, 'dx': site_dx, 'xref': xref}
            
            # construct pft distribution
            segments = []
            pfts = []
            orders = []
            for jj in range(npft):
                assert (pft_segments[jj][iid]>0)==(pft_orders[jj][iid]>=0), \
                    "inconsistent pft segment found"
                if pft_orders[jj][iid]>=0:
                    segments.append(pft_segments[jj][iid])
                    orders.append(pft_orders[jj][iid])
                    pfts.append(jj)
            indices = sorted(range(len(orders)), key=lambda k: orders[k])
            segments = np.array(segments)[indices]
            pfts = np.array(pfts)[indices]
            site_pft = utils.construct_platform_pft(segments, pfts, site_x)
            
            # instantiate hydrodynamics model
            nvar = len(namelist['HYDRO_TOL'])
            taihydro.inithydromod(site_x, site_zh, nvar, npft)
            taihydro.setmodelparams(hydro_params['d50'], hydro_params['Cz0'], 
                                    hydro_params['Kdf'], hydro_params['cbc'], 
                                    hydro_params['fr'], hydro_params['alphaA'], 
                                    hydro_params['betaA'], hydro_params['alphaD'], 
                                    hydro_params['betaD'], hydro_params['cD0'], 
                                    hydro_params['ScD'])
        
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
            site_OM = np.zeros((nx,npool), dtype=np.float64, order='F')
            tai_state = {'pft': site_pft, 'zh': site_zh, 'Bag': site_Bag, 
                         'Bbg': site_Bbg, 'OM': site_OM}
            
            rslr = SLR[iid] - site_uplift[iid]
            Cs0 = np.array([site_TSM[iid], site_sal[iid]], dtype=np.float64, 
                           order='F')
            forcings = {'U10': U10[iid], 'Tair': Tair[iid], 'h0': h0[iid], 
                        'U0': U0[iid], 'Hwav0': Hwav0[iid], 'Twav': Twav[iid], 
                        'Cs0': Cs0, 'rslr': rslr, 'trng': site_trng[iid], 
                        'mhws': site_mhws[iid], 'mhwn': site_mhwn[iid]}
            
            input_data = {'coord': coords, 'state': tai_state, 
                          'forcings': forcings, 'namelist': namelist, 
                          'id': site_id}
            tai_state, __, __ = cpl.run_tai_maces(input_data, models, True)
            
            # then do the formal run
            input_data = {'coord': coords, 'state': tai_state, 
                          'forcings': forcings, 'namelist': namelist, 
                          'id': site_id}
            __, uhydro_out, ecogeom_out = cpl.run_tai_maces(input_data, \
                models, False)
            
        except AssertionError as errstr:
            # print error message and exit the program
            print("Model spinup stops due to that", errstr)
            taihydro.finalizehydromod()
            sys.exit()
            
        # deallocate
        taihydro.finalizehydromod()
        
        # gather data to master
        uhydro_out_gather = {}
        ecogeom_out_gather = {}
        if master_process:
            sids = np.zeros(numprocs, dtype=np.int32)
            for okey in uhydro_out:
                oshape = list(np.shape(uhydro_out[okey]))
                oshape[0] = numprocs
                oshape = tuple(oshape)
                dtype = uhydro_out[okey].dtype
                uhydro_out_gather[okey] = np.zeros(oshape,dtype=dtype)
            for okey in ecogeom_out:
                oshape = list(np.shape(ecogeom_out[okey]))
                oshape[0] = numprocs
                oshape = tuple(oshape)
                dtype = ecogeom_out[okey].dtype
                ecogeom_out_gather[okey] = np.zeros(oshape,dtype=dtype)
        else:
            sids = None
            for okey in uhydro_out:
                uhydro_out_gather[okey] = None
            for okey in ecogeom_out:
                ecogeom_out_gather[okey] = None
        comm.Gather(np.array([iid],dtype=np.int32), sids, root=0)
        for okey in uhydro_out:
            counts = ()
            dspls = ()
            nsize = np.size(uhydro_out[okey])
            for jj in range(numprocs):
                counts = counts + (nsize,)
                dspls = dspls + (jj*nsize,)
            mpi_dtype = utils.get_mpi_dtype(uhydro_out[okey].dtype)
            sendbuf = [uhydro_out[okey], nsize]
            recvbuf = [uhydro_out_gather[okey], counts, dspls, mpi_dtype]
            comm.Gatherv(sendbuf, recvbuf, root=0)
        for okey in ecogeom_out:
            counts = ()
            dspls = ()
            nsize = np.size(ecogeom_out[okey])
            for jj in range(numprocs):
                counts = counts + (nsize,)
                dspls = dspls + (jj*nsize,)
            mpi_dtype = utils.get_mpi_dtype(ecogeom_out[okey].dtype)
            sendbuf = [ecogeom_out[okey], nsize]
            recvbuf = [ecogeom_out_gather[okey], counts, dspls, mpi_dtype]
            comm.Gatherv(sendbuf, recvbuf, root=0)
            
        # archive outputs
        if master_process:
            if ii==0:
                to_create = True
            else:
                to_create = False
            utils.write_hydro_outputs(namelist['FILE_HYDRO'], site_ids, 
                                      sids, namelist['HYDRO_TSTEP'], 
                                      uhydro_out_gather, to_create)
            utils.write_ecogeom_outputs(namelist['FILE_ECOGEOM'], site_ids, 
                                        sids, namelist['ECOGEOM_TSTEP'], 
                                        ecogeom_out_gather, to_create)
