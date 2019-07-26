#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 12 22:08:10 2019

Main program to run the MACES model

@author: Zeli Tan
"""

import pandas as pd
import numpy as np
import TAIMODSuper as TAIMOD
from mpi4py import MPI
from optparse import OptionParser
from scipy.io import netcdf
from TAIHydroMOD import taihydromod
from TAIHydroMOD import rungekutta4 as rk4

# model simulation constants
MAX_OF_STEP = 1800  # maximum simulation time step (s)
NXLAYER = 1000      # transect cell grid number
rk4_mode = 101      # adaptive mode
uhydro_tol = np.array([1e-6,1e-6,1e-6,1e-6,1e-6,1e-6], 
                      dtype=np.float64, order='F')  # toleratance

def run_tai_maces(params, input_data, verbose):
    # create state variables
    site_id = input_data['site_id']
    site_lon = input_data['site_lon']
    site_lat = input_data['site_lat']
    nday = np.size(input_data['wind'])
    nhour = 24 * nday
    ncell = NXLAYER
    uhydro_out = {}
    uhydro_out['h'] = np.zeros((nhour,ncell), dtype=np.float32, order='F')
    uhydro_out['U'] = np.zeros((nhour,ncell), dtype=np.float32, order='F')
    uhydro_out['Hw'] = np.zeros((nhour,ncell), dtype=np.float32, order='F')
    uhydro_out['tau'] = np.zeros((nhour,ncell), dtype=np.float32, order='F')
    uhydro_out['css'] = np.zeros((nhour,ncell), dtype=np.float32, order='F')
    uhydro_out['sal'] = np.zeros((nhour,ncell), dtype=np.float32, order='F')
    ecogeo_out = {}
    ecogeo_out['m_acc'] = np.zeros((nday,ncell), dtype=np.float32, order='F')
    ecogeo_out['om_acc'] = np.zeros((nday,ncell), dtype=np.float32, order='F')
    ecogeo_out['shore_ero'] = np.zeros((nday,ncell), dtype=np.float32, order='F')
    ecogeo_out['lnd_elev'] = np.zeros((nday,ncell), dtype=np.float32, order='F')
    ecogeo_out['abg_bio'] = np.zeros((nday,ncell), dtype=np.float32, order='F')
    ecogeo_out['bg_bio'] = np.zeros((nday,ncell), dtype=np.float32, order='F')
    ecogeo_out['coast_veg'] = np.zeros((nday,ncell), dtype=np.float32, order='F')
    # create a model object and initialize RungeKutta allocatables
    nhydro = len(uhydro_out)
    out_uhydro = np.zeros((nhydro,ncell), dtype=np.float64, order='F')
    sed_ero = np.zeros(ncell, dtype=np.float64, order='F')
    sed_dep = np.zeros(ncell, dtype=np.float64, order='F')
    om_dep = np.zeros(ncell, dtype=np.float64, order='F')
    site_zh = TAIMOD.construct_platform_elev(input_data['diva_topo'])
    taihydromod.initialize(site_x, site_zh)
    taihydromod.setmodelparameters(d50, )
    # initialize model run
    t = 0
    tf = 8.64e4 * nday
    hindx = 0
    dindx = 0
    ncount = 0
    curstep = np.array([50], dtype=np.float64, order='F')
    nextstep = np.array([MAX_OF_STEP], dtype=np.float64, order='F')
    try:
        while t < tf:
            if t>=3.6e3*hindx and hindx<nhour:
                if verbose:
                    print('Id', int(site_id), ': time step', int(hindx))
                isHourNode = True
                hindx = hindx + 1
                if np.mod(hindx,24)==0:
                    dindx = dindx + 1
            if isHourNode:
                uhydro_out['h'][hindx] = model.hydro_states[0,:]
                uhydro_out['U'][hindx] = model.U_arr
                uhydro_out['Hw'][hindx] = model.hydro_states[2,:]
                uhydro_out['tau'][hindx] = model.tau_arr
                uhydro_out['css'][hindx] = model.hydro_states[3,:]
                uhydro_out['sal'][hindx] = model.hydro_states[4,:]
                # calculate daily mean eco-geomorphology variables
            
            # simulate hydrodynamics
            taihydromod.modelsetup(zh, sed_ero, sed_dep)
            error = np.array([0], dtype=np.int32, order='F')
            rk4.rk4fehlberg(model.taihydroequations, model.m_uhydro, out_uhydro, 
                            rk4_mode, uhydro_tol, curstep, nextstep, error)
            assert error[0]==0, "Runge-Kutta iteration is more than MAXITER"
            taihydromod.callback(out_uhydro)
            # simulate eco-geomorphology
            sed_ero = 0.0
            sed_dep = 0.0
            om_dep = 0.0
            # update platform elevation
            site_zh = site_zh + (sed_dep+om_dep-sed_ero-rslr)*curstep
            # check small time step
            isHourNode = False
            if curstep<0.1:
                ncount = ncount + 1
                err_msg = 'Site' + '{:d}'.format(site_id) + ': run diverge at step' + \
                    '{:d}'.format(hindx)
                assert ncount<=100, err_msg
                nextstep = 50.0
            else:
                ncount = 0
            t = t + curstep
            curstep = nextstep
            nextstep = MAX_OF_STEP
            
    except AssertionError as error:
        # print error message
        print(error)
    finally:
        # deallocate
        taihydromod.destruct()
    # save the simulation (only short-term hydrodynamics are saved)
    
    
        
if __name__=='__main__':
    
    # use the OptionParser to parse the command line options
    parser = OptionParser()
    parser.add_option("-b", "--begin", type="int", dest="begin", default=0)
    parser.add_option("-e", "--end", type="int", dest="end", default=12148)
    parser.add_option("-v", "--verbose", type="logical", dest="verbose", default=False)
    (options, args) = parser.parse_args()
    
    site_1 = options.begin
    site_n = options.end
    verbose = options.verbose
    
    # MPI commands
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    numprocs = comm.Get_size()
    
    if rank == 0:
        master_process = True
    else:
        master_process = False
        
    # read parameters and site information
    if master_process:
        filename = 'diva_site_database.xlsx'
        df = pd.read_excel(filename, sheetname="diva", header=0, 
                           parse_cols="A,C,AC,AD")
        site_ids = np.array(df['ID'], dtype=np.int32)
        model_params = None
    else:
        site_ids = None
        model_params = None
    site_ids = comm.bcast(site_ids, root=0)
    model_params = comm.bcast(model_params, root=0)
    nsite = np.size(site_ids)
    site_n = min(site_n, nsite)
    nrun = site_n - site_1
        
    # run simulations (in each iteration, a processor tests the sensitivity of 
    # one site for all parameters)
    niter = int( np.ceil( float(nrun) / float(numprocs) ) )
    for ii in range(niter):
        iid = np.mod( ii*numprocs + rank, nrun ) + site_1
        print( "Simulate site ", site_ids[iid] )
        forcing_data = {'site_id': site_ids[iid]}
        run_tai_maces(model_params, forcing_data, verbose)
    
