#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 30 22:03:14 2019

Utility functions for the MACES

@author: Zeli Tan
"""

import netCDF4
import numpy as np
import xml.etree.ElementTree as ET
from scipy import constants
from datetime import date
from mpi4py import MPI

NTOPSEG = 17

G = constants.g
Karman = 0.42
Roul = 1028.0
visc = 1e-6     # kinematic viscosity of seawater (m2/s)
TOL = 1e-6      # tolerance for near-zero state variable

def get_date_from_julian(julian):
    """Get date from Julian day number
    Arguments:
        julian : Julian day number
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
    Arguments:
        year, month, day
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

def get_spinup_stop_date(date0, tstep, nstep):
    """Get the stop date object of spinup
    Arguments:
        date0 : start date object
        tstep : time step string
        nstep : number of time steps of spin-up
    Returns : date object
    """
    if tstep=='ndays':
        jdn = get_julian_from_date(date0.year, date0.month, date0.day)
        year, month, day = get_date_from_julian(jdn+nstep)
    elif tstep=='nmonths':
        day = date0.day
        month = np.mod(date0.month+nstep-1,12) + 1
        year = date0.year + int((date0.month+nstep-1)/12)
    elif tstep=='nyears':
        day = date0.day
        month = date0.month
        year = date0.year + nstep
    return date(year, month, day)
        
def construct_tai_platform(diva_segments, xRes, nmax):
    """Construct the MACES TAI platform.
       DIVA elevations are fixed at -12.5, -8.5, -5.5, -4.5, -3.5, -2.5, -1.5,
       -0.5, 0, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 8.5, 12.5, 16.5 msl.
    Arguments:
        diva_segments : DIVA segment length (km)
        xRes : reference node cell length (m)
        nmax : maximum cell number in a segment
    Returns :
        x_tai   : platform grid coordinate (m)
        zh_tai  : platform grid elevation (m)
    """
    zhs = [-12.5, -8.5, -5.5, -4.5, -3.5, -2.5, -1.5, -0.5, 0, 0.5, 1.5, 
           2.5, 3.5, 4.5, 5.5, 8.5, 12.5, 16.5]
    assert len(zhs)-1 == len(diva_segments), \
        "DIVA segments do not match with elevation nodes"
    Nx = 0
    for length in diva_segments:
        nnode = int( 1e3 * length / xRes )
        if nnode>0:
            Nx = Nx + min( max(nnode,2), nmax )
    Nx = Nx + 1     # the end node
    x_tai = np.zeros(Nx, dtype=np.float64, order='F')
    zh_tai = np.zeros(Nx, dtype=np.float64, order='F')
    indx = 0
    x0 = 0.0
    for ii, length in enumerate(diva_segments):
        nnode = int( 1e3 * length / xRes )
        if nnode>0:
            nnode = min( max( nnode, 2 ), nmax )
            x_tai[indx:indx+nnode] = x0 + 1e3*length*np.arange(nnode)/nnode
            zh_tai[indx:indx+nnode] = zhs[ii] + (zhs[ii+1]-zhs[ii])* \
                np.arange(nnode)/nnode
            indx = indx + nnode
            x0 = x0 + 1e3*length
            # the end node
            x_tai[-1] = x0
            zh_tai[-1] = zhs[ii+1]
    return x_tai, zh_tai

def construct_platform_pft(segments, pfts, x_tai):
    """Construct the pft on the MACES platform.
    Arguments:
        segments : segment length (km)
        pfts : pft on each segment 
        x_tai : coordinate of MACES platform nodes (m)
    Returns : The pft on the MACES platform
    """
    nseg = len(segments)
    pft_tai = np.zeros_like(x_tai, dtype=np.int8, order='F')
    for ii in range(nseg):
        x0 = 1e3 * np.sum(segments[:ii])
        x1 = x0 + 1e3 * segments[ii]
        indice = np.logical_and(x_tai>=x0, x_tai<x1)
        pft_tai[indice] = pfts[ii]
    indice = x_tai>=1e3*np.sum(segments)
    pft_tai[indice] = 0
    return pft_tai
            
def get_refshore_coordinate(x, zh):
    """Get the coordinate of the shore at msl.
    Arguments:
        x : longitudinal coordinate (m)
        zh : platform surface elevation (msl)
    Returns : the coordinate of shore at msl (km)
    """
    indx = np.argmin(np.abs(zh))
    return 1e-3 * x[indx]

def get_platform_slope(x, zh):
    """Get the platform slope.
    Arguments:
        x : longitudinal coordinate (m)
        zh : platform surface elevation (msl)
    Returns : the platform cell slope (m/m)
    """
    Nx = np.size(x)
    dzh = np.zeros(Nx, dtype=np.float64, order='F')
    for ii in range(Nx):
        if ii>0 and ii<Nx-1:
            dzh[ii] = (zh[ii+1] - zh[ii-1]) / (x[ii+1] - x[ii-1])
        elif ii==0:
            dzh[ii] = (zh[ii+1] - zh[ii]) / (x[ii+1] - x[ii])
        else:
            dzh[ii] = (zh[ii] - zh[ii-1]) / (x[ii] - x[ii-1])
    return dzh

def update_platform_elev(zh, Esed, Dsed, Lbed, DepOM, rhoSed, rhoOM, 
                         porSed, rslr, dt):
    """Get the platform slope.
    Arguments:
        zh : platform surface elevation (msl)
        Esed : sediment erosion rate (kg/m2/s)
        Dsed : sediment deposition rate (kg/m2/s)
        Lbed : bedload rate (kg/m2/s)
        DepOM : OM deposition rate (kg/m2/s)
        rhoSed : sediment density (kg/m3)
        rhoOM : OM density (kg/m3)
        porSed : sediment porosity
        rslr : relative sea level rise (mm/yr)
        dt : time step (s)
    Returns : new platform surface elevation (msl)
    """
    return zh + ((Dsed/rhoSed + Lbed/rhoSed + DepOM/rhoOM - Esed/rhoSed)/ \
                 (1.0-porSed) - rslr/3.1536e10) * dt

#def estimate_Hwav_seaward(U10, h0):
#    """Estimate the significant wave height at the seaward side.
#    Arguments:
#        U10 : 10-m wind speed (m/s)
#        h0 : water depth at the seaward side (m)
#    Returns : significant wave height (m)
#    """
#    Hwav_ocean = 0.27 * U10**2 / G
#    Hwav0 = np.sinh(Karman*h0) / np.sinh(Karman*30.0) * Hwav_ocean
#    return Hwav0

def parseXML_namelist(xmlfile):
    """Read MACES settings from a xml file.
    https://docs.python.org/2/library/xml.etree.elementtree.html
    Arguments:
        xmlfile : the file name string
    Returns : a model setting dictionary
    """
    tree = ET.parse(xmlfile)
    root = tree.getroot()
    namelist = {}
    keys_str = []
    for entry in root.findall('./group/entry'):
        key = entry.get('id')
        dtype = entry.find('type').text
        if dtype=='char':
            keys_str.append(key)
        if 'value' in entry.keys():
            # single value setting
            if dtype=='integer':
                value = int(entry.get('value'))
            elif dtype=='real':
                value = float(entry.get('value'))
            elif dtype=='logical':
                if entry.get('value')=='TRUE':
                    value = True
                else:
                    value = False
            else:
                value = entry.get('value')
            namelist[key] = value
        else:
            # list value setting
            values = []
            for value in entry.findall('values/value'):
                if dtype=='integer':
                    values.append(int(value.text))
                elif dtype=='real':
                    values.append(float(value.text))
                elif dtype=='logical':
                    if value.text=='TRUE':
                        values.append(True)
                    else:
                        values.append(False)
                else:
                    values.append(value.text)
            namelist[key] = np.array(values, dtype=np.float64, order='F')
    # fill environment variable values
    for key in keys_str:
        keys_str_iter = list(set(keys_str)-set([key]))
        for okey in keys_str_iter:
            string = namelist[okey]
            if string.find('$'+key)>=0:
                new_string = string.replace('$'+key, namelist[key])
                namelist[okey] = new_string
    return namelist
    
def parseXML_params(xmlfile, model):
    """Read model parameters from a xml file.
    https://docs.python.org/2/library/xml.etree.elementtree.html
    Arguments:
        xmlfile : the file name string
        model : model name string
    Returns : a model parameter dictionary
    """
    tree = ET.parse(xmlfile)
    root = tree.getroot()
    param_dict = {}
    findstr = "./group/[@id='" + model + "']/entry"
    for entry in root.findall(findstr):
        key = entry.get('id')
        dtype = entry.find('type').text
        if 'value' in entry.keys():
            # single value setting
            if dtype=='integer':
                value = int(entry.get('value'))
            elif dtype=='real':
                value = float(entry.get('value'))
            elif dtype=='logical':
                value = bool(entry.get('value'))
            else:
                value = entry.get('value')
            param_dict[key] = value
        else:
            # list value setting
            values = []
            for value in entry.findall('values/value'):
                if dtype=='integer':
                    values.append(int(value.text))
                elif dtype=='real':
                    values.append(float(value.text))
                elif dtype=='logical':
                    values.append(bool(value.text))
                else:
                    values.append(value.text)
            param_dict[key] = np.array(values, dtype=np.float64, order='F')
    return param_dict

def parseXML_hydro_params(xmlfile):
    """Read model parameters from a xml file.
    https://docs.python.org/2/library/xml.etree.elementtree.html
    Arguments:
        xmlfile : the file name string
    Returns : a model parameter dictionary
    """
    tree = ET.parse(xmlfile)
    root = tree.getroot()
    param_dict = {}
    for entry in root.findall("./entry"):
        key = entry.get('id')
        dtype = entry.find('type').text
        if 'value' in entry.keys():
            # single value setting
            if dtype=='integer':
                value = int(entry.get('value'))
            elif dtype=='real':
                value = float(entry.get('value'))
            elif dtype=='logical':
                value = bool(entry.get('value'))
            else:
                value = entry.get('value')
            param_dict[key] = value
        else:
            # list value setting
            values = []
            for value in entry.findall('values/value'):
                if dtype=='integer':
                    values.append(int(value.text))
                elif dtype=='real':
                    values.append(float(value.text))
                elif dtype=='logical':
                    values.append(bool(value.text))
                else:
                    values.append(value.text)
            param_dict[key] = np.array(values, dtype=np.float64, order='F')
    return param_dict

def get_forcing_index(t, tstep, ntstep):
    """Get the time index of forcing.
    Arguments:
        t : time in seconds
        tstep : time step string
        ntstep : number of time steps
    Returns : forcing time index
    """
    if tstep=='minute':
        index = int( t/(60*ntstep) )
    elif tstep=='hour':
        index = int( t/(3600*ntstep) )
    else:
        index = -1
    return index

def read_force_data(filename, varname, date0, date1, ntstep, 
                    tstep, id_range):
    """Read forcing data from a nc file.
    Arguments:
        filename : forcing data file
        varname : variable name
        date0 : date object of the first record
        date1 : date object of the last record
        ntstep : number of record time step
        tstep : string of record time step
        id_range : [first site id, last side id]
    Returns : forcing data array
    """
    try:
        nc = netCDF4.Dataset(filename, 'r')
        dateint = int(nc.variables['date'][:])
        year = int(dateint/1e4)
        month = int((dateint-1e4*year)/1e2)
        day = int(dateint - 1e4*year - 1e2*month)
        refdate = date(year, month, day)
        id0 = id_range[0]
        id1 = id_range[1] + 1
        nday = (date1 - date0).days
        day0 = (date0 - refdate).days
        nyear = date1.year - date0.year
        year0 = date0.year - refdate.year
        if tstep=='hour':
            nstart = int( 24*day0/ntstep )
            ntime = max( int( 24*nday/ntstep ), 1 )
        elif tstep=='minute':
            nstart = int( 24*60*day0/ntstep )
            ntime = max( int( 24*60*nday/ntstep ), 1 )
        elif tstep=='year':
            nstart = int( year0/ntstep )
            ntime = max( int( nyear/ntstep ), 1 )
        data = np.array(nc.variables[varname][:][id0:id1,nstart:nstart+ntime])
    finally:
        nc.close()
    return data

def get_shr_output_index(t, tstep):
    """Get the time index of short term outputs.
    Arguments:
        t : time in seconds
        tstep : time step string
    Returns : time index
    """
    if tstep=='minute':
        index = int( t/60 )
    elif tstep=='hour':
        index = int( t/3600 )
    else:
        index = -1
    return index

def get_shr_output_num(date0, date1, tstep):
    """Get the total number of short term outputs.
    Arguments:
        date0 : run start date
        date1 : run stop date
        tstep : time step string
    Returns : output record number
    """
    nday = (date1 - date0).days
    nhour = 24 * nday
    nmin = 60 * nhour
    if tstep=='hour':
        ntime = nhour
    elif tstep=='minute':
        ntime = nmin
    else:
        ntime = -1
    return ntime

def get_lng_output_index(date0, date_cur, tstep):
    """Get the time index of long term outputs.
    Arguments:
        date0 : simulation start date
        date_cur : simulation current date
        tstep : time step string
    Returns : time index
    """
    if tstep=='day':
        index = (date_cur - date0).days
    elif tstep=='month':
        index = 12*(date_cur.year-date0.year) + (date_cur.month-date0.month)
    elif tstep=='year':
        index = date_cur.year - date0.year
    elif tstep=='decade':
        index = int( (date_cur.year-date0.year)/10 )
    elif tstep=='century':
        index = int( (date_cur.year-date0.year)/100 )
    else:
        index = -1
    return index

def get_lng_output_num(date0, date1, tstep):
    """Get the total number of long term outputs.
    Arguments:
        date0 : run start date
        date1 : run stop date
        tstep : time step string
    Returns : output record number
    """
    nday = (date1 - date0).days
    nmonth = 12*(date1.year-date0.year) + (date1.month-date0.month)
    nyear = date1.year - date0.year
    if tstep=='day':
        ntime = nday
    elif tstep=='month':
        ntime = nmonth
    elif tstep=='year':
        ntime = nyear
    elif tstep=='decade':
        ntime = int(nyear/10)
    elif tstep=='century':
        ntime = int(nyear/100)
    else:
        ntime = -1
    return ntime

def get_mpi_dtype(dtype):
    """Get the mpi data type corresponding to numpy data type.
    Arguments:
        dtype : numpy array data type
    Returns : mpi data type
    """
    if dtype==np.dtype('float64'):
        return MPI.DOUBLE
    elif dtype==np.dtype('float32'):
        return MPI.FLOAT
    elif dtype==np.dtype('int32'):
        return MPI.INT
    elif dtype==np.dtype('int8'):
        return MPI.BYTE

def write_hydro_outputs(filename, all_ids, sids, tstep, uhydro_out, 
                        to_create):
    """Write model hydrodynamics outputs into a nc file.
    Arguments:
        filename : output file name
        all_ids : all site ids
        sids : site index range of the current iteration
        tstep : time step type string
        uhydro_out  : hydrodynamic model outputs
        to_create : True if create the new file otherwise False
    Returns : 
    """
    nid = len(all_ids)
    nt = np.shape(uhydro_out['h'])[1]
    nx = np.shape(uhydro_out['h'])[2]
    # create output file if needed
    if to_create:
        try:
            nc = netCDF4.Dataset(filename, 'w', format='NETCDF4_CLASSIC')
            nc.history = 'MACES simulated ' + tstep + ' hydrodynamics'
            nc.contact = r'Please contact zeli.tan@pnnl.gov for more information'
            nc.createDimension('site', None)
            nc.createDimension('time', nt)
            nc.createDimension('x', nx)
            # create and write variables
            site_var = nc.createVariable('site', 'i4', ('site',))
            site_var.long_name = r'site DIVA id'
            site_var[:] = all_ids
            x_var = nc.createVariable('x', 'f4', ('site','x',))
            x_var.long_name = r'platform transect coordinate'
            x_var.units = 'm'
            x_var[:] = 1e20*np.ones((nid,nx),dtype=np.float32)
            h_var = nc.createVariable('h', 'f4', ('site','time','x',), 
                                      fill_value=1e20)
            h_var.long_name = r'water depth'
            h_var.units = 'm'
            h_var[:] = 1e20*np.ones((nid,nt,nx),dtype=np.float32)
            U_var = nc.createVariable('U', 'f4', ('site','time','x',), 
                                      fill_value=1e20)
            U_var.long_name = r'tide signed flow velocity'
            U_var.units = 'm/s'
            U_var[:] = 1e20*np.ones((nid,nt,nx),dtype=np.float32)
            Hwav_var = nc.createVariable('Hwav', 'f4', ('site','time','x',), 
                                         fill_value=1e20)
            Hwav_var.long_name = r'significant wave height'
            Hwav_var.units = 'm'
            Hwav_var[:] = 1e20*np.ones((nid,nt,nx),dtype=np.float32)
            Uwav_var = nc.createVariable('Uwav', 'f4', ('site','time','x',), 
                                         fill_value=1e20)
            Uwav_var.long_name = r'wave velocity'
            Uwav_var.units = 'm/s'
            Uwav_var[:] = 1e20*np.ones((nid,nt,nx),dtype=np.float32)            
            tau_var = nc.createVariable('tau', 'f4', ('site','time','x',), 
                                        fill_value=1e20)
            tau_var.long_name = r'bottom shear stress'
            tau_var.units = 'Pa'
            tau_var[:] = 1e20*np.ones((nid,nt,nx),dtype=np.float32)
            Css_var = nc.createVariable('TSM', 'f4', ('site','time','x',), 
                                        fill_value=1e20)
            Css_var.long_name = r'suspended sediment concentration'
            Css_var.units = 'kg/m3'
            Css_var[:] = 1e20*np.ones((nid,nt,nx),dtype=np.float32)
            Cj_var = nc.createVariable('sal', 'f4', ('site','time','x',), 
                                       fill_value=1e20)
            Cj_var.long_name = r'water salinity'
            Cj_var.units = 'PSU'
            Cj_var[:] = 1e20*np.ones((nid,nt,nx),dtype=np.float32)
        finally:
            nc.close()
    # write data only
    try:
        nc = netCDF4.Dataset(filename, 'r+')
        for ii, iid in enumerate(sids):
            x_var = nc.variables['x']
            x_var[iid] = uhydro_out['x'][ii]
            h_var = nc.variables['h']
            h_var[iid] = uhydro_out['h'][ii]
            U_var = nc.variables['U']
            U_var[iid] = uhydro_out['U'][ii]
            Hwav_var = nc.variables['Hwav']
            Hwav_var[iid] = uhydro_out['Hwav'][ii]
            Uwav_var = nc.variables['Uwav']
            Uwav_var[iid] = uhydro_out['Uwav'][ii]
            tau_var = nc.variables['tau']
            tau_var[iid] = uhydro_out['tau'][ii]
            Css_var = nc.variables['TSM']
            Css_var[iid] = uhydro_out['Css'][ii]
            Cj_var = nc.variables['sal']
            Cj_var[iid] = uhydro_out['Cj'][ii]
    finally:
        nc.close()
        
def write_ecogeom_outputs(filename, all_ids, sids, tstep, ecogeom_out, 
                          to_create):
    """Write model outputs into a nc file.
    Arguments:
        filename : output file name
        all_ids : all site ids
        sids : site index range of the current iteration
        tstep : time step type string
        ecogeom_out  : ecogeomorphology model outputs
        to_create : True if create the new file otherwise False
    Returns : 
    """
    nid = len(all_ids)
    nt = np.shape(ecogeom_out['OM'])[1]
    nx = np.shape(ecogeom_out['OM'])[2]
    npool = np.shape(ecogeom_out['OM'])[3]
    # create output file if needed
    if to_create:
        try:
            nc = netCDF4.Dataset(filename, 'w', format='NETCDF4_CLASSIC')
            nc.history = 'MACES simulated ' + tstep + ' eco-geomorphology'
            nc.contact = r'Please contact zeli.tan@pnnl.gov for more information'
            nc.createDimension('site', None)
            nc.createDimension('time', nt)
            nc.createDimension('x', nx)
            nc.createDimension('pool',npool)
            # create and write variables
            site_var = nc.createVariable('site', 'i4', ('site',))
            site_var.long_name = r'site DIVA id'
            site_var[:] = all_ids
            x_var = nc.createVariable('x', 'f4', ('site','x',))
            x_var.long_name = r'platform transect coordinate'
            x_var.units = 'm'
            x_var[:] = 1e20*np.ones((nid,nx),dtype=np.float32)
            pft_var = nc.createVariable('pft', 'i1', ('site','time','x',), 
                                        fill_value=-1)
            pft_var.long_name = r'platform plant function type'
            pft_var.units = '0 to 8'
            pft_var[:] = -1*np.ones((nid,nt,nx),dtype=np.int8)
            zh_var = nc.createVariable('zh', 'f4', ('site','time','x',), 
                                       fill_value=1e20)
            zh_var.long_name = r'platform surface elevation'
            zh_var.units = 'msl'
            zh_var[:] = 1e20*np.ones((nid,nt,nx),dtype=np.float32)
            Esed_var = nc.createVariable('Esed', 'f4', ('site','time','x',), 
                                         fill_value=1e20)
            Esed_var.long_name = r'sediment erosion rate'
            Esed_var.units = 'kg/m2/s'
            Esed_var[:] = 1e20*np.ones((nid,nt,nx),dtype=np.float32)
            Dsed_var = nc.createVariable('Dsed', 'f4', ('site','time','x',), 
                                         fill_value=1e20)
            Dsed_var.long_name = r'suspended sediment deposition rate'
            Dsed_var.units = 'kg/m2/s'
            Dsed_var[:] = 1e20*np.ones((nid,nt,nx),dtype=np.float32)
            Lbed_var = nc.createVariable('Lbed', 'f4', ('site','time','x',), 
                                         fill_value=1e20)
            Lbed_var.long_name = r'sand bed load rate'
            Lbed_var.units = 'kg/m2/s'
            Lbed_var[:] = 1e20*np.ones((nid,nt,nx),dtype=np.float32)
            DepOM_var = nc.createVariable('DepOM', 'f4', ('site','time','x',), 
                                          fill_value=1e20)
            DepOM_var.long_name = r'Organic matter deposition rate'
            DepOM_var.units = 'kg/m2/s'
            DepOM_var[:] = 1e20*np.ones((nid,nt,nx),dtype=np.float32)
            Bag_var = nc.createVariable('Bag', 'f4', ('site','time','x',), 
                                        fill_value=1e20)
            Bag_var.long_name = r'platform aboveground biomass'
            Bag_var.units = 'kg/m2'
            Bag_var[:] = 1e20*np.ones((nid,nt,nx),dtype=np.float32)
            Bbg_var = nc.createVariable('Bbg', 'f4', ('site','time','x',), 
                                        fill_value=1e20)
            Bbg_var.long_name = r'platform belowground biomass'
            Bbg_var.units = 'kg/m2'
            Bbg_var[:] = 1e20*np.ones((nid,nt,nx),dtype=np.float32)
            OM_var = nc.createVariable('OM', 'f4', ('site','time','x','pool',), 
                                       fill_value=1e20)
            OM_var.long_name = r'platform column-integrated soil organic matter'
            OM_var.units = 'kg/m2'
            OM_var[:] = 1e20*np.ones((nid,nt,nx,npool),dtype=np.float32)
        finally:
            nc.close()
    # write data only
    try:
        nc = netCDF4.Dataset(filename, 'r+')
        for ii, iid in enumerate(sids):
            x_var = nc.variables['x']
            x_var[iid] = ecogeom_out['x'][ii]
            pft_var = nc.variables['pft']
            pft_var[iid] = ecogeom_out['pft'][ii]
            zh_var = nc.variables['zh']
            zh_var[iid] = ecogeom_out['zh'][ii]
            Esed_var = nc.variables['Esed']
            Esed_var[iid] = ecogeom_out['Esed'][ii]
            Dsed_var = nc.variables['Dsed']
            Dsed_var[iid] = ecogeom_out['Dsed'][ii]
            Lbed_var = nc.variables['Lbed']
            Lbed_var[iid] = ecogeom_out['Lbed'][ii]
            DepOM_var = nc.variables['DepOM']
            DepOM_var[iid] = ecogeom_out['DepOM'][ii]
            Bag_var = nc.variables['Bag']
            Bag_var[iid] = ecogeom_out['Bag'][ii]
            Bbg_var = nc.variables['Bbg']
            Bbg_var[iid] = ecogeom_out['Bbg'][ii]
            OM_var = nc.variables['OM']
            OM_var[iid] = ecogeom_out['OM'][ii]
    finally:
        nc.close()