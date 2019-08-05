#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 30 22:03:14 2019

Utility functions for the MACES

@author: Zeli Tan
"""

import numpy as np
import pandas as pd
from scipy import constants

G = constants.g
Karman = 0.42

#    def biomass(self, TR):
#        """update plant-induced aboveground biomass
#        Arguments:
#            TR : tidal range
#        """
#        self.m_params['marsh']['aa'] = B_marsh['aa']
#        self.m_params['marsh']['bb'] = B_marsh['bb']
#        self.m_params['marsh']['cc'] = B_marsh['cc']
#        self.m_params['mangrove']['aa'] = B_mangrove['aa']
#        self.m_params['mangrove']['bb'] = B_mangrove['bb']
#        self.m_params['mangrove']['cc'] = B_mangrove['cc']
#        hMHT = TR
#        DMHT = hMHT - self.zh_arr
#        pft = self.m_pft[ii]
#        if pft in ['marsh','mangrove']:
#            aa = self.m_params[pft]['aa']
#            alpha_a = self.m_params[pft]['alpha_a']
#            beta_a = self.m_params[pft]['beta_a']
#            alpha_d = self.m_params[pft]['alpha_d']
#            beta_d = self.m_params[pft]['beta_d']
#            cD0 = self.m_params[pft]['cD0']
#            ScD = self.m_params[pft]['ScD']
#        else:
#            aa, bb, cc = [0.0, 0.0, 0.0]
#            alpha_a, beta_a, alpha_d, beta_d = [0.0, 0.0, 0.0, 0.0]
#            cD0, ScD = [0.0, 0.0]
#        self.m_Bag[ii] = aa*DMHT[ii] + bb*(DMHT[ii])**2 + cc
        
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
    