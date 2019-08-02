#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 30 22:03:14 2019

Utility functions for the MACES

@author: Zeli Tan
"""

import numpy as np
import pandas as pd

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
       elevation nodes are fixed at -12.5, -8.5, -5.5, -4.5, -3.5, -2.5, -1.5,
       0, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 8.5, 12.5, 16.5 msl.
    Arguments:
        diva_segments['length'] : DIVA segment length (km)
        diva_segments['pop'] : DIVA segment population density (people/km^2)
    Returns :
        x_tai   : platform grid coordinate (m)
        zh_tai  : platform grid elevation (m)
        pop_tai : platform grid population density (people/km^2)
    """
    zhs = [-12.5, -8.5, -5.5, -4.5, -3.5, -2.5, -1.5, -0.5, 0, 0.5, 1.5, 
           2.5, 3.5, 4.5, 5.5, 8.5, 12.5, 16.5]
    assert len(zhs)-1 == len(diva_segments), \
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

def nearest_search(array, li, ri, val):
    """ Search the index of the grid that has the nearest to the target
    Arguments:
        array : a sorted array
        li : left index of the sorting range
        ri : right index of the sorting range
        val : the target value
    """
    
    while li <= ri:
        mid = int( (ri+li)/2 )
        if array[mid]==val:
            return mid
        elif array[mid]<        
        
    return -1

def construct_platform_pft(pft_grids, x_tai, x0):
    """Construct the pft on the MACES platform.
    Arguments:
        pft_grids['x'] : grid cell coordinate (m)
        pft_grids['pft'] : grid cell pft
        x_tai : coordinate of MACES platform nodes (m)
        x0 : coordinate of the grid node at msl (m)
    Returns : The pft on the MACES platform
    """
    
    