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
    """Construct the TAI platform coordinates and elevation profile.
       elevation nodes are fixed at -12.5, -8.5, -5.5, -4.5, -3.5, -2.5, -1.5,
       0, 1.5, 2.5, 3.5, 4.5, 5.5, 8.5, 12.5, 16.5 msl.
    Arguments:
        diva_segments : DIVA segment length (km)
    Returns : platform coordinate and elevation
    """
    zhs_pnts = [-12.5, -8.5, -5.5, -4.5, -3.5, -2.5, -1.5, 0, 1.5, 2.5, 3.5, 
                4.5, 5.5, 8.5, 12.5, 16.5]
    assert len(zhs_pnts)-1 == len(diva_segments), \
        "DIVA segments do not match with elevation nodes"
    Nx = 0
    xRes = 50.0
    nmax = 50
    for length in diva_segments:
        nnode = int( 1e3 * length / xRes )
        Nx = Nx + min( max(nnode,2), nmax )
    Nx = Nx + 1     # the end node
    xcor = np.zeros(Nx, dtype=np.float64)
    zh = np.zeros(Nx, dtype=np.float64)
    indx = 0
    x0 = 0.0
    for ii, length in enumerate(diva_segments):
        nnode = min( max( int(1e3*length/xRes), 2 ), nmax )
        xcor[indx:indx+nnode] = x0 + 1e3*length*np.arange(0,nnode)/nnode
        zh[indx:indx+nnode] = zhs_pnts[ii] + (zhs_pnts[ii+1]-zhs_pnts[ii])* \
            np.arange(0,nnode)/nnode
        indx = indx + nnode
        x0 = x0 + 1e3*length
    # the end node
    xcor[-1] = x0 
    zh[-1] = zhs_pnts[-1]
    return xcor, zh