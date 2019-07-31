#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 30 22:03:14 2019

Utility functions for the MACES

@author: Zeli Tan
"""

import numpy as np

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
        
def construct_platform_elev(diva_topo):
    """Construct the TAI platform elevation profile
    Arguments:
        diva_topo : DIVA topography input
    Returns : new platform elevation
    """
    ntopo = len(diva_topo)
    N = 
    zh = np.zeros(N, dtype=np.float64)
    return zh