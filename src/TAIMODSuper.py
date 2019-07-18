#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  5 21:53:49 2019

Base class for coastal wetland models (geomorphology + biogeochemistry)

@author: Zeli Tan
"""

import numpy as np
from abc import ABCMeta, abstractmethod

class TAIMODSuper(object):
    """Abstract base class for TAI models.

    Attributes:
        sediment_arr[:]: sediment height
        biomass_arr[:]: plant biomass

    """
    
    sediment_arr = None
    biomass_arr = None
    
    __metaclass__ = ABCMeta
    
    @abstractmethod 
    def mineral_accretion(self, sed, shear):
        """"Calculate mineral accretion rate.
        Parameters:
            sed : suspended sediment concentration (g m-3)
            shear : water flow shear stress (Pa)
        Returns: mineral accretion rate (g m-2 s-1)
        """
        pass
    
    @abstractmethod 
    def organic_accretion(self, tide):
        """"Calculate organic matter accretion rate.
        Parameters:
            tide : inundation depth (m)
        Returns: organic matter accretion rate (g m-2 s-1)
        """
        pass

    @abstractmethod 
    def wind_erosion(self, wave):
        """"Calculate storm surge erosion rate.
        Parameters:
            wave : wave height
        Returns: storm surge erosion rate (m s-1)
        """
        pass  
    
    @abstractmethod 
    def landward_migration(self):
        """"Calculate coastal wetland landward migration.
        Parameters:
            
        Returns: coastal wetland landward migration
        """
        pass

    def biomass(self, TR):
        """update plant-induced aboveground biomass
        Arguments:
            TR : tidal range
        """
        self.m_params['marsh']['aa'] = B_marsh['aa']
        self.m_params['marsh']['bb'] = B_marsh['bb']
        self.m_params['marsh']['cc'] = B_marsh['cc']
        self.m_params['mangrove']['aa'] = B_mangrove['aa']
        self.m_params['mangrove']['bb'] = B_mangrove['bb']
        self.m_params['mangrove']['cc'] = B_mangrove['cc']
        hMHT = TR
        DMHT = hMHT - self.zh_arr
        pft = self.m_pft[ii]
        if pft in ['marsh','mangrove']:
            aa = self.m_params[pft]['aa']
            alpha_a = self.m_params[pft]['alpha_a']
            beta_a = self.m_params[pft]['beta_a']
            alpha_d = self.m_params[pft]['alpha_d']
            beta_d = self.m_params[pft]['beta_d']
            cD0 = self.m_params[pft]['cD0']
            ScD = self.m_params[pft]['ScD']
        else:
            aa, bb, cc = [0.0, 0.0, 0.0]
            alpha_a, beta_a, alpha_d, beta_d = [0.0, 0.0, 0.0, 0.0]
            cD0, ScD = [0.0, 0.0]
        self.m_Bag[ii] = aa*DMHT[ii] + bb*(DMHT[ii])**2 + cc
        
def construct_platform_elev(diva_topo):
    """Construct the TAI platform elevation profile
    Arguments:
        diva_topo : DIVA topography input
    """
    N = 1000
    zh = np.zeros(N, dtype=np.float64, order='F')
    return zh