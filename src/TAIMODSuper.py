#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  5 21:53:49 2019

Base class for coastal wetland models (geomorphology + biogeochemistry)

@author: Zeli Tan
"""

"""
Platform plant function type
0 => uninhabitable barriers (such as roads, dyks, harbors and etc)
1 => tidal flats
2 => salt marshes
3 => brackish marshes
4 => freshwater marshes
5 => mangroves
6 => needleleaf evergreen tree
7 => needleleaf deciduous tree
7 => broadleaf evergreen tree
8 => broadleaf deciduous tree
"""

import numpy as np
from abc import ABCMeta, abstractmethod

class MACMODSuper(object):
    """Abstract base class for TAI mineral accretion models.

    Attributes:
        m_params : model parameters

    """
    
    m_params = {}
    
    __metaclass__ = ABCMeta
    
    @abstractmethod 
    def mineral_suspend(self, inputs):
        """"Calculate mineral suspension rate.
        Arguments:
            inputs : driving data for mineral suspension calculation
            inputs['x']   : platform coordinate (m)
            inputs['zh']  : platform elevation relative to MSL (m)
            inputs['pft'] : platform vegetation cover pft
            inputs['Css'] : suspended sediment concentration (g m-3)
            inputs['tau'] : bottom shear stress (Pa)
            inputs['Bag'] : aboveground biomass (gC m-2)
            inputs['U']   : water flow velocity (m s-1)
            inputs['h']   : water depth (m)
            inputs['TR']  : tidal range (m)
            inputs['d50'] : sediment median diameter (m)
        Returns: mineral suspension rate (g m-2 s-1)
        """
        pass
    
    @abstractmethod
    def mineral_deposition(self, inputs):
        """"Calculate mineral deposition rate.
        Arguments:
            inputs : driving data for mineral deposition calculation
            inputs['x']   : platform coordinate (m)
            inputs['zh']  : platform elevation relative to MSL (m)
            inputs['pft'] : platform vegetation cover pft
            inputs['Css'] : suspended sediment concentration (g m-3)
            inputs['tau'] : bottom shear stress (Pa)
            inputs['Bag'] : aboveground biomass (gC m-2)
            inputs['U']   : water flow velocity (m s-1)
            inputs['h']   : water depth (m)
            inputs['TR']  : tidal range (m)
            inputs['d50'] : sediment median diameter (m)
        Returns: mineral deposition rate (g m-2 s-1)
        """
        pass
    
class OMACMODSuper(object):
    """Abstract base class for TAI organic matter accretion models.

    Attributes:
        m_params : model parameters

    """
    
    m_params = {}
    
    __metaclass__ = ABCMeta    
    
    @abstractmethod 
    def organic_accretion(self, inputs):
        """"Calculate organic matter accretion rate.
        Arguments:
            inputs : driving data for OM accretion calculation
            inputs['x']     : platform coordinate (m)
            inputs['zh']    : platform elevation relative to MSL (m)
            inputs['pft']   : platform vegetation cover pft
            inputs['h']     : water depth (m)
            inputs['TR']    : tidal range (m)
            inputs['Tair']  : annual mean air temperature (K)
            inputs['month'] : month
            inputs['day']   : day
        Returns: organic matter accretion rate (g m-2 s-1)
        """
        pass
    
class WINDEROMODSuper(object):
    """Abstract base class for TAI storm surge erosion models.

    Attributes:
        m_params : model parameters

    """
    
    m_params = {}
    
    __metaclass__ = ABCMeta

    @abstractmethod 
    def wind_erosion(self, inputs):
        """"Calculate storm surge erosion rate. This should be used to get the
        new values of grid cell coordinate and elevation.
        Arguments:
            inputs : driving data for storm surge erosion calculation
            inputs['x']    : platform coordinate (m)
            inputs['zh']   : platform elevation relative to MSL (m)
            inputs['pft']  : platform vegetation cover pft
            inputs['Hwav'] : significant wave height (m)
        Returns: storm surge erosion rate (m s-1)
        """
        pass
    
class LNDMGMODSuper(object):
    """Abstract base class for TAI landward migration models.

    Attributes:
        m_params : model parameters

    """
    
    m_params = {}
    
    __metaclass__ = ABCMeta
    
    @abstractmethod 
    def landward_migration(self, inputs):
        """"Calculate coastal wetland landward migration at the end of each year.
        Arguments:
            inputs : driving data for landward migration calculation
            inputs['x']      : platform coordinate (m)
            inputs['zh']     : platform elevation relative to MSL (m)
            inputs['pft']    : platform vegetation cover pft
            inputs['h_hist'] : daily water depth (m) of the last five years
            inputs['Cj_hist']: daily water salinity (PSU) of the last five years
        Returns: the new pft on the coastal platform
        """
        pass

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
    """
    N = 1000
    zh = np.zeros(N, dtype=np.float64, order='F')
    return zh