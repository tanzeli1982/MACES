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
import maces_utilities as utils
from abc import ABCMeta, abstractmethod

###############################################################################
class MACMODSuper(object):
    """Abstract base class for TAI mineral accretion models.

    Attributes:
        m_params : model parameters

    """
    
    m_params = {}
    
    __metaclass__ = ABCMeta
    
    @abstractmethod 
    def mineral_suspension(self, inputs):
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
        Returns: mineral suspension rate (kg m-2 s-1)
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
        Returns: mineral deposition rate (kg m-2 s-1)
        """
        pass
    
    def settling_velocity(self, tau, d50, Rous):
        """"Calculate effective sediment settling velocity (Mudd et al., 2010).
        Arguments:
            tau : bottom shear stress (Pa)
            d50 : sediment median diameter (m)
            Rous : sediment density (kg/m3)
        Returns: sediment settling velocity (m s-1)
        """
        # parameters for cohesive sediment (clay and silt)
        A = 38.0
        F = 3.55
        m = 1.2
        S = Rous / utils.Roul
        nv = utils.visc
        G = utils.G
        ws = (( np.sqrt(0.25*(A/F)**(2/m)+(4./3.*d50**3*G*(S-1)/F/nv**2)**(1/m)) \
              - 0.5*(A/F)**(1/m))**m) * nv / d50
        ustar = np.sqrt(tau/utils.Roul)
        wup = utils.Karman * ustar
        return ws - wup

###############################################################################    
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
        Returns: 
            omac : organic matter accretion rate (kg m-2 s-1)
            Bag  : aboveground biomass (kg m-2)
        """
        pass

###############################################################################    
class WINDEROMODSuper(object):
    """Abstract base class of models for TAI storm surge erosion at wetland edges.

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

###############################################################################    
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