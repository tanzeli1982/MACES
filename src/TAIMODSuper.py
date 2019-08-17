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
1 => tidal flats or open water
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
    
    def bed_loading(self, inputs):
        """"Calculate sand bed loading rate.
        Arguments:
            inputs : driving data for bed loading calculation
        Returns: bed loading rate (kg m-2 s-1)
        """
        x = inputs['x']
        return np.zeros_like(x, dtype=np.float64, order='F')
    
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
        return np.maximum(ws - wup, 0.0)

###############################################################################    
class OMACMODSuper(object):
    """Abstract base class for TAI organic matter accretion models.

    Attributes:
        m_params : model parameters

    """
    
    m_params = {}
    
    __metaclass__ = ABCMeta    
    
    @abstractmethod
    def organic_deposition(self, inputs):
        """"Calculate organic matter deposition rate.
        Arguments:
            inputs : driving data for OM deposition calculation
            inputs['x']     : platform coordinate (m)
            inputs['zh']    : platform elevation relative to MSL (m)
            inputs['pft']   : platform vegetation cover pft
            inputs['h']     : water depth (m)
            inputs['TR']    : tidal range (m)
            inputs['Tair']  : annual mean air temperature (K)
            inputs['month'] : month
            inputs['day']   : day
        Returns: organic matter deposition rate (kg m-2 s-1)
        """
        pass
    
    @abstractmethod
    def aboveground_biomass(self, inputs):
        """"Calculate aboveground biomass.
        Arguments:
            inputs : driving data for aboveground biomass calculation
        Returns: aboveground biomass (kg m-2)
        """
        pass
    
    def belowground_biomass(self, inputs):
        """"Calculate belowground biomass.
        Arguments:
            inputs : driving data for belowground biomass calculation
        Returns: belowground biomass (kg m-2)
        """
        phi = self.m_params['phi']  # the root:shoot quotient
        Bag = inputs['Bag']         # aboveground biomass (kg/m2)
        return phi*Bag
    
    def soilcarbon_decay(self, inputs):
        """"Calculate soil OC mineralization rate.
        Arguments:
            inputs : driving data for SOC decay rate calculation
        Returns: SOC decay rate (kg m-2 s-1) of two pools
        """
        SOM = inputs['SOM']     # soil organic matter pools (kg/m2)
        Nx = np.shape(SOM)[0]
        npool = np.shape(SOM)[1]
        return np.zeros((Nx,npool), dtype=np.float64, order='F')

###############################################################################    
class WAVEROMODSuper(object):
    """Abstract base class of models for TAI storm surge erosion at wetland edges.

    Attributes:
        m_params : model parameters

    """
    
    m_params = {}
    
    __metaclass__ = ABCMeta

    @abstractmethod 
    def wave_erosion(self, inputs):
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