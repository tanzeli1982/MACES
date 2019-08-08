#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 29 23:42:08 2019

Derived class for organic matter accretion algorithms

@author: Zeli Tan
"""

import numpy as np
from TAIMODSuper import OMACMODSuper

###############################################################################
class NULLMOD(OMACMODSuper):
    """Realization of the null organic matter accretion model.

    Attributes:
        Parameters : aa, bb, cc
    Constants:
        
    """ 
    
    # constructor
    def __init__(self, params):
        self.m_params = params
    
    def organic_accretion(self, inputs):
        """"Calculate organic matter accretion rate.
        Arguments:
            inputs : driving data for OM accretion calculation
        Returns: organic matter accretion rate (kg m-2 s-1)
        """
        x = inputs['x']
        return np.zeros_like(x, dtype=np.float64)
    
###############################################################################
class VDK05MOD(OMACMODSuper):
    """Realization of the null organic matter accretion model with the 
       van de Koppel et al. (2005) Bag scheme.

    Attributes:
        Parameters : rB0, dP, dB, Bmax, czh
    Constants:
        
    """ 
    
    # constructor
    def __init__(self, params):
        self.m_params = params
    
    def organic_accretion(self, inputs):
        """"Calculate organic matter accretion rate.
        Arguments:
            inputs : driving data for OM accretion calculation
        Returns: organic matter accretion rate (kg m-2 s-1)
        """
        x = inputs['x']
        return np.zeros_like(x, dtype=np.float64)
    
    def aboveground_biomass(self, inputs):
        """"Calculate aboveground biomass.
        Arguments:
            inputs : driving data for OM accretion calculation
        Returns: aboveground biomass (kg m-2)
        """
        rB0 = self.m_params['rB0']      # intrinsic growth rate (yr-1)
        Bmax = self.m_params['Bmax']    # maximal standing biomass (kg/m2)
        czh = self.m_params['czh']      # a half-saturation elev constant (m)
        dP = self.m_params['dP']        # plant mortality due to senescence (yr-1)
        dB = self.m_params['dB']        # plant mortality due to wave damage (yr-1)
        zh = inputs['zh']       # platform surface elevation (msl)
        S = inputs['S']         # platform surface slope (m/m)
        Bag_old = inputs['Bag'] # aboveground biomass of last time step (kg/m2)
        pft = inputs['pft']     # platform pft
        dt = inputs['dt']       # time step (s)
        A = rB0*(1-Bag_old/Bmax)*(zh/(zh+czh)) - dP - dB*S
        Bag = Bag_old * (1.0 + A*dt) / (1.0 - A*dt)
        indice = np.logical_or(pft<2, pft>5)
        Bag[indice] = 0.0
        return Bag

###############################################################################
class M12MOD(OMACMODSuper):
    """Realization of the Morris et al. (2012) organic matter accretion model.

    Attributes:

    Constants:
        
    """
    
    # constructor
    def __init__(self, params):
        self.m_params = params
        
    def organic_accretion(self, inputs):
        """"Calculate organic matter accretion rate.
        Arguments:
            inputs : driving data for OM accretion calculation
        Returns: 
            omac : organic matter accretion rate (kg m-2 s-1)
            Bag  : aboveground biomass (kg m-2)
        """

###############################################################################