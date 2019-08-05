#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 29 23:42:08 2019

Derived class for organic matter accretion algorithms

@author: Zeli Tan
"""

import numpy as np
from TAIMODSuper import OMACMODSuper

class NULLMOD(OMACMODSuper):
    """Realization of the null organic matter accretion model.

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
        Returns: organic matter accretion rate (g m-2 s-1) and aboveground
            biomass (gC m-2)
        """
        x = inputs['x']
        h = inputs['h']
        hMHT = inputs['DMHT']
        DMHT = hMHT - h
        aa = self.m_params['aa']
        bb = self.m_params['bb']
        cc = self.m_params['cc']
        Bag = aa * DMHT + bb * DMHT**2 + cc
        omac = np.zeros_like(x, dtype=np.float64)
        return omac, Bag

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
        Returns: organic matter accretion rate (g m-2 s-1)
        """
        