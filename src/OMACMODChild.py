#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 29 23:42:08 2019

Derived class for organic matter accretion algorithms

@author: Zeli Tan
"""

import numpy as np
from TAIMODSuper import OMACMODSuper

class NULLMAC(OMACMODSuper):
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
        Returns: organic matter accretion rate (g m-2 s-1)
        """
        x = inputs['x']
        return np.zeros_like(x, dtype=np.float64)

class M12OMAC(OMACMODSuper):
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
        