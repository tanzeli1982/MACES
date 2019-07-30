#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 29 23:55:57 2019

Derived class for landward migration algorithms

@author: Zeli Tan
"""

import numpy as np
from TAIMODSuper import LNDMGMODSuper

class NULLLNDMG(LNDMGMODSuper):
    """Realization of the null landward migration model.

    Attributes:

    Constants:
        
    """ 
    
    # constructor
    def __init__(self, params):
        self.m_params = params
    
    def landward_migration(self, inputs):
        """"Calculate coastal wetland landward migration at the end of each year.
        Arguments:
            inputs : driving data for landward migration calculation
        Returns: the new pft on the coastal platform
        """
        pft = inputs['pft']
        return pft
    
class R20LNDMG(LNDMGMODSuper):
    """Realization of the Reyes et al. (2000) landward migration model.

    Attributes:

    Constants:
        
    """ 
    
    # constructor
    def __init__(self, params):
        self.m_params = params
    
    def landward_migration(self, inputs):
        """"Calculate coastal wetland landward migration at the end of each year.
        Arguments:
            inputs : driving data for landward migration calculation
        Returns: the new pft on the coastal platform
        """
        pft = inputs['pft']
        return pft 