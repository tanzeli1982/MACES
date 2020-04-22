#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 29 23:55:57 2019

Derived class for landward migration algorithms

@author: Zeli Tan
"""

import numpy as np
import maces_utilities as utils
from TAIMODSuper import LNDMGMODSuper

###############################################################################
class NULLMOD(LNDMGMODSuper):
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
        Bag = inputs['Bag']
        
        wtlnd_indice = np.logical_and(pft>=2,pft<=5) 
        if np.size(pft[wtlnd_indice])>0:
            ref_pft = pft[wtlnd_indice][0]
            indice = np.logical_and(pft==1, Bag>utils.TOL)
            pft[indice] = ref_pft
            indice = np.logical_and(wtlnd_indice, Bag<utils.TOL)
            pft[indice] = 1
        return pft

###############################################################################    
class R20MOD(LNDMGMODSuper):
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
    
###############################################################################
class S18MOD(LNDMGMODSuper):
    """Realization of the Schuerch et al. (2018) landward migration model.

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