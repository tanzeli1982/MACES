#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 29 23:51:00 2019

Derived class for storm surge erosion algorithms

@author: Zeli Tan
"""

import numpy as np
from TAIMODSuper import WINDEROMODSuper

class NULLWINDERO(WINDEROMODSuper):
    """Realization of the null storm surge erosion model.

    Attributes:

    Constants:
        
    """ 
    
    # constructor
    def __init__(self, params):
        self.m_params = params
        
    def wind_erosion(self, inputs):
        """"Calculate storm surge erosion rate. This should be used to get the
        new values of grid cell coordinate and elevation.
        Arguments:
            inputs : driving data for storm surge erosion calculation
        Returns: storm surge erosion rate (m s-1)
        """
        x = inputs['x']
        return np.zeros_like(x, dtype=np.float64)
    
class L16WINDERO(WINDEROMODSuper):
    """Realization of the Leonardi et al. (2016) storm surge erosion model.

    Attributes:

    Constants:
        
    """ 
    
    # constructor
    def __init__(self, params):
        self.m_params = params
        
    def wind_erosion(self, inputs):
        """"Calculate storm surge erosion rate. This should be used to get the
        new values of grid cell coordinate and elevation.
        Arguments:
            inputs : driving data for storm surge erosion calculation
        Returns: storm surge erosion rate (m s-1)
        """