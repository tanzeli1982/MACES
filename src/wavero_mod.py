#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 29 23:51:00 2019

Derived class for storm surge erosion algorithms

@author: Zeli Tan
"""

import numpy as np
from TAIMODSuper import WAVEROMODSuper

###############################################################################
class NULLMOD(WAVEROMODSuper):
    """Realization of the null storm surge erosion model.

    Attributes:

    Constants:
        
    """ 
    
    # constructor
    def __init__(self, params):
        self.m_params = params
        
    def wave_erosion(self, inputs):
        """"Calculate storm surge erosion rate. This should be used to get the
        new values of grid cell coordinate and elevation.
        Arguments:
            inputs : driving data for storm surge erosion calculation
        Returns: storm surge erosion rate (m s-1)
        """
        x = inputs['x']
        return np.zeros_like(x, dtype=np.float64)

###############################################################################    
class L16MOD(WAVEROMODSuper):
    """Realization of the Leonardi et al. (2016) storm surge erosion model.

    Attributes:

    Constants:
        
    """ 
    
    # constructor
    def __init__(self, params):
        self.m_params = params
        
    def wave_erosion(self, inputs):
        """"Calculate storm surge erosion rate. This should be used to get the
        new values of grid cell coordinate and elevation.
        Arguments:
            inputs : driving data for storm surge erosion calculation
        Returns: storm surge erosion rate (m s-1)
        """
        x = inputs['x']
        return np.zeros_like(x, dtype=np.float64)
    
    def calculate_wave_erosion(self, inputs); 

        dWave_power = self.calculate_wave_power()
        #dRate_erosion units m / s
        #dKe is the erodibility coefficient
        dRate_erosion = dKe * dWave_power / dDistance


        print('Finished')

    #calculate wave power
    def calculate_wave_power(self):

        dCelerity_wave = self.m_params['dCelerity_wave']

        dWave_power = 0.125 * dDensity_water * g * dHeight_wave * dHeight_wave * dCelerity_wave

        return dWave_power