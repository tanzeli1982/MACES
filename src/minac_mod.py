#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 29 23:06:01 2019

Derived class for mineral accretion algorithms

@author: Zeli Tan
"""

import numpy as np
from TAIMODSuper import MACMODSuper

class F06MAC(MACMODSuper):
    """Realization of the French (2006) mineral accretion model.

    Attributes:

    Constants:
        
    """
    
    # constructor
    def __init__(self, params):
        self.m_params = params
        
    def mineral_suspend(self, inputs):
        """"Calculate mineral suspension rate.
        Arguments:
            inputs : driving data for mineral suspension calculation
        Returns: mineral suspension rate (g m-2 s-1)
        """
        Css = inputs['Css']     # sediment conc (g/m3)
        tau = inputs['tau']     # bottom shear stress (Pa)
        d50 = inputs['d50']     # sediment median diameter (m)
        rsuspend = ws * Css
        return rsuspend
    
    def mineral_deposition(self, inputs):
        """"Calculate mineral deposition rate.
        Arguments:
            inputs : driving data for mineral deposition calculation
        Returns: mineral deposition rate (g m-2 s-1)
        """
        x = inputs['x']
        rdeposit = np.zeros_like(x, dtype=np.float64)
        return rdeposit