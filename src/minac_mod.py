#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 29 23:06:01 2019

Derived class for mineral accretion algorithms

@author: Zeli Tan
"""

import numpy as np
import maces_utilities as utils
from TAIMODSuper import MACMODSuper

###############################################################################
class F06MOD(MACMODSuper):
    """Realization of the French (2006) mineral accretion model.

    Attributes:

    Constants:
        
    """
    
    # constructor
    def __init__(self, params):
        self.m_params = params
        
    def mineral_suspension(self, inputs):
        """"Calculate mineral suspension rate.
        Arguments:
            inputs : driving data for mineral suspension calculation
        Returns: mineral suspension rate (kg m-2 s-1)
        """
        x = inputs['x']
        return np.zeros_like(x, dtype=np.float64)
    
    def mineral_deposition(self, inputs):
        """"Calculate mineral deposition rate.
        Arguments:
            inputs : driving data for mineral deposition calculation
        Returns: mineral deposition rate (kg m-2 s-1)
        """
        Css = inputs['Css']     # sediment conc (kg/m3)
        tau = inputs['tau']     # bottom shear stress (Pa)
        d50 = inputs['d50']     # sediment median diameter (m)
        Rous = inputs['Rous']   # sediment density (kg/m3)
        ws = self.settling_velocity(tau, d50, Rous)
        return ws * Css

###############################################################################    
class T03MOD(MACMODSuper):
    """Realization of the Temmerman et al. (2003) mineral accretion model.

    Attributes:
        parameters k > 0, l < 0, m < 0
    Constants:
        
    """
    
    # constructor
    def __init__(self, params):
        self.m_params = params
        
    def mineral_suspension(self, inputs):
        """"Calculate mineral suspension rate.
        Arguments:
            inputs : driving data for mineral suspension calculation
        Returns: mineral suspension rate (kg m-2 s-1)
        """
        x = inputs['x']
        return np.zeros_like(x, dtype=np.float64)
    
    def mineral_deposition(self, inputs):
        """"Calculate mineral deposition rate.
        Arguments:
            inputs : driving data for mineral deposition calculation
        Returns: mineral deposition rate (kg m-2 s-1)
        """
        k = self.m_params['k']
        l = self.m_params['l']
        m = self.m_params['m']
        H = inputs['zh']
        Dc = inputs['x'] - inputs['xref']
        return 1e-3 * k * np.exp(l*H) * np.exp(m*Dc) / 27.32 / 8.64e4

###############################################################################    
class KM12MOD(MACMODSuper):
    """Realization of the Morris et al. (2012) mineral accretion model.

    Attributes:

    Constants:
        
    """
    
    # constructor
    def __init__(self, params):
        self.m_params = params
        
    def mineral_suspension(self, inputs):
        """"Calculate mineral suspension rate.
        Arguments:
            inputs : driving data for mineral suspension calculation
        Returns: mineral suspension rate (kg m-2 s-1)
        """
        x = inputs['x']
        rsuspend = np.zeros_like(x, dtype=np.float64)
        return rsuspend
    
    def mineral_deposition(self, inputs):
        """"Calculate mineral deposition rate.
        Arguments:
            inputs : driving data for mineral deposition calculation
        Returns: mineral deposition rate (kg m-2 s-1)
        """
        Css = inputs['Css']     # sediment conc (kg/m3)
        tau = inputs['tau']     # bottom shear stress (Pa)
        d50 = inputs['d50']     # sediment median diameter (m)
        Rous = inputs['Rous']   # sediment density (kg/m3)
        ws = self.settling_velocity(tau, d50, Rous)
        return ws * Css
    
    def settling_velocity(self, d50, Rous, Bag, U):
        """"Calculate effective sediment settling velocity (Morris et al., 2012).
        Arguments:
            tau : bottom shear stress (Pa)
            d50 : sediment median diameter (m)
            Rous : sediment density (kg/m3)
            Bag : aboveground biomass (kg/m2)
            U : tide flow velocity (m/s)
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
        # parameters
        ak = 0.9
        a0 = 11.0
        
        wup = 
        return ws - wup

###############################################################################    
class M12MOD(MACMODSuper):
    """Realization of the Morris et al. (2012) mineral accretion model.

    Attributes:

    Constants:
        
    """
    
    # constructor
    def __init__(self, params):
        self.m_params = params
        
    def mineral_suspension(self, inputs):
        """"Calculate mineral suspension rate.
        Arguments:
            inputs : driving data for mineral suspension calculation
        Returns: mineral suspension rate (kg m-2 s-1)
        """
        x = inputs['x']
        return np.zeros_like(x, dtype=np.float64)
    
    def mineral_deposition(self, inputs):
        """"Calculate mineral deposition rate.
        Arguments:
            inputs : driving data for mineral deposition calculation
        Returns: mineral deposition rate (kg m-2 s-1)
        """
        Css = inputs['Css']     # sediment conc (kg/m3)
        tau = inputs['tau']     # bottom shear stress (Pa)
        d50 = inputs['d50']     # sediment median diameter (m)
        Rous = inputs['Rous']   # sediment density (kg/m3)
        ws = self.settling_velocity(tau, d50, Rous)
        return ws * Css
    