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
        parameters : k > 0, l < 0, m < 0
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
        parameters : alphaA, betaA, alphaD, betaD
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
        Roul = utils.Roul
        Karman = utils.Karman
        ws = (( np.sqrt(0.25*(A/F)**(2/m)+(4./3.*d50**3*G*(S-1)/F/nv**2)**(1/m)) \
              - 0.5*(A/F)**(1/m))**m) * nv / d50
        # parameters
        ak = 0.9
        a0 = 11.0
        chi = 0.46  # 0.46+/-0.11 (Tanino & Nepf, 2008)
        xi = 3.8    # 3.8+/-0.5 (Tanino & Nepf, 2008)
        alphaA = self.m_params['alphaA']
        betaA = self.m_params['betaA']
        alphaD = self.m_params['alphaD']
        betaD = self.m_params['betaD']
        aps = alphaA * Bag**betaA
        dps = alphaD * Bag**betaD
        cD = 2.0*(a0*nv/U/dps + chi + xi*0.25*np.pi*aps*dps)
        wup = Karman * np.sqrt(0.2*(ak**2)*(U**2)*(cD*aps*dps)**(2/3)/Roul)
        return ws - wup

###############################################################################    
class M12MOD(MACMODSuper):
    """Realization of the Morris et al. (2012) mineral accretion model.

    Attributes:
        parameters : E0, tauE_cr, ks, alphaSG, betaSG
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
        E0 = self.m_params['E0']    # the reference erosion rate (m/s)
        tauE_cr = self.m_params['tauE_cr']  # critical shear stress (Pa)
        Rous = inputs['Rous']       # sediment density (kg/m3)
        tau = inputs['tau']         # bottom shear stress (Pa)
        pft = inputs['pft']         # platform pft
        rsuspend = np.zeros_like(tau, dtype=np.float64)
        indice = np.logical_and( tau>tauE_cr, pft==1 )
        rsuspend[indice] = E0 * Rous * (tau[indice]/tauE_cr - 1.0)
        return rsuspend
    
    def mineral_deposition(self, inputs):
        """"Calculate mineral deposition rate.
        Arguments:
            inputs : driving data for mineral deposition calculation
        Returns: mineral deposition rate (kg m-2 s-1)
        """
        ks = self.m_params['ks']
        alphaSG = self.m_params['alphaSG']  # m2/yr
        betaSG = self.m_params['betaSG']    # m4/yr
        Css = inputs['Css']     # sediment conc (kg/m3)
        tau = inputs['tau']     # bottom shear stress (Pa)
        d50 = inputs['d50']     # sediment median diameter (m)
        Rous = inputs['Rous']   # sediment density (kg/m3)
        Bag = inputs['Bag']     # aboveground biomass (kg/m2)
        pft = inputs['pft']     # platform pft
        zh = inputs['zh']       # platform surface elevation (msl)
        MHT = inputs['MHT']     # mean high tide level (msl)
        TR = inputs['TR']       # tidal range (m)
        S = inputs['S']         # platform slope (m/m)
        ws = self.settling_velocity(tau, d50, Rous)
        DMHT = MHT - zh
        indice = np.logical_and(DMHT<=0, pft==1)
        DMHT[indice] = 0.0
        Qm = Rous * Css * (ws + ks*Bag) * (DMHT**2) / TR
        Qsg = Rous * (alphaSG/3.1536e7 - betaSG/3.1536e7*Bag) * S
        return Qm + Qsg
    
###############################################################################
class F07MOD(MACMODSuper):
    """Realization of the Fagherazzi et al. (2007) mineral accretion model.

    Attributes:
        parameters : E0, tauE_cr, tauD_cr, gamma
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
        E0 = self.m_params['E0']    # reference erosion rate (kg/m2/s)
        tauE_cr = self.m_params['tauE_cr']  # critical shear stress (Pa)
        gamma = self.m_params['gamma']  # the increase of critical shear stress with depth
        tau = inputs['tau']         # bottom shear stress (Pa)
        dtau = inputs['dtau']       # bottom shear stress gradient (Pa/s)
        dt = inputs['dt']           # time step (s)
        Eold = inputs['Eold']       # Esilt at the last time step (kg/m2/s)
        Esand = np.zeros_like(tau, dtype=np.float64)
        indice = tau > tauE_cr
        Esand[indice] = E0 * (tau[indice] - tauE_cr)**1.5
        Esilt = ((1-0.5*E0*gamma*dt)*Eold + E0*dtau*dt) / (1 + 0.5*E0*gamma*dt)
        return Esand + Esilt
        
    def mineral_deposition(self, inputs):
        """"Calculate mineral deposition rate.
        Arguments:
            inputs : driving data for mineral deposition calculation
        Returns: mineral deposition rate (kg m-2 s-1)
        """
        tauD_cr = self.m_params['tauD_cr']  # critical shear stress (Pa)
        Css = inputs['Css']     # sediment conc (kg/m3)
        tau = inputs['tau']     # bottom shear stress (Pa)
        Rous = inputs['Rous']   # sediment density (kg/m3) 
        ws = self.settling_velocity(tau, 1.25e-4, Rous)
        Qd_sand = 0.3 * Css * ws     # sand deposition (assume 30% of sediment)
        KD = 0.00077            # m5 s-1 kg-4/3
        Qd_silt = np.zeros_like(Css, dtype=np.float64)
        indice = tau < tauD_cr
        Qd_silt[indice] = KD * (0.7*Css[indice])**(7/3) * (1 - tau/tauD_cr)
        return Qd_sand + Qd_silt
    
###############################################################################
class VDK05MOD(MACMODSuper):
    """Realization of the van de Koppel et al. (2005) mineral accretion model.

    Attributes:
        parameters : Dmax, Emax, aNv, bNv
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
        Emax = self.m_params['Emax']    # maximum sediment erosion by water flow (yr-1)
        ds = self.m_params['ds']        # a conversion factor for wave erosion (yr-1)
        aNv = self.m_params['aNv']      # sediment erosion increase rate by plant (kg/m2)
        bNv = self.m_params['bNv']
        Bag = inputs['Bag']     # aboveground biomass (kg/m2)
        tau = inputs['tau']     # bottom shear stress (Pa)
        S = inputs['S']         # platform slope (m/m)
        tau_max = np.max(tau)
        
    def mineral_deposition(self, inputs):
        """"Calculate mineral deposition rate.
        Arguments:
            inputs : driving data for mineral deposition calculation
        Returns: mineral deposition rate (kg m-2 s-1)
        """
        Dmax = self.m_params['Dmax']    # maximum deposition rate (m/yr)
        zh = inputs['zh']           # platform surface elevation (msl)
        Rous = inputs['Rous']       # sediment density (kg/m3)
        Ks = inputs['MHT']          # mean high water level (msl)
        rdeposit = np.zeros_like(zh, dtype=np.float64)
        indice = np.logical_and(zh>=0, zh<=Ks)
        rdeposit[indice] = Rous * Dmax/3.1536e7 * (1.0 - zh[indice]/Ks)
        return rdeposit
        
        