#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 29 23:42:08 2019

Derived class for organic matter accretion algorithms

@author: Zeli Tan
"""

import numpy as np
import maces_utilities as utils
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
    
    def organic_deposition(self, inputs):
        """"Calculate organic matter deposition rate.
        Arguments:
            inputs : driving data for OM deposition calculation
        Returns: organic matter deposition rate (kg m-2 s-1)
        """
        DepOM = inputs['DepOM']     # OM deposition (kg/m2/s)
        DepOM[:] = 0.0
        return DepOM
    
    def aboveground_biomass(self, inputs):
        """"Calculate aboveground biomass (Morris et al., 2012).
        Arguments:
            inputs : driving data for OM accretion calculation
        Returns: aboveground biomass (kg m-2)
        """
        aa = self.m_params['aa']
        bb = self.m_params['bb']
        cc = self.m_params['cc']
        zh = inputs['zh']           # platform surface elevation (msl)
        #MHT = 0.5*inputs['TR']      # mean high tide water level (msl)
        MHT = inputs['MHHW']        # mean high high water level (msl)
        pft = inputs['pft']         # platform pft
        Bag = inputs['Bag']         # aboveground biomass (kg/m2)
        
        Bag[:] = 0.0
        indice = np.logical_and(np.logical_and(zh>=0,zh<=MHT), 
                                np.logical_and(pft>=1,pft<=5))
        DMHT = MHT - zh
        Bag[indice] = np.maximum(1e-3, aa[pft[indice]]*DMHT[indice]+ \
           bb[pft[indice]]*(DMHT[indice]**2)+cc[pft[indice]])
        return Bag
    
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
    
    def organic_deposition(self, inputs):
        """"Calculate organic matter deposition rate.
        Arguments:
            inputs : driving data for OM deposition calculation
        Returns: organic matter deposition rate (kg m-2 s-1)
        """
        DepOM = inputs['DepOM']     # OM deposition (kg/m2/s)
        DepOM[:] = 0.0
        return DepOM
    
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
        Bag = inputs['Bag']     # aboveground biomass (kg/m2)
        pft = inputs['pft']     # platform pft
        #MHT = 0.5*inputs['TR']      # mean high tide water level (msl)
        MHT = inputs['MHHW']        # mean high high water level (msl)
        
        Bag[:] = 0.0
        indice = np.logical_and(np.logical_and(pft>=2,pft<=5), 
                                np.logical_and(zh>=0,zh<=MHT))
        Bag[indice] = np.maximum(Bmax[pft[indice]] * (1.0 - (dP[pft[indice]]+ \
           dB[pft[indice]]*S[indice])*(1.0+czh[pft[indice]]/(0.1+zh[indice])) / \
           rB0[pft[indice]]), 1e-3)
        indice = np.logical_and(np.logical_and(zh>=0,zh<=MHT), pft==1)
        Bag[indice] = 1e-3
        return Bag
    
###############################################################################
class M12MOD(OMACMODSuper):
    """Realization of the Morris et al. (2012) organic matter accretion model.

    Attributes:
        parameters : Kr, Tr, phi, aa, bb, cc
    Constants:
        
    """
    
    # constructor
    def __init__(self, params):
        self.m_params = params
        
    def organic_deposition(self, inputs):
        """"Calculate organic matter deposition rate.
        Arguments:
            inputs : driving data for OM deposition calculation
        Returns: organic matter deposition rate (kg m-2 s-1)
        """
        Kr = self.m_params['Kr']    # the refractory fraction of root and rhizome biomass
        Tr = self.m_params['Tr']    # the root and rhizome turnover time (yr)
        phi = self.m_params['phi']  # the root:shoot quotient
        Bag = inputs['Bag']         # aboveground biomass (kg/m2)
        pft = inputs['pft']         # platform pft
        DepOM = inputs['DepOM']     # OM deposition (kg/m2/s)
        
        DepOM[:] = 0.0
        indice = np.logical_and(Bag>0, Tr[pft]>0)
        DepOM[indice] = Kr[pft[indice]]*(phi[pft[indice]]*Bag[indice])/ \
            (Tr[pft[indice]]*3.1536e7)
        return DepOM
    
    def aboveground_biomass(self, inputs):
        """"Calculate aboveground biomass (Morris et al., 2012).
        Arguments:
            inputs : driving data for OM accretion calculation
        Returns: aboveground biomass (kg m-2)
        """    
        aa = self.m_params['aa']
        bb = self.m_params['bb']
        cc = self.m_params['cc']
        zh = inputs['zh']           # platform surface elevation (msl)
        #MHT = 0.5*inputs['TR']      # mean high tide water level (msl)
        MHT = inputs['MHHW']        # mean high high water level (msl)
        pft = inputs['pft']         # platform pft
        Bag = inputs['Bag']         # aboveground biomass (kg/m2)
        
        Bag[:] = 0.0
        indice = np.logical_and(np.logical_and(zh>=0,zh<=MHT), 
                                np.logical_and(pft>=1,pft<=5))
        DMHT = MHT - zh
        Bag[indice] = np.maximum(1e-3, aa[pft[indice]]*DMHT[indice]+ \
           bb[pft[indice]]*(DMHT[indice]**2)+cc[pft[indice]])
        return Bag

###############################################################################
class DA07MOD(OMACMODSuper):
    """Realization of the D'Alpaos et al. (2007) organic matter accretion model.

    Attributes:
        parameters : Qom0, Bmax, omega, mps
    Constants:
        
    """
    
    # constructor
    def __init__(self, params):
        self.m_params = params
        
    def organic_deposition(self, inputs):
        """"Calculate organic matter deposition rate.
        Arguments:
            inputs : driving data for OM deposition calculation
        Returns: organic matter deposition rate (kg m-2 s-1)
        """
        Qom0 = self.m_params['Qom0']    # a typical OM deposition rate (m/yr)
        Bmax = self.m_params['Bmax']    # maximum Bag (kg/m2)
        rhoOM = self.m_params['rhoOM']  # OM density (kg/m3)
        Bag = inputs['Bag']             # aboveground biomass (kg/m2)
        pft = inputs['pft']             # platform pft
        DepOM = inputs['DepOM']         # OM deposition (kg/m2/s)
        
        DepOM[:] = 0.0
        indice = np.logical_and(Bag>utils.TOL, Bmax[pft]>0)
        DepOM[indice] = Qom0/3.1536e7 * rhoOM * Bag[indice] / Bmax[pft[indice]]
        return DepOM
        
    def aboveground_biomass(self, inputs):
        """"Calculate aboveground biomass.
        Arguments:
            inputs : driving data for OM accretion calculation
        Returns: aboveground biomass (kg m-2)
        """
        Bmax = self.m_params['Bmax']    # maximum Bag (kg/m2)
        omega = self.m_params['omega']  # the ratio of winter Bag to Bps 
        mps = self.m_params['mps']      # month of Bag at its peak
        Bag = inputs['Bag']         # aboveground biomass (kg/m2)
        zh = inputs['zh']           # platform surface elevation (msl)
        pft = inputs['pft']         # platform pft
        #MHT = 0.5*inputs['TR']      # mean high tide water level (msl)
        MHT = inputs['MHHW']        # mean high high water level (msl)
        m = inputs['month']         # month (1 to 12)
        
        Bag[:] = 0.0
        indice = np.logical_and(np.logical_and(pft>=2,pft<=5), 
                                np.logical_and(zh>=0,zh<=MHT))
        Bps = (MHT - zh[indice]) / MHT * Bmax[pft[indice]]   # peak season Bag
        Bag[indice] = np.maximum(1e-3, \
           0.5*Bps*(1-omega)*(np.sin(np.pi*m/6-mps*np.pi/12)+1) + omega*Bps)
        indice = np.logical_and(np.logical_and(zh>=0,zh<=MHT), pft==1)
        Bag[indice] = 1e-3
        return Bag

###############################################################################
class KM12MOD(OMACMODSuper):
    """Realization of the Kirwan & Mudd (2012) organic matter accretion model.

    Attributes:
        parameters : Bmax, sigmaB, rBmin, jdps, thetaBG, Dmbm, 
                     rGmin, rGps, sigmaOM, kl0, kr0
    Constants:
        
    """
    
    # constructor
    def __init__(self, params):
        self.m_params = params
        
    def organic_deposition(self, inputs):
        """"Calculate organic matter deposition rate.
        Arguments:
            inputs : driving data for OM deposition calculation
        Returns: organic matter deposition rate (kg m-2 s-1)
        """
        Kr = self.m_params['Kr']    # the refractory fraction of root and rhizome biomass
        Tr = self.m_params['Tr']    # the root and rhizome turnover time (yr)
        sigmaOM = self.m_params['sigmaOM']  # decay increase due to temperature (K-1)
        thetaBG = self.m_params['thetaBG']  # coef for the root:shoot quotient
        Dmbm = self.m_params['Dmbm']        # coef for the root:shoot quotient 
        Bmax = self.m_params['Bmax']        # maximum Bag (kg/m2)
        rBmin = self.m_params['rBmin']      # the ratio of winter Bag to Bps
        sigmaB = self.m_params['sigmaB']    # biomass increase due to temperature (K-1)
        rGmin = self.m_params['rGmin']      # the ratio of winter growth rate to Bps (day-1)
        rGps = self.m_params['rGps']        # the ratio of peak growth rate to Bps (day-1)
        jdps = self.m_params['jdps']        # the DOY when Bag is at its peak
        Tsummer = inputs['Tsummer']     # summer temperature (K)
        Tmean = inputs['Tmean']         # annual mean temperature (K)
        zh = inputs['zh']               # platform surface elevation (msl)
        pft = inputs['pft']             # platform pft
        MHHW = inputs['MHHW']           # mean high high water level (msl)
        jd = inputs['doy']              # day (1 to 365)
        DepOM = inputs['DepOM']         # OM deposition (kg/m2/s)
        
        DepOM[:] = 0.0
        jd_phi = 56     # the phase shift (in days) between Gps and Bps
        indice = np.logical_and(np.logical_and(zh>=0, zh<=MHHW), 
                                np.logical_and(pft>=2, pft<=5))
        # the root:shoot quotient
        phi = thetaBG[pft[indice]]*(MHHW-zh[indice]) + Dmbm[pft[indice]]
        # peak season Bag
        Tref = 293.15   # K
        Bps = Bmax[pft[indice]]*(MHHW-zh[indice])/MHHW* \
            (1+(Tsummer-Tref)*sigmaB[pft[indice]])
        Bmin = rBmin * Bps          # winter Bag
        Gmin = rGmin/8.64e4 * Bps   # winter growth rate (kg/m2/s)
        Gps = rGps/8.64e4 * Bps     # peak growth rate (kg/m2/s)
        # the mortality rate (kg/m2/s) of aboveground biomass
        TrefOM = 273.15     # K
        Tr_adjust = Tr[pft[indice]] / max(1.0+(Tmean-TrefOM)*sigmaOM, 0.1)
        Mag = ( 0.5*(Gmin+Gps+(Gps-Gmin)*np.cos(2.0*np.pi*(jd-jdps+jd_phi)/365)) + \
            np.pi*(Bps-Bmin)/3.1536e7*np.sin(2.0*np.pi*(jd-jdps)/365) ) * \
            Kr[pft[indice]] / Tr_adjust
        DepOM[indice] = np.maximum(phi,0.0) * np.maximum(Mag,0.0)
        return DepOM
        
    def aboveground_biomass(self, inputs):
        """"Calculate aboveground biomass.
        Arguments:
            inputs : driving data for OM accretion calculation
        Returns: aboveground biomass (kg m-2)
        """
        Bmax = self.m_params['Bmax']        # maximum Bag (kg/m2)
        rBmin = self.m_params['rBmin']      # the ratio of winter Bag to Bps
        sigmaB = self.m_params['sigmaB']    # biomass increase due to temperature (K-1)
        jdps = self.m_params['jdps']        # the DOY when Bag is at its peak
        Tsummer = inputs['Tsummer']     # summer temperature (K)
        zh = inputs['zh']               # platform surface elevation (msl)
        pft = inputs['pft']             # platform pft
        Bag = inputs['Bag']             # aboveground biomass (kg/m2)
        MHHW = inputs['MHHW']           # mean high high water level (msl)
        jd = inputs['doy']              # day (1 to 365)
        
        Bag[:] = 0.0
        Tref = 293.15   # K
        indice = np.logical_and(np.logical_and(pft>=2,pft<=5), 
                                np.logical_and(zh>=0,zh<=MHHW))
        Bps = Bmax[pft[indice]]*(MHHW-zh[indice])/MHHW* \
            (1+(Tsummer-Tref)*sigmaB[pft[indice]])
        Bmin = rBmin * Bps          # winter Bag
        Bag[indice] = np.maximum(0.5*(Bmin+Bps+(Bps-Bmin)* \
           np.cos(2*np.pi*(jd-jdps)/365)), 1e-3)
        indice = np.logical_and(np.logical_and(zh>=0,zh<=MHHW), pft==1)
        Bag[indice] = 1e-3
        return Bag
    
    def belowground_biomass(self, inputs):
        """"Calculate belowground biomass.
        Arguments:
            inputs : driving data for belowground biomass calculation
        Returns: belowground biomass (kg m-2)
        """
        thetaBG = self.m_params['thetaBG']  # coef for the root:shoot quotient
        Dmbm = self.m_params['Dmbm']        # coef for the root:shoot quotient    
        Bag = inputs['Bag']         # aboveground biomass (kg/m2)
        Bbg = inputs['Bbg']         # belowground biomass (kg/m2)
        zh = inputs['zh']           # platform surface elevation (msl)
        pft = inputs['pft']         # platform pft
        MHHW = inputs['MHHW']       # mean high high water level (msl)
        
        indice = np.logical_and(np.logical_and(zh>=0, zh<=MHHW), 
                                np.logical_and(pft>=2, pft<=5))
        phi = thetaBG[pft[indice]]*(MHHW-zh[indice]) + Dmbm[pft[indice]]
        Bbg[indice] = np.maximum(phi,0.0) * np.maximum(Bag[indice],0.0)
        return Bbg
        
###############################################################################
class K16MOD(OMACMODSuper):
    """Realization of the Kakeh et al. (2016) organic matter accretion model.

    Attributes:
        parameters : gammaB, Bmax, Gmgv, b2mgv, b3mgv, Mdmax, Mhmax, phi
    Constants:
        
    """
    
    Md = None       # individual mangrove tree diameter (cm)
    Mh = None       # individual mangrove tree height (cm)
    
    # constructor
    def __init__(self, params):
        self.m_params = params
        # solve Md and Mh
        b2mgv = self.m_params['b2mgv']  # coef for Md vs Mh equation (dimensionless)
        b3mgv = self.m_params['b3mgv']  # coef for Md vs Mh equation (cm-1)
        Mhmax = self.m_params['Mhmax']  # mangrove maximum height (cm)
        Mdmax = self.m_params['Mdmax']  # mangrove maximum diameter (cm)
        coefs = [b3mgv, b2mgv, 137, -Mhmax*Mdmax]
        roots = np.roots(coefs)
        self.Md = roots[roots>0][0]
        self.Mh = 137 + b2mgv*self.Md + b3mgv*(self.Md)**2
        
    def organic_deposition(self, inputs):
        """"Calculate organic matter deposition rate.
        gammaB within [1e-3,3e-3]
        Arguments:
            inputs : driving data for OM deposition calculation
        Returns: organic matter deposition rate (kg m-2 s-1)
        """
        gammaB = self.m_params['gammaB']    # m yr-1 m2 kg-1
        rhoOM = self.m_params['rhoOM']      # OM density (kg/m3)
        DepOM = inputs['DepOM']         # OM deposition (kg/m2/s)
        Bag = inputs['Bag']         # aboveground biomass (kg/m2)
        
        DepOM[:] = rhoOM * gammaB/3.1536e7 * Bag
        return DepOM
        
    def aboveground_biomass(self, inputs):
        """"Calculate aboveground biomass.
        Arguments:
            inputs : driving data for OM accretion calculation
        Returns: aboveground biomass (kg m-2)
        """
        Bmax = self.m_params['Bmax']    # maximum Bag (kg/m2)
        b2mgv = self.m_params['b2mgv']  # coef for Md vs Mh equation (dimensionless)
        b3mgv = self.m_params['b3mgv']  # coef for Md vs Mh equation (cm-1)
        Bag = inputs['Bag']         # aboveground biomass (kg/m2)
        zh = inputs['zh']           # platform surface elevation (msl)
        pft = inputs['pft']         # platform pft
        #MHT = 0.5*inputs['TR']      # mean high tide water level (msl)
        MHT = inputs['MHHW']        # mean high high water level (msl)
        
        Bag[:] = 0.0
        indice_zh = np.logical_and(zh>=0, zh<=MHT)
        # tidal flats
        indice = np.logical_and(indice_zh, pft==1)
        Bag[indice] = 1e-3
        # Spartina alterniflora dominated marshes
        indice = np.logical_and(indice_zh, pft==2)
        rz = (1-0.5*zh[indice]/MHT)   # Bag production rate
        mz = 0.5*zh[indice]/MHT       # Bag mortality rate
        Bag[indice] = np.maximum( Bmax[pft[indice]]*(1.0-mz/rz), 1e-3 )
        # multi-species marshes
        indice = np.logical_and(np.logical_or(pft==3,pft==4), indice_zh)
        rz = 0.5*(1+zh[indice]/MHT)
        mz = 0.5*(1-zh[indice]/MHT)
        Bag[indice] = np.maximum( Bmax[pft[indice]]*(1.0-mz/rz), 1e-3 )
        # mangroves
        indice = np.logical_and(indice_zh, pft==5)
        Md = np.zeros_like(zh[indice])
        Mh = np.zeros_like(zh[indice])
        P = 1 - zh[indice]/MHT
        I = 4*P - 8*P**2 + 0.5
        Md[I<=0] = self.Md
        Mh[I<=0] = self.Mh
        Md[I>0] = self.Md * (1-0.5*zh[indice][I>0]/MHT)
        Mh[I>0] = 137 + b2mgv*Md[I>0] + b3mgv*(Md[I>0])**2
        #rout = 0.5/Md   # tree density (tree/m2)
        Bag[indice] = 0.154 * Md**1.11    # kg/m2
        return Bag
        
    def belowground_biomass(self, inputs):
        """"Calculate belowground biomass.
        Arguments:
            inputs : driving data for belowground biomass calculation
        Returns: belowground biomass (kg m-2)
        """
        phi = self.m_params['phi']  # the root:shoot quotient
        Bbg = inputs['Bbg']         # belowground biomass (kg/m2)
        Bag = inputs['Bag']         # aboveground biomass (kg/m2)
        pft = inputs['pft']         # platform pft
        zh = inputs['zh']           # platform surface elevation (msl)
        #MHT = 0.5*inputs['TR']      # mean high tide water level (msl)
        MHT = inputs['MHHW']        # mean high high water level (msl)
        
        indice_zh = np.logical_and(zh>=0, zh<=MHT)
        # marshes
        indice = np.logical_and(np.logical_and(pft>=2, pft<=4), 
                                indice_zh)
        Bbg[indice] = phi[pft[indice]]*Bag[indice]
        # mangroves
        indice = np.logical_and(indice_zh, pft==5)
        Bbg[indice] = 0.64 * (Bag[indice]/0.154)**(0.17/1.11)
        return Bbg
         
