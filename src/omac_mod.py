#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 29 23:42:08 2019

Derived class for organic matter accretion algorithms

@author: Zeli Tan
"""

import numpy as np
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
        x = inputs['x']
        return np.zeros_like(x, dtype=np.float64, order='F')
    
    def aboveground_biomass(self, inputs):
        """"Calculate aboveground biomass (Morris et al., 2012).
        Arguments:
            inputs : driving data for OM accretion calculation
        Returns: aboveground biomass (kg m-2)
        """
        aa = self.m_params['aa']
        bb = self.m_params['bb']
        cc = self.m_params['cc']
        zh = inputs['zh']       # platform surface elevation (msl)
        MHT = inputs['MHT']     # mean high tide water level (msl)
        pft = inputs['pft']     # platform pft
        DMHT = MHT - zh
        Bag = aa * DMHT + bb * DMHT**2 + cc
        indice = np.logical_and(np.logical_and(zh>=0, zh<=MHT), 
                                np.logical_and(pft>=2, pft<=5))
        Bag[np.logical_not(indice)] = 0.0
        Bag[Bag<0] = 0.0
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
        x = inputs['x']
        return np.zeros_like(x, dtype=np.float64, order='F')
    
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
        Bag_old = inputs['Bag'] # aboveground biomass of last time step (kg/m2)
        pft = inputs['pft']     # platform pft
        dt = inputs['dt']       # time step (s)
        A = rB0*(1-Bag_old/Bmax)*(zh/(zh+czh)) - dP - dB*S
        Bag = Bag_old * (1.0 + A*dt) / (1.0 - A*dt)
        indice = np.logical_or(Bag<0, np.logical_or(pft<2, pft>5))
        Bag[indice] = 0.0
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
        Tr = self.m_params['Tr']    # the root aand rhizome turnover time (yr)
        phi = self.m_params['phi']  # the root:shoot quotient
        Bag = inputs['Bag']         # aboveground biomass (kg/m2)
        return Kr * (phi*Bag) / (Tr*3.1536e7)
    
    def aboveground_biomass(self, inputs):
        """"Calculate aboveground biomass (Morris et al., 2012).
        Arguments:
            inputs : driving data for OM accretion calculation
        Returns: aboveground biomass (kg m-2)
        """
        aa = self.m_params['aa']
        bb = self.m_params['bb']
        cc = self.m_params['cc']
        zh = inputs['zh']       # platform surface elevation (msl)
        MHT = inputs['MHT']     # mean high tide water level (msl)
        pft = inputs['pft']     # platform pft
        DMHT = MHT - zh
        Bag = aa * DMHT + bb * DMHT**2 + cc
        indice = np.logical_and(np.logical_and(zh>=0, zh<=MHT), 
                                np.logical_and(pft>=2, pft<=5))
        Bag[np.logical_not(indice)] = 0.0
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
        Qom0 = self.m_params['Qom0']    # a typical OM deposition rate (kg/m2/s)
        Bmax = self.m_params['Bmax']    # maximum Bag (kg/m2)
        Bag = inputs['Bag']             # aboveground biomass (kg/m2)
        return Qom0 * Bag / Bmax
        
    def aboveground_biomass(self, inputs):
        """"Calculate aboveground biomass.
        Arguments:
            inputs : driving data for OM accretion calculation
        Returns: aboveground biomass (kg m-2)
        """
        Bmax = self.m_params['Bmax']    # maximum Bag (kg/m2)
        omega = self.m_params['omega']  # the ratio of winter Bag to Bps 
        mps = self.m_params['mps']      # month of Bag at its peak
        zh = inputs['zh']       # platform surface elevation (msl)
        MHT = inputs['MHT']     # mean high tide water level (msl)
        m = inputs['MONTH']     # month (1 to 12)
        Bps = (MHT - zh) / MHT * Bmax   # peak season Bag
        Bps[np.logical_or(zh>MHT,zh<0)] = 0.0
        return 0.5*Bps*(1-omega)*(np.sin(np.pi*m/6-mps*np.pi/12)+1) + omega*Bps

###############################################################################
class KM12MOD(OMACMODSuper):
    """Realization of the Kirwan & Mudd (2012) organic matter accretion model.

    Attributes:
        parameters : Bmax, Tref, sigmaB, rBmin, jdps, thetaBG, Dmbm, 
                     rGmin, rGps, sigmaOM, TrefOM, kl0, kr0
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
        thetaBG = self.m_params['thetaBG']  # coef for the root:shoot quotient
        Dmbm = self.m_params['Dmbm']        # coef for the root:shoot quotient 
        Bmax = self.m_params['Bmax']        # maximum Bag (kg/m2)
        rBmin = self.m_params['rBmin']      # the ratio of winter Bag to Bps
        Tref = self.m_params['Tref']        # reference temperature for veg growth (K)
        sigmaB = self.m_params['sigmaB']    # biomass increase due to temperature (K-1)
        rGmin = self.m_params['rGmin']      # the ratio of winter growth rate to Bps (day-1)
        rGps = self.m_params['rGps']        # the ratio of peak growth rate to Bps (day-1)
        jdps = self.m_params['jdps']        # the DOY when Bag is at its peak
        Tair = inputs['Tair']       # air temperature (K)
        zh = inputs['zh']           # platform surface elevation (msl)
        MHHW = inputs['MHHW']       # mean high high water level (msl)
        jd = inputs['day']          # day (1 to 365)
        phi = thetaBG*(MHHW-zh) + Dmbm      # the root:shoot quotient
        phi[np.logical_or(zh<0,zh>MHHW)] = 0.0
        Bps = Bmax * (MHHW-zh) / MHHW * (1 + (Tair-Tref)*sigmaB)   # peak season Bag
        Bmin = rBmin * Bps          # winter Bag
        Gmin = rGmin/8.64e4 * Bps   # winter growth rate (kg/m2/s)
        Gps = rGps/8.64e4 * Bps     # peak growth rate (kg/m2/s)
        jd_phi = 56     # the phase shift (in days) between Gps and Bps
        # the mortality rate (kg/m2/s) of aboveground biomass
        Mag = 0.5*(Gmin+Gps+(Gps-Gmin)*np.cos(2.0*np.pi*(jd-jdps+jd_phi)/365)) + \
            np.pi/365*(Bps-Bmin)*np.sin(2.0*np.pi*(jd-jdps)/365)
        return phi*Mag
        
    def aboveground_biomass(self, inputs):
        """"Calculate aboveground biomass.
        Arguments:
            inputs : driving data for OM accretion calculation
        Returns: aboveground biomass (kg m-2)
        """
        Bmax = self.m_params['Bmax']        # maximum Bag (kg/m2)
        rBmin = self.m_params['rBmin']      # the ratio of winter Bag to Bps
        Tref = self.m_params['Tref']        # reference temperature for veg growth (K)
        sigmaB = self.m_params['sigmaB']    # biomass increase due to temperature (K-1)
        jdps = self.m_params['jdps']        # the DOY when Bag is at its peak
        Tair = inputs['Tair']       # air temperature (K)
        zh = inputs['zh']           # platform surface elevation (msl)
        MHHW = inputs['MHHW']       # mean high high water level (msl)
        jd = inputs['day']          # day (1 to 365)
        Bps = Bmax * (MHHW-zh) / MHHW * (1 + (Tair-Tref)*sigmaB)   # peak season Bag
        Bps[np.logical_or(zh<0,zh>MHHW)] = 0.0
        Bmin = rBmin * Bps          # winter Bag
        return 0.5*(Bmin+Bps+(Bps-Bmin)*np.cos(2*np.pi*(jd-jdps)/365))
    
    def belowground_biomass(self, inputs):
        """"Calculate belowground biomass.
        Arguments:
            inputs : driving data for belowground biomass calculation
        Returns: belowground biomass (kg m-2)
        """
        thetaBG = self.m_params['thetaBG']  # coef for the root:shoot quotient
        Dmbm = self.m_params['Dmbm']        # coef for the root:shoot quotient    
        Bag = inputs['Bag']         # aboveground biomass (kg/m2)
        zh = inputs['zh']           # platform surface elevation (msl)
        MHHW = inputs['MHHW']       # mean high high water level (msl)
        phi = thetaBG*(MHHW-zh) + Dmbm      # the root:shoot quotient
        phi[np.logical_or(zh<0,zh>MHHW)] = 0.0
        return phi*Bag
    
    def soilcarbon_decay(self, inputs):
        """"Calculate soil OC mineralization rate.
        Arguments:
            inputs : driving data for SOC decay rate calculation
        Returns: SOC decay rate (kg m-2 s-1) of two pools
        """
        kl0 = self.m_params['kl0']  # column-integrated decay rate of labile pool (yr-1)
        kr0 = self.m_params['kr0']  # column-integrated decay rate of refractory pool (yr-1)
        TrefOM = self.m_params['TrefOM']    # reference temperature for decay (K)
        sigmaOM = self.m_params['sigmaOM']  # decay increase due to temperature (K-1)
        SOM = inputs['SOM']     # soil organic matter pools (kg/m2)
        Tsoi = inputs['Tsoi']   # soil temperature (K)
        Cl = SOM[:,0]           # labile belowground SOM pool
        Cr = SOM[:,1]           # refractory belowground SOM pool
        Nx = np.shape(SOM)[0]
        npool = np.shape(SOM)[1]
        rdC = np.zeros((Nx,npool), dtype=np.float64, order='F')
        rdC[:,0] = ((1.0+(Tsoi-TrefOM)*sigmaOM)*kl0/3.1536e7) * Cl
        rdC[:,1] = ((1.0+(Tsoi-TrefOM)*sigmaOM)*kr0/3.1536e7) * Cr
        return rdC
        
###############################################################################
class K16MOD(OMACMODSuper):
    """Realization of the Kakeh et al. (2016) organic matter accretion model.

    Attributes:
        parameters : gammaB, Bmax, Gmgv, b2mgv, b3mgv, Mdmax, Mhmax, phi
    Constants:
        
    """
    
    tmp_Md = None       # individual mangrove tree diameter (cm)
    tmp_Mh = None       # individual mangrove tree height (cm)
    
    # constructor
    def __init__(self, params, Md0, Mh0):
        self.m_params = params
        self.tmp_Md = Md0
        self.tmp_Mh = Mh0
        
    def organic_deposition(self, inputs):
        """"Calculate organic matter deposition rate.
        Arguments:
            inputs : driving data for OM deposition calculation
        Returns: organic matter deposition rate (kg m-2 s-1)
        """
        gammaB = self.m_params['gammaB']    # m yr-1 m2 kg-1
        Bag = inputs['Bag']         # aboveground biomass (kg/m2)
        rhoOM = inputs['rhoOM']     # OM density (kg/m3)
        return rhoOM * gammaB/3.1536e7 * Bag
        
    def aboveground_biomass(self, inputs):
        """"Calculate aboveground biomass.
        Arguments:
            inputs : driving data for OM accretion calculation
        Returns: aboveground biomass (kg m-2)
        """
        Bmax = self.m_params['Bmax']    # maximum Bag (kg/m2)
        Gmgv = self.m_params['Gmgv']    # stem diameter growth rate (cm s-1)
        b2mgv = self.m_params['b2mgv']  # coef for Md vs Mh equation (dimensionless)
        b3mgv = self.m_params['b3mgv']  # coef for Md vs Mh equation (cm-1)
        Mhmax = self.m_params['Mhmax']  # mangrove maximum height (cm)
        Mdmax = self.m_params['Mdmax']  # mangrove maximum diameter (cm)
        Bag_old = inputs['Bag']     # Bag at the last time step (kg/m2)
        zh = inputs['zh']           # platform surface elevation (msl)
        pft = inputs['pft']         # platform pft
        MHT = inputs['MHT']         # mean high tide water level (msl)
        dt = inputs['dt']           # time step (s)
        Nx = np.size(zh)
        Bag = np.zeros(Nx, dtype=np.float64, order='F')
        for ii in range(Nx):
            if pft[ii]==2:
                # Spartina alterniflora dominated marshes
                if zh[ii]>=0 and zh[ii]<=MHT:
                    rz = (1-0.5*zh[ii]/MHT)/3.1536e7   # Bag production rate (s-1)
                    mz = 0.5*zh[ii]/MHT/3.1536e7       # Bag mortality rate (s-1)
                    Bag[ii] = max( (1+0.5*(rz*(1-Bag_old[ii]/Bmax)-mz)*dt)* \
                       Bag_old[ii]/(1-0.5*(rz*(1-Bag_old[ii]/Bmax)-mz)*dt), 0.0 )
                else:
                    Bag[ii] = Bag_old[ii]
            elif pft[ii]==3 or pft[ii]==4:
                # multi-species marshes
                if zh[ii]>=0 and zh[ii]<=MHT:
                    rz = 0.5*(1+zh[ii]/MHT)
                    mz = 0.5*(1-zh[ii]/MHT)
                    Bag[ii] = max( (1+0.5*(rz*(1-Bag_old[ii]/Bmax)-mz)*dt)* \
                       Bag_old[ii]/(1-0.5*(rz*(1-Bag_old[ii]/Bmax)-mz)*dt), 0.0 )
                else:
                    Bag[ii] = Bag_old[ii]
            elif pft[ii]==5:
                # mangroves (Avicennia marina)
                if zh[ii]>=0 and zh[ii]<=MHT:
                    P = 1 - zh[ii]/MHT
                    I = max(4*P-8*P**2+0.5, 0.0)
                    self.tmp_Md[ii] = self.tmp_Md[ii] + I*Gmgv/Mdmax/Mhmax* \
                        self.tmp_Md[ii]*(1-self.tmp_Md[ii]*self.tmp_Mh[ii])/ \
                        (274+3*b2mgv*self.tmp_Md[ii]-4*b3mgv*self.tmp_Md[ii]**2)
                    self.tmp_Mh[ii] = 137 + b2mgv*self.tmp_Md[ii] + \
                        b3mgv*self.tmp_Md[ii]**2
                    rout = 0.5/self.tmp_Md[ii]   # tree density (tree/m2)
                    Bag[ii] = rout * 0.308*self.tmp_Md[ii]**2.11    # kg/m2
                else:
                    Bag[ii] = Bag_old[ii]
        return Bag
        
    def belowground_biomass(self, inputs):
        """"Calculate belowground biomass.
        Arguments:
            inputs : driving data for belowground biomass calculation
        Returns: belowground biomass (kg m-2)
        """
        phi = self.m_params['phi']  # the root:shoot quotient
        Bbg_old = inputs['Bbg']     # belowground biomass of the last time step (kg/m2)
        Bag = inputs['Bag']         # aboveground biomass (kg/m2)
        pft = inputs['pft']         # platform pft
        zh = inputs['zh']           # platform surface elevation (msl)
        MHT = inputs['MHT']         # mean high tide water level (msl)
        Nx = np.size(zh)
        Bbg = np.zeros(Nx, dtype=np.float64, order='F')
        for ii in range(Nx):
            if pft[ii]>=2 and pft[ii]<=4:
                if zh[ii]>=0 and zh[ii]<=MHT:
                    Bbg[ii] = phi*Bag[ii]
                else:
                    Bbg[ii] = Bbg_old[ii]
            elif pft[ii]==5:
                # mangroves (Avicennia marina)
                if zh[ii]>=0 and zh[ii]<=MHT:
                    rout = 0.5/self.tmp_Md[ii]  # tree density (tree/m2)
                    Bbg[ii] = rout * 1.28*self.tmp_Md[ii]**1.17
                else:
                    Bbg[ii] = Bbg_old[ii]
        return Bbg
         