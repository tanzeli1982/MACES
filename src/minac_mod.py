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
        parameters : d50, rhoSed, porSed
    Constants:
        
    """
    
    # constructor
    def __init__(self, params):
        self.m_params = params
        self.m_update_Css = False
        
    def mineral_suspension(self, inputs):
        """"Calculate mineral suspension rate.
        Arguments:
            inputs : driving data for mineral suspension calculation
        Returns: mineral suspension rate (kg m-2 s-1)
        """
        Esed = inputs['Esed']   # sediment erosion (kg/m2/s)
        Esed[:] = 0.0
        return Esed
    
    def mineral_deposition(self, inputs):
        """"Calculate mineral deposition rate.
        Arguments:
            inputs : driving data for mineral deposition calculation
        Returns: mineral deposition rate (kg m-2 s-1)
        """
        pft = inputs['pft']     # platform pft
        Css0 = inputs['refCss'] # reference sediment conc (kg/m3)
        tau = inputs['tau']     # bottom shear stress (Pa)
        Dsed = inputs['Dsed']   # sediment deposition (kg/m2/s)
        
        ws = self.settling_velocity(tau)
        Dsed[:] = ws * Css0
        Dsed[pft==1] = 0.0
        return Dsed

###############################################################################    
class T03MOD(MACMODSuper):
    """Realization of the Temmerman et al. (2003) mineral accretion model.

    Attributes:
        parameters : k > 0, l < 0, m < 0, d50, rhoSed, porSed
    Constants:
        
    """
    
    # constructor
    def __init__(self, params):
        self.m_params = params
        self.m_update_Css = False
        
    def mineral_suspension(self, inputs):
        """"Calculate mineral suspension rate.
        Arguments:
            inputs : driving data for mineral suspension calculation
        Returns: mineral suspension rate (kg m-2 s-1)
        """
        Esed = inputs['Esed']   # sediment erosion (kg/m2/s)
        Esed[:] = 0.0
        return Esed
    
    def mineral_deposition(self, inputs):
        """"Calculate mineral deposition rate.
        Arguments:
            inputs : driving data for mineral deposition calculation
        Returns: mineral deposition rate (kg m-2 s-1)
        """
        k = self.m_params['k']  # reference deposition (kg/m2/spring cycle)
        l = self.m_params['l']  # elevation factor (m-1)
        m = self.m_params['m']  # distance factor (m-1)
        pft = inputs['pft']     # platform pft
        H = inputs['zh']        # platform surface elevation (msl)
        Dsed = inputs['Dsed']   # sediment deposition (kg/m2/s)
        
        Dsed[:] = 0.0
        Dc = inputs['x'] - inputs['xref']   # distance to wetland edge (m)
        indice = np.logical_and(np.logical_and(pft>1,pft<=9), 
                                np.logical_and(H>=0,Dc>=0))
        Dsed[indice] = k * np.exp(l*H[indice]) * np.exp(m*Dc[indice]) / 27.32 / 8.64e4
        return Dsed

###############################################################################    
class KM12MOD(MACMODSuper):
    """Realization of the Morris et al. (2012) mineral accretion model.

    Attributes:
        parameters : d50, rhoSed, porSed, alphaA, betaA, alphaD, betaD
    Constants:
        
    """
    
    # constructor
    def __init__(self, params):
        self.m_params = params
        self.m_update_Css = False
        
    def mineral_suspension(self, inputs):
        """"Calculate mineral suspension rate.
        Arguments:
            inputs : driving data for mineral suspension calculation
        Returns: mineral suspension rate (kg m-2 s-1)
        """
        Esed = inputs['Esed']   # sediment erosion (kg/m2/s)
        Esed[:] = 0.0
        return Esed
    
    def mineral_deposition(self, inputs):
        """"Calculate mineral deposition rate.
        Arguments:
            inputs : driving data for mineral deposition calculation
        Returns: mineral deposition rate (kg m-2 s-1)
        """
        pft = inputs['pft']     # platform pft
        Css0 = inputs['refCss'] # reference sediment conc (kg/m3)
        Bag = inputs['Bag']     # aboveground biomass (kg/m2)
        U = inputs['U']         # flow velocity (m/s)
        Dsed = inputs['Dsed']   # sediment deposition (kg/m2/s)
        
        Dsed[:] = 0.0
        ws = self.settling_velocity(pft, 1e3*Bag, U)
        indice = np.logical_and(pft>1,pft<=9)
        Dsed[indice] = ws[indice] * Css0
        return Dsed
    
    def settling_velocity(self, pft, Bag, U):
        """"Calculate effective sediment settling velocity (Morris et al., 2012).
        Arguments:
            d50 : sediment median diameter (m)
            Rous : sediment density (kg/m3)
            pft : platform pft
            Bag : aboveground biomass (gC/m2)
            U : tide flow velocity (m/s)
        Returns: sediment settling velocity (m s-1)
        """
        d50 = self.m_params['d50']          # sediment median diameter (m)
        Rous = self.m_params['rhoSed']      # sediment density (kg/m3)
        alphaA = self.m_params['alphaA']
        betaA = self.m_params['betaA']
        alphaD = self.m_params['alphaD']
        betaD = self.m_params['betaD']
        
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
        wup = Karman * np.sqrt(0.2*(ak**2)*(2*a0*nv)**(2/3)*np.abs(U)**(4/3)/Roul)
        wup[np.abs(U)<=utils.TOL] = 0.0
        alphaA_x = alphaA[pft]
        betaA_x = betaA[pft]
        alphaD_x = alphaD[pft]
        betaD_x = betaD[pft]
        aps = np.zeros_like(Bag)
        dps = np.zeros_like(Bag)
        cD = np.zeros_like(Bag)
        indice = np.logical_and(alphaA_x>0,Bag>0)
        aps[indice] = np.exp( np.log(alphaA_x[indice]) + betaA_x[indice]*np.log(Bag[indice]) )
        indice = np.logical_and(alphaD_x>0,Bag>0)
        dps[indice] = np.exp( np.log(alphaD_x[indice]) + betaD_x[indice]*np.log(Bag[indice]) )
        indice = np.logical_and(np.logical_and(aps>0, dps>0), np.abs(U)>0)
        cD[indice] = 2.0*(a0*nv/np.abs(U[indice])/dps[indice] + chi + \
          xi*0.25*np.pi*aps[indice]*dps[indice])
        indice = np.logical_and(Bag>utils.TOL, np.abs(U)>utils.TOL)
        wup[indice] = Karman*np.sqrt(0.2*(ak**2)*(U[indice]**2)* \
           (cD[indice]*aps[indice]*dps[indice])**(2/3)/Roul)
        return np.maximum( ws-wup, 0.0 )

###############################################################################    
class M12MOD(MACMODSuper):
    """Realization of the Morris et al. (2012) mineral accretion model.

    Attributes:
        parameters : d50, rhoSed, porSed, E0, tauE_cr, ks, alphaSG, betaSG
    Constants:
        
    """
    
    # constructor
    def __init__(self, params):
        self.m_params = params
        self.m_update_Css = True
        
    def mineral_suspension(self, inputs):
        """"Calculate mineral suspension rate.
        Arguments:
            inputs : driving data for mineral suspension calculation
        Returns: mineral suspension rate (kg m-2 s-1)
        """
        E0 = self.m_params['E0']    # the reference erosion rate (mm/yr)
        tauE_cr = self.m_params['tauE_cr']  # critical shear stress (Pa)
        Rous = 2650.0               # sediment density (kg/m3)
        tau = inputs['tau']         # bottom shear stress (Pa)
        pft = inputs['pft']         # platform pft
        Esed = inputs['Esed']       # sediment erosion (kg/m2/s)
        
        Esed[:] = 0.0
        indice = np.logical_and( tau>tauE_cr, pft==1 )
        Esed[indice] = 1e-3 * E0 * Rous * ((tau[indice]-tauE_cr)/0.25) / 3.1536e7
        return Esed
    
    def mineral_deposition(self, inputs):
        """"Calculate mineral deposition rate.
        Arguments:
            inputs : driving data for mineral deposition calculation
        Returns: mineral deposition rate (kg m-2 s-1)
        """
        ks = self.m_params['ks']        # vegetation trapping efficiency (m3/s/kg)
        Dsed = inputs['Dsed']   # sediment deposition (kg/m2/s)
        Css = inputs['Css']     # sediment conc (kg/m3)
        tau = inputs['tau']     # bottom shear stress (Pa)
        Bag = inputs['Bag']     # aboveground biomass (kg/m2)
        pft = inputs['pft']     # platform pft
        Css0 = inputs['refCss'] # reference sediment conc (kg/m3)
        
        ws = self.settling_velocity(tau)
        Dsed[:] = 0.0
        indice = np.logical_and(np.logical_and(pft>=1,pft<=9), Css>Css0)
        Dsed[indice] = np.maximum( Css[indice]*(ws[indice]+ks[pft[indice]]* \
            Bag[indice]), 0.0 )
        return Dsed
    
#    def bed_loading(self, inputs):
#        """"Calculate sand bed loading rate.
#        Arguments:
#            inputs : driving data for bed loading calculation
#        Returns: bed loading rate (kg m-2 s-1)
#        """
#        alphaSG = self.m_params['alphaSG']  # m2/yr (default: 3.65)
#        betaSG = self.m_params['betaSG']    # m4/yr (default: 0.0019)
#        Rous = 2650.0               # sediment density (kg/m3)
#        Lbed = inputs['Lbed']   # sediment bed load (kg/m2/s)
#        pft = inputs['pft']     # platform pft
#        Bbg = inputs['Bbg']     # belowground biomass (kg/m2)
#        S = inputs['S']         # platform slope (m/m)
#
#        Lbed[:] = 0.0
#        Nx = len(pft)
#        for ii in range(Nx):
#            if pft[ii]>=1 and pft[ii]<=5:
#                Lbed[ii] = -Rous*max(alphaSG-betaSG*Bbg[ii],0.0)/3.1536e7*S[ii]
#                if ii<Nx-1:
#                    Sup = S[ii+1]
#                else:
#                    Sup = 0.0
#                if pft[ii+1]>=1 and pft[ii+1]<=5:
#                    Lbed[ii] = Lbed[ii] + Rous*max(alphaSG-betaSG*Bbg[ii],0.0)/ \
#                        3.1536e7*Sup
#        return Lbed
    
###############################################################################
class F07MOD(MACMODSuper):
    """Realization of the Fagherazzi et al. (2007) mineral accretion model.

    Attributes:
        parameters : d50, Rous, E0, tauE_cr, gamma, KD, tauD_cr
    Constants:
        
    """
    
    # constructor
    def __init__(self, params):
        self.m_params = params
        self.m_update_Css = True
        
    def mineral_suspension(self, inputs):
        """"Calculate mineral suspension rate.
        Arguments:
            inputs : driving data for mineral suspension calculation
        Returns: mineral suspension rate (kg m-2 s-1)
        """
        E0 = self.m_params['E0']    # reference erosion rate (kg/m2/s/Pa)
        tauE_cr = self.m_params['tauE_cr']  # critical shear stress (Pa)
        gamma = self.m_params['gamma']  # the increase of critical shear stress with depth
        d50 = self.m_params['d50']  # sediment median diameter (m)
        tau = inputs['tau']         # bottom shear stress (Pa)
        dtau = inputs['dtau']       # bottom shear stress diff (Pa)
        dt = inputs['dt']           # time step (s)
        pft = inputs['pft']         # platform pft
        Esed = inputs['Esed']       # sediment erosion (kg/m2/s)
        
        Esed[:] = 0.0
        if d50>6.25e-5: 
            indice = np.logical_and( tau>tauE_cr, pft==1 )
            Esed[indice] = E0*(tau[indice]-tauE_cr)**1.5
        else:
            indice = pft==1
            Esed[indice] = Esed[indice] + np.maximum(E0*(dtau[indice]- \
                gamma*Esed[indice]*dt), 0.0)
        return Esed
        
    def mineral_deposition(self, inputs):
        """"Calculate mineral deposition rate.
        Arguments:
            inputs : driving data for mineral deposition calculation
        Returns: mineral deposition rate (kg m-2 s-1)
        """
        tauD_cr = self.m_params['tauD_cr']  # critical shear stress (Pa)
        KD = self.m_params['KD']    # cohesive sed settling rate (m5 s-1 kg-4/3)
        d50 = self.m_params['d50']  # sediment median diameter (m)
        Css = inputs['Css']     # sediment conc (kg/m3)
        tau = inputs['tau']     # bottom shear stress (Pa)
        Css0 = inputs['refCss'] # reference sediment conc (kg/m3)
        pft = inputs['pft']     # platform pft
        Dsed = inputs['Dsed']   # sediment deposition (kg/m2/s)
        
        Dsed[:] = 0.0
        if d50>6.25e-5:
            ws = self.settling_velocity(tau)
            indice = np.logical_and(np.logical_and(pft>=1,pft<=9), Css>Css0)
            Dsed[indice] = np.maximum( Css[indice]*ws[indice], 0.0 )  
        else:
            indice = np.logical_and(np.logical_and(pft>=1,pft<=9), 
                                    np.logical_and(Css>Css0,tau<tauD_cr))
            Dsed[indice] = KD*Css[indice]**(7./3.)*(tauD_cr-tau[indice])/0.1
        return Dsed
    
###############################################################################
class VDK05MOD(MACMODSuper):
    """Realization of the van de Koppel et al. (2005) mineral accretion model.

    Attributes:
        parameters : Dmax, Emax, ds, aNv, bNv, d50, Rous
    Constants:
        
    """
    
    # constructor
    def __init__(self, params):
        self.m_params = params
        self.m_update_Css = True
        
    def mineral_suspension(self, inputs):
        """"Calculate mineral suspension rate.
        Arguments:
            inputs : driving data for mineral suspension calculation
        Returns: mineral suspension rate (kg m-2 s-1)
        """
        Emax = self.m_params['Emax']    # maximum sediment erosion by water flow (yr-1)
        ds = self.m_params['ds']        # a conversion coefficient for wave erosion (yr-1)
        aNv = self.m_params['aNv']      # sediment erosion increase rate by plant (kg/m2)
        bNv = self.m_params['bNv']      # sediment erosion increase rate by plant (kg/m2)
        Rous = 2650.0               # sediment density (kg/m3)
        Bag = inputs['Bag']         # aboveground biomass (kg/m2)
        tau = inputs['tau']         # bottom shear stress (Pa)
        h = inputs['h']             # water depth (m)
        S = inputs['S']             # platform surface slope (m/m)
        pft = inputs['pft']         # platform pft
        Esed = inputs['Esed']       # sediment erosion (kg/m2/s)
        
        Esed[:] = 0.0
        tau_max = np.max(tau)
        if tau_max>utils.TOL:
            # tide driven erosion
            indice = np.logical_and(pft>=1, pft<=9)
            Esed[indice] = Rous * Emax/3.1536e7 * (aNv/(aNv+Bag[indice])) * \
                (tau[indice]/tau_max)
            # wave driven erosion
            indice = np.logical_and(np.logical_and(pft>=1, pft<=9), h>utils.TOL)
            Esed[indice] = Esed[indice] + ds/3.1536e7 * (bNv/(bNv+Bag[indice])) * \
                np.maximum(S[indice],0.0)
        return Esed
        
    def mineral_deposition(self, inputs):
        """"Calculate mineral deposition rate.
        Arguments:
            inputs : driving data for mineral deposition calculation
        Returns: mineral deposition rate (kg m-2 s-1)
        """
        Dmax = self.m_params['Dmax']    # maximum deposition rate (m/yr)
        Rous = self.m_params['rhoSed']  # sediment density (kg/m3)
        Dsed = inputs['Dsed']           # sediment deposition (kg/m2/s)
        Css = inputs['Css']             # sediment conc (kg/m3)
        zh = inputs['zh']               # platform surface elevation (msl)
        pft = inputs['pft']             # platform pft
        Ks = inputs['MHHW']             # mean high high water level (msl)
        
        Dsed[:] = 0.0
        indice = np.logical_and(np.logical_and(zh>=0, zh<=Ks), 
                                np.logical_and(Css>utils.TOL,pft>=1))
        Dsed[indice] = Rous * Dmax/3.1536e7 * (1.0 - zh[indice]/Ks)
        return Dsed
    
###############################################################################
class DA07MOD(MACMODSuper):
    """Realization of the D'Alpaos et al. (2007) mineral accretion model.

    Attributes:
        parameters : E0, tauE_cr, tauD_cr, Kveg, Bmax, alphaN, betaN,
                     alphaH, betaH, alphaD, betaD, alphaE, betaE, gammaE
    Constants:
        
    """
    
    # constructor
    def __init__(self, params):
        self.m_params = params
        self.m_update_Css = True
        
    def mineral_suspension(self, inputs):
        """"Calculate mineral suspension rate.
        Arguments:
            inputs : driving data for mineral suspension calculation
        Returns: mineral suspension rate (kg m-2 s-1)
        """
        E0 = self.m_params['E0']    # reference erosion rate (mm/yr)
        tauE_cr0 = self.m_params['tauE_cr']  # critical shear stress (Pa)
        Kveg = self.m_params['Kveg']    # shear stress increase rate with Bag
        Bmax = self.m_params['Bmax']    # maximum Bag (kg/m2)
        Rous = 2650.0               # sediment density (kg/m3)
        tau = inputs['tau']         # bottom shear stress (Pa)
        Bag = inputs['Bag']         # aboveground biomass (kg/m2)
        pft = inputs['pft']         # platform pft
        Esed = inputs['Esed']       # sediment erosion (kg/m2/s)
        
        Esed[:] = 0.0
        Bmax_x = Bmax[pft]
        Kveg_x = Kveg[pft]
        indice = np.logical_and(np.logical_and(pft>=1,pft<=9), Bmax_x>utils.TOL)
        tauE_cr = tauE_cr0 * (1.0 + Kveg_x[indice]*Bag[indice]/Bmax_x[indice])
        Esed[indice] = np.maximum( 0.0, 1e-3*E0*Rous*(tau[indice]/tauE_cr-1.0)/3.1536e7 )
        indice = np.logical_and(np.logical_and(pft>=1,pft<=9), Bmax_x<=utils.TOL)
        Esed[indice] = np.maximum( 0.0, 1e-3*E0*Rous*(tau[indice]/tauE_cr0-1.0)/3.1536e7 )
        return Esed
        
    def mineral_deposition(self, inputs):
        """"Calculate mineral deposition rate.
        Arguments:
            inputs : driving data for mineral deposition calculation
        Returns: mineral deposition rate (kg m-2 s-1)
        """
        tauD_cr = self.m_params['tauD_cr']  # critical shear stress (Pa)
        alphaN = self.m_params['alphaN']    # coef for stem density per unit area
        betaN = self.m_params['betaN']      # coef for stem density per unit area
        alphaH = self.m_params['alphaH']    # coef for stem height
        betaH = self.m_params['betaH']      # coef for stem height
        alphaD = self.m_params['alphaD']    # coef for stem diameter
        betaD = self.m_params['betaD']      # coef for stem diameter
        alphaE = self.m_params['alphaE']    # coef for trapping efficiency
        betaE = self.m_params['betaE']      # coef for trapping efficiency
        gammaE = self.m_params['gammaE']    # coef for trapping efficiency
        d50 = self.m_params['d50']          # sediment median diameter (m)
        Css0 = inputs['refCss']     # reference sediment conc (kg/m3)
        Css = inputs['Css']         # sediment conc (kg/m3)
        U = inputs['U']             # water flow velocity (m/s)
        h = inputs['h']             # water depth (m)
        tau = inputs['tau']         # bottom shear stress (Pa)
        Bag = inputs['Bag']         # aboveground biomass (kg/m2)
        pft = inputs['pft']         # platform pft
        Dsed = inputs['Dsed']       # sediment deposition (kg/m2/s)
        
        Dsed[:] = 0.0
        nv = utils.visc
        # direct deposition
        ws = self.settling_velocity(tau)
        indice = np.logical_and(np.logical_and(pft>=1,pft<=9), 
                                np.logical_and(Css>Css0,tau<tauD_cr))
        Dsed[indice] = 2.0*ws[indice]*Css[indice]*(tauD_cr-tau[indice])/0.1
        # plant trapping
        alphaN_x = alphaN[pft]
        betaN_x = betaN[pft]
        alphaH_x = alphaH[pft]
        betaH_x = betaH[pft]
        alphaD_x = alphaD[pft]
        betaD_x = betaD[pft]
        ns = np.zeros_like(Bag)
        hs = np.zeros_like(Bag)
        ds = np.zeros_like(Bag)
        eps = np.zeros_like(Bag)
        indice = np.logical_and(alphaN_x>0, Bag>0)
        ns[indice] = np.exp( np.log(alphaN_x[indice]) + betaN_x[indice]*np.log(1e3*Bag[indice]))
        indice = np.logical_and(alphaH_x>0, Bag>0)
        hs[indice] = np.exp( np.log(alphaH_x[indice]) + betaH_x[indice]*np.log(1e3*Bag[indice]))
        indice = np.logical_and(alphaD_x>0, Bag>0)
        ds[indice] = np.exp( np.log(alphaD_x[indice]) + betaD_x[indice]*np.log(1e3*Bag[indice]))
        indice = ds>0
        eps[indice] = alphaE * (np.abs(U[indice])*ds[indice]/nv)**betaE * \
            (d50/ds[indice])**gammaE
        indice = np.logical_and(Bag>utils.TOL, Css>Css0)
        Dsed[indice] = Dsed[indice] + Css[indice]*np.abs(U[indice])* \
            eps[indice]*ds[indice]*ns[indice]*np.minimum(hs[indice],h[indice])
        return Dsed
    
#    def bed_loading(self, inputs):
#        """"Calculate sand bed loading rate (Zhou et al., 2016).
#        Arguments:
#            inputs : driving data for bed loading calculation
#        Returns: bed loading rate (kg m-2 s-1)
#        """
#        d50 = self.m_params['d50sand']      # sediment median diameter (m) default: 200e-6
#        Rous = 2650.0               # sediment density (kg/m3)
#        Lbed = inputs['Lbed']   # sediment bed load (kg/m2/s)
#        x = inputs['x']         # coordinate (m)
#        h = inputs['h']         # water depth (m)
#        U = inputs['U']         # water flow velocity (m/s)
#        Uwav = inputs['Uwav']   # wave speed (m/s)
#        G = utils.G
#        Roul = utils.Roul
#        Karman = utils.Karman
#        zr = 6e-3   # bed roughness length (m)
#        Fsand = np.zeros_like(h, dtype=np.float64)
#        indice = h >= 0.1
#        As = 0.005*h[indice]*(d50/h[indice])**1.2/((Rous/Roul-1)*G*d50)**1.2
#        # critical velocity for initiation of motion
#        Ucr = 0.19*(d50**0.1)*np.log(2*h[indice]/d50)   
#        Urms2 = 2.0 * Uwav**2  # wave orbital velocity
#        # non-dimensional drag coefficient due to current
#        Cdc = (Karman/(np.log(h[indice]/zr)-1))**2
#        Fsand[indice] = Rous*As*np.maximum(U[indice],0.0)* \
#            np.maximum((U[indice]**2+0.018/Cdc*Urms2[indice])**0.5-Ucr,0.0)**2.4
#        Nx = np.size(x)
#        for ii in range(Nx):
#            if ii==0:
#                Lbed[ii] = -Fsand[ii+1]/(x[ii+1]-x[ii])
#            elif ii==Nx-1:
#                Lbed[ii] = Fsand[ii-1]/(x[ii]-x[ii-1])
#            else:
#                Lbed[ii] = (Fsand[ii-1]-Fsand[ii+1])/(x[ii+1]-x[ii-1])
#        return Lbed
        
        
