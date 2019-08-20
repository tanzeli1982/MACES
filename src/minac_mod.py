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
        d50 = self.m_params['d50']      # sediment median diameter (m)
        Rous = self.m_params['rhoSed']  # sediment density (kg/m3)
        pft = inputs['pft']     # platform pft
        Css = inputs['Css']     # sediment conc (kg/m3)
        tau = inputs['tau']     # bottom shear stress (Pa)
        Dsed = inputs['Dsed']   # sediment deposition (kg/m2/s)
        ws = self.settling_velocity(d50, Rous, tau)
        Dsed[:] = ws * Css
        Dsed[np.logical_or(Css<utils.TOL,pft==1)] = 0.0
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
        Dsed = inputs['Dsed']   # sediment deposition (kg/m2/s)
        pft = inputs['pft']     # platform pft
        H = inputs['zh']        # platform surface elevation (msl)
        Css = inputs['Css']     # sediment conc (kg/m3)
        Dc = inputs['x'] - inputs['xref']   # distance to wetland edge (m)
        Dsed[:] = k * np.exp(l*H) * np.exp(m*Dc) / 27.32 / 8.64e4
        Dsed[np.logical_or(H<0,Dc<0)] = 0.0
        Dsed[np.logical_or(Css<utils.TOL,pft==1)] = 0.0
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
        d50 = self.m_params['d50']      # sediment median diameter (m)
        Rous = self.m_params['rhoSed']  # sediment density (kg/m3)
        Dsed = inputs['Dsed']   # sediment deposition (kg/m2/s)
        pft = inputs['pft']     # platform pft
        Css = inputs['Css']     # sediment conc (kg/m3)
        tau = inputs['tau']     # bottom shear stress (Pa)
        Bag = inputs['Bag']     # aboveground biomass (kg/m2)
        U = inputs['U']         # flow velocity (m/s)
        ws = self.settling_velocity(d50, Rous, pft, tau, Bag, U)
        Dsed[:] = ws * Css
        Dsed[np.logical_or(Css<utils.TOL,pft==1)] = 0.0
        return Dsed
    
    def settling_velocity(self, d50, Rous, pft, tau, Bag, U):
        """"Calculate effective sediment settling velocity (Morris et al., 2012).
        Arguments:
            d50 : sediment median diameter (m)
            Rous : sediment density (kg/m3)
            pft : platform pft
            tau : bottom shear stress (Pa)
            Bag : aboveground biomass (kg/m2)
            U : tide flow velocity (m/s)
        Returns: sediment settling velocity (m s-1)
        """
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
        wup = np.zeros_like(Bag, dtype=np.float64, order='F')
        for ii, Bag_ii in enumerate(Bag):
            aps = alphaA[pft[ii]] * Bag_ii**betaA[pft[ii]]
            dps = alphaD[pft[ii]] * Bag_ii**betaD[pft[ii]]
            if Bag_ii>utils.TOL and abs(U[ii])>utils.TOL:
                cD = 2.0*(a0*nv/np.abs(U[ii])/dps + chi + xi*0.25*np.pi*aps*dps)
                wup[ii] = Karman * np.sqrt(0.2*(ak**2)*(U[ii]**2)* \
                   (cD*aps*dps)**(2/3)/Roul)
            elif abs(U[ii])<=utils.TOL:
                wup[ii] = 0.0
            else:
                wup[ii] = Karman * np.sqrt(0.2*(ak**2)*(2.0*a0*nv)**(2/3)* \
                   abs(U[ii])**(4/3)/Roul)
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
        
    def mineral_suspension(self, inputs):
        """"Calculate mineral suspension rate.
        Arguments:
            inputs : driving data for mineral suspension calculation
        Returns: mineral suspension rate (kg m-2 s-1)
        """
        E0 = self.m_params['E0']    # the reference erosion rate (mm/yr)
        tauE_cr = self.m_params['tauE_cr']  # critical shear stress (Pa)
        Rous = self.m_params['rhoSed']  # sediment density (kg/m3)
        tau = inputs['tau']         # bottom shear stress (Pa)
        pft = inputs['pft']         # platform pft
        Esed = inputs['Esed']       # sediment erosion (kg/m2/s)
        Esed[:] = 0.0
        indice = np.logical_and( tau>tauE_cr, pft==1 )
        Esed[indice] = 1e-3 * E0 * Rous * (tau[indice]/tauE_cr-1.0) / 3.1536e7
        return Esed
    
    def mineral_deposition(self, inputs):
        """"Calculate mineral deposition rate.
        Arguments:
            inputs : driving data for mineral deposition calculation
        Returns: mineral deposition rate (kg m-2 s-1)
        """
        ks = self.m_params['ks']        # vegetation trapping efficiency (m3/s/kg)
        d50 = self.m_params['d50']      # sediment median diameter (m)
        Rous = self.m_params['rhoSed']  # sediment density (kg/m3)
        Dsed = inputs['Dsed']   # sediment deposition (kg/m2/s)
        Css = inputs['Css']     # sediment conc (kg/m3)
        tau = inputs['tau']     # bottom shear stress (Pa)
        Bag = inputs['Bag']     # aboveground biomass (kg/m2)
        pft = inputs['pft']     # platform pft
        Css = inputs['Css']     # sediment conc (kg/m3)
        zh = inputs['zh']       # platform surface elevation (msl)
        TR = inputs['TR']       # tidal range (m)
        MHT = 0.5*TR            # mean high tide (msl)
        ws = self.settling_velocity(d50, Rous, tau)
        Dsed[:] = 0.0
        indice = np.logical_and(np.logical_and(zh>=0,zh<=MHT), 
                                np.logical_and(pft>=2,pft<=5))
        indice = np.logical_and(indice, Css>utils.TOL)
        DMHT = MHT - zh[indice]
        Dsed[indice] = np.maximum( Rous*Css[indice]*(ws[indice]+ \
            ks[pft[indice]]*Bag[indice])*(DMHT**2)/TR, 0.0 )
        return Dsed
    
    def bed_loading(self, inputs):
        """"Calculate sand bed loading rate.
        Arguments:
            inputs : driving data for bed loading calculation
        Returns: bed loading rate (kg m-2 s-1)
        """
        alphaSG = self.m_params['alphaSG']  # m2/yr
        betaSG = self.m_params['betaSG']    # m4/yr
        Rous = self.m_params['rhoSed']      # sediment density (kg/m3)
        Lbed = inputs['Lbed']   # sediment bed load (kg/m2/s)
        pft = inputs['pft']     # platform pft
        Bbg = inputs['Bbg']     # belowground biomass (kg/m2)
        S = inputs['S']         # platform slope (m/m)
        Lbed[:] = 0.0
        Nx = len(pft)
        for ii in range(Nx):
            if pft[ii]>=1 and pft[ii]<=5:
                Lbed[ii] = -Rous*max(alphaSG-betaSG*Bbg[ii],0.0)/3.1536e7*S[ii]
                if ii<Nx-1:
                    Sup = S[ii+1]
                else:
                    Sup = 0.0
                if pft[ii+1]>=1 and pft[ii+1]<=5:
                    Lbed[ii] = Lbed[ii] + Rous*max(alphaSG-betaSG*Bbg[ii],0.0)/ \
                        3.1536e7*Sup
        return Lbed
    
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
        
    def mineral_suspension(self, inputs):
        """"Calculate mineral suspension rate.
        Arguments:
            inputs : driving data for mineral suspension calculation
        Returns: mineral suspension rate (kg m-2 s-1)
        """
        E0 = self.m_params['E0']    # reference erosion rate (kg/m2/s/Pa)
        gamma = self.m_params['gamma']  # the increase of critical shear stress with depth
        dtau = inputs['dtau']       # bottom shear stress diff (Pa)
        dt = inputs['dt']           # time step (s)
        Esed = inputs['Esed']       # sediment erosion (kg/m2/s)
        Esed[:] = Esed + E0*(dtau-gamma*Esed*dt)
        Esed[Esed<0] = 0.0
        return Esed
        
    def mineral_deposition(self, inputs):
        """"Calculate mineral deposition rate.
        Arguments:
            inputs : driving data for mineral deposition calculation
        Returns: mineral deposition rate (kg m-2 s-1)
        """
        tauD_cr = self.m_params['tauD_cr']  # critical shear stress (Pa)
        KD = self.m_params['KD']    # cohesive sed settling rate (m5 s-1 kg-4/3)
        Css = inputs['Css']     # sediment conc (kg/m3)
        tau = inputs['tau']     # bottom shear stress (Pa)
        Dsed = inputs['Dsed']   # sediment deposition (kg/m2/s)
        Dsed[:] = np.maximum( KD*Css**(7/3)*(1-tau/tauD_cr), 0.0 )
        Dsed[:10] = 0.0   # avoid weird deposition at the seaward node
        return Dsed
    
    def bed_loading(self, inputs):
        """"Calculate sand bed loading rate.
        Arguments:
            inputs : driving data for bed loading calculation
        Returns: bed loading rate (kg m-2 s-1)
        """
        E0 = self.m_params['E0']    # reference erosion rate (kg/m2/s)
        tauE_cr = self.m_params['tauE_cr']  # critical shear stress (Pa)
        d50 = self.m_params['d50sand']  # sediment median diameter (m)
        Rous = self.m_params['rhoSed']  # sediment density (kg/m3)
        Lbed = inputs['Lbed']   # sediment bed load (kg/m2/s)
        Css = inputs['Css']     # sediment conc (kg/m3)
        tau = inputs['tau']     # bottom shear stress (Pa)
        ws = self.settling_velocity(d50, Rous, tau)
        # assume 50% of cohesive sediment
        Lbed[:] = 0.5*Css*ws
        indice = tau>tauE_cr
        Lbed[indice] = Lbed[indice] - E0*(tau[indice]-tauE_cr)**1.5
        return Lbed
    
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
        Rous = self.m_params['rhoSed']  # sediment density (kg/m3)
        Bag = inputs['Bag']         # aboveground biomass (kg/m2)
        tau = inputs['tau']         # bottom shear stress (Pa)
        h = inputs['h']             # water depth (m)
        zh = inputs['zh']           # platform surface elevation (m)
        S = inputs['S']             # platform surface slope (m/m)
        Esed = inputs['Esed']       # sediment erosion (kg/m2/s)
        tau_max = np.max(tau)
        Esed[:] = 0.0
        indice = np.logical_and(zh>=0, h>0)
        # tide driven erosion
        Esed[indice] = Rous * Emax/3.1536e7 * (aNv/(aNv+Bag[indice])) * \
            (tau[indice]/tau_max) * zh[indice]
        # wave driven erosion
        Esed[indice] = Esed[indice] + ds/3.1536e7 * (bNv/(bNv+Bag[indice])) * \
            S[indice] * zh[indice]
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
        Ks = 0.5*inputs['TR']           # mean high water level (msl)
        Dsed[:] = 0.0
        indice = np.logical_and(np.logical_and(zh>=0, zh<=Ks), Css>utils.TOL)
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
        Rous = self.m_params['rhoSed']  # sediment density (kg/m3)
        tau = inputs['tau']         # bottom shear stress (Pa)
        Bag = inputs['Bag']         # aboveground biomass (kg/m2)
        pft = inputs['pft']         # platform pft
        Esed = inputs['Esed']       # sediment erosion (kg/m2/s)
        Esed[:] = 0.0
        indice = Bmax[pft]>utils.TOL
        tauE_cr = tauE_cr0 * (1.0 + Kveg[pft[indice]]*Bag[indice]/Bmax[pft[indice]])
        Esed[indice] = np.maximum( 0.0, 1e-3*E0*Rous*(tau[indice]/tauE_cr-1.0)/3.1536e7 )
        indice = Bmax[pft]<=utils.TOL
        tauE_cr = tauE_cr0
        Esed[indice] = np.maximum( 0.0, 1e-3*E0*Rous*(tau[indice]/tauE_cr-1.0)/3.1536e7 )
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
        Rous = self.m_params['rhoSed']      # sediment density (kg/m3)
        Dsed = inputs['Dsed']       # sediment deposition (kg/m2/s)
        Css = inputs['Css']         # sediment conc (kg/m3)
        U = inputs['U']             # water flow velocity (m/s)
        h = inputs['h']             # water depth (m)
        tau = inputs['tau']         # bottom shear stress (Pa)
        Bag = inputs['Bag']         # aboveground biomass (kg/m2)
        pft = inputs['pft']         # platform pft
        nv = utils.visc
        # direct deposition
        ws = self.settling_velocity(d50, Rous, tau)
        Dsed[:] = np.maximum( 2.0*ws*Css*(1.0-tau/tauD_cr), 0.0 )
        Dsed[:10] = 0.0   # avoid weird deposition at the seaward node
        # plant trapping
        for ii, Bag_ii in enumerate(Bag):
            if Bag_ii>utils.TOL and Css[ii]>utils.TOL:
                ns = alphaN[pft[ii]] * Bag_ii**betaN[pft[ii]]
                hs = alphaH[pft[ii]] * Bag_ii**betaH[pft[ii]]
                ds = alphaD[pft[ii]] * Bag_ii**betaD[pft[ii]]
                eps = alphaE * (abs(U[ii])*ds/nv)**betaE * (d50/ds)**gammaE
                Dsed[ii] = Dsed[ii] + Css[ii]*abs(U[ii])*eps*ds*ns*min(hs,h[ii])
        return Dsed
    
#    def bed_loading(self, inputs):
#        """"Calculate sand bed loading rate (Zhou et al., 2016).
#        Arguments:
#            inputs : driving data for bed loading calculation
#        Returns: bed loading rate (kg m-2 s-1)
#        """
#        d50 = self.m_params['d50sand']      # sediment median diameter (m)
#        Rous = self.m_params['rhoSed']  # sediment density
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
        
        