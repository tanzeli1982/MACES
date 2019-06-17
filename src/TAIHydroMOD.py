#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  4 14:39:22 2019

1D transect-based hydrodynamics model for coastal TAI ecogeomorphology 
algorithm comparison

@author: Zeli Tan
"""

import numpy as np
from scipy import constants
from scipy import optimize

roul = 1e3          # density of water (kg/m^3)
roua = 1.225        # density of air (kg/m^3)
karman = 0.4        # von Karman's constant
Cd = 1.3e-3         # drag coefficient for wind at 10-m height

class TAIHydroMOD:
    """Class of the 1D hydrodynamics model for coastal wetland.

    State variables:
        m_lat   : latitude (rad)
        m_lon   : longitude (rad)
        m_ncell : the number of transect cells
        m_d50   : sediment particle median diameter (m)
        m_x[:]  : transect cell coordinate (m)
        m_zh[:] : transect cell elevation relative to MSL (m)
        m_Cf[:] : flow conductance ()
        m_Bag[:]: aboveground biomass (kg/m^2)
        m_pft[:]: plant function type
        m_h[:]  : water flow depth (m)
        m_U[:]  : water flow velocity (m/s)
        m_E[:]  : wave energy density per unit area (J/m^2)
        m_Hw[:] : significant wave height (m)
        m_tau[:]: bottom shear stress (Pa)
        m_css[:]: suspended sediment concentration (kg/m^3)
        m_cj[:] : conceptual salinity or nutrients ()
        
    External drivings:
        m_ancillary: a set of ancillary data (wind, tide, temperature, ...)
        
    Parameters:
        m_params: a set of free parameters
        
    """
    
    # class member state variables
    m_ncell = 0
    m_lat = None
    m_lon = None
    m_d50 = 20e-6
    m_x = None
    m_zh = None
    m_Cf = None
    m_Bag = None
    m_pft = None
    m_h = None
    m_U = None
    m_E = None
    m_Hw = None
    m_tau = None
    m_css = None
    m_cj = None
    
    # class member ancillary data
    m_ancillary = {}
    
    # class free parameters
    m_params = {}
    
    # constructor
    def __init__(self, lat, lon, d50, x, zh):
        self.m_lat = lat
        self.m_lon = lon
        self.m_d50 = d50
        self.m_x = x
        self.m_zh = zh
        self.m_ncell = np.size(x)
        
    def load_ancillary_data(self, arrU0, arrTR, arrh0, arrU10, arrT, 
                            arrTa):
        """load ancillary data, including hydrodynamic forcings
        Arguments:
            arrU0 : tide speed at the seaward boundary (m/s)
            arrTR : tidal range (m)
            arrh0 : sea water level at the seaward boundary (m)
            arrU10: wind speed at 10-m height over the domain (m/s)
            arrT  : wave period (s)
            arrTa : air temperature (celsius)
        """
        self.m_ancillary = {}   # remove previous loadings
        
        self.m_ancillary['U0'] = arrU0
        self.m_ancillary['TR'] = arrTR
        self.m_ancillary['h0'] = arrh0
        self.m_ancillary['U10'] = arrU10
        self.m_ancillary['T'] = arrT
        self.m_ancillary['Ta'] = arrTa
    
    def init_state_variables(self, pft):
        """initialize the state variables
        Arguments:
            pft : platform initial plant function type
        """
        self.m_h = -self.m_zh
        self.m_h[self.m_h<0] = 0.0
        self.m_U = np.zeros(self.m_ncell)
        self.m_E = np.zeros(self.m_ncell)
        self.m_Hw = np.zeros(self.m_ncell)
        self.m_tau = np.zeros(self.m_ncell)
        self.m_css = np.zeros(self.m_ncell)
        self.m_cj = np.zeros(self.m_ncell)
        self.m_Bag = np.zeros(self.m_ncell)
        self.m_pft = pft
        
    def set_model_parameters(self, B_marsh, B_mangrove):
        """set model free parameters
        Arguments:
            
        """
        self.m_params['marsh']['aa'] = B_marsh['aa']
        self.m_params['marsh']['bb'] = B_marsh['bb']
        self.m_params['marsh']['cc'] = B_marsh['cc']
        self.m_params['mangrove']['aa'] = B_mangrove['aa']
        self.m_params['mangrove']['bb'] = B_mangrove['bb']
        self.m_params['mangrove']['cc'] = B_mangrove['cc']
        
    def update_ground_roughness(self, yindx, C0):
        """update plant-induced surface roughness for water flow
        Arguments:
            yindx : year index
            C0 : flow conductance when no vegetation (van Rijn's (1984) formula)
        """
        hMHT = self.m_ancillary['TR'][yindx]
        DMHT = hMHT - self.m_zh
        cb = self.m_params['cb']
        for ii in range(self.m_ncell):
            h = self.m_h[ii]
            if h<0.1:
                self.m_Bag[ii] = 1e20
                self.m_Cf[ii] = 1e-20
                continue
            pft = self.m_pft[ii]
            if pft in ['marsh','mangrove']:
                aa = self.m_params[pft]['aa']
                alpha_a = self.m_params[pft]['alpha_a']
                beta_a = self.m_params[pft]['beta_a']
                alpha_d = self.m_params[pft]['alpha_d']
                beta_d = self.m_params[pft]['beta_d']
                cD0 = self.m_params[pft]['cD0']
                ScD = self.m_params[pft]['ScD']
            else:
                aa, bb, cc = [0.0, 0.0, 0.0]
                alpha_a, beta_a, alpha_d, beta_d = [0.0, 0.0, 0.0, 0.0]
                cD0, ScD = [0.0, 0.0]
            self.m_Bag[ii] = aa*DMHT[ii] + bb*(DMHT[ii])**2 + cc
            asb = alpha_a*(self.m_Bag[ii]**beta_a)
            dsb = alpha_d*(self.m_Bag[ii]**beta_d)
            cD = cD0 + ScD*self.m_Bag[ii]
            self.m_Cf[ii] = C0*np.sqrt(2.0/(cD*(cb**2)*asb*h+2.0*(1-asb*dsb)))
        
    def update_shear_stress(self, ti):
        """update bottom shear stress
        Arguments:
            ti : time index
        """
        T = self.m_ancillary['T'][ti]
        indice = self.m_h > 0
        self.m_tau[self.m_h<=0] = 0.0
        # bottom shear stress induced by currents
        fcurr = 0.24 / (np.log(12*self.m_h[indice]/2.5/self.m_d50))**2
        tau_curr = 0.125*roul*fcurr*(self.m_U[indice])**2
        # bottom shear stress induced by wind
        Um = np.pi*self.m_Hw[indice]/T/np.sinh(karman*self.m_h[indice])
        fwave = 1.39*(6.0*Um*T/np.pi/self.m_d50)**(-0.52)
        tau_wave = 0.5*fwave*roul*Um**2
        self.m_tau[indice] = tau_curr*(1+1.2*(tau_wave/(tau_curr+tau_wave))**3.2)

    def wave_number(self, ti):
        """calculate wave number
        Arguments:
            ti : time index
        """
        T = self.m_ancillary['T'][ti]
        indice = self.m_h > 0
        sigma = 2.0*np.pi/T
        k = np.ones_like(self.m_h[indice])
        for ii, h in enumerate(self.m_h[indice]):
            a = sigma**2/constants.g
            b = sigma / np.sqrt(constants.g*h)
            k[ii] = optimize.root_scalar(lambda x: constants.g*x*np.tanh(x*h)-sigma**2, \
                 bracket=[a,b], method='brentq').root
        return k

    def breaking_possibility(self):
        """calculate wave breaking possibility
        Arguments:

        """
        indice = self.m_h > 0
        Hmax = 0.78 * self.m_h[indice]
        Qb = np.ones_like(Hmax)
        for ii, h in enumerate(self.m_h[indice]):
            Hrms = self.m_Hw[indice][ii]
            Hmax = 0.78 * h
            Qb[ii] = optimize.root_scalar(lambda x: 1-x+np.log(x)*(Hrms/Hmax)**2, \
                  bracket=[0,1], method='brentq').root
        return Qb
        
    def wave_generation(self, ti, k):
        """calculate wave generation by wind (J/m^2)
        Arguments:
            ti : time index
            k  : wave number
        """
        T = self.m_ancillary['T'][ti]
        U10 = self.m_ancillary['U10'][ti]
        indice = self.m_h > 0
        sigma = 2.0*np.pi/T
        alpha = 80.0*(roua**2)*sigma*(Cd**2)*(U10**2)/(roul**2)/ \
            (constants.g**2)/(k**2)
        beta = 5.0*roua/roul/T*(U10*k/sigma-0.9)
        Swg = np.zeros(self.m_ncell)
        Swg[indice] = alpha + beta * self.m_E[indice]
        return Swg
    
    def wave_bottom_friction(self, ti, k, Qb):
        """calculate wave dissipation through bottom friction (J/m^2)
        Arguments:
            ti : time index
            k  : wave number
            Qb : breaking possibility
        """
        T = self.m_ancillary['T'][ti]
        indice = self.m_h > 0
        Sbf = np.zeros(self.m_ncell)
        Cf = 2*0.015*np.pi*self.m_Hw[indice]/T/np.sinh(k*self.m_h[indice])
        Sbf[indice] = (1-Qb)*2.0*Cf*k*self.m_E[indice]/np.sinh(2*k*self.m_h[indice])
        return Sbf
        
    def wave_white_capping(self, ti):
        """calculate wave dissipation through white capping (J/m^2)
        Arguments:
            ti : time index
        """
        T = self.m_ancillary['T'][ti]
        indice = self.m_h > 0
        Swc = np.zeros(self.m_ncell)
        sigma = 2.0*np.pi/T
        gamma = self.m_E[indice]*sigma**4/constants.g**2
        gamma_PM = 4.57e-3
        Swc[indice] = 1.0/3.0*1e-5*sigma*((gamma/gamma_PM)**2)*self.m_E[indice]
        return Swc
        
    def wave_depth_breaking(self, ti, k, Qb):
        """calculate wave dissipation through depth-induced breaking (J/m^2)
        Arguments:
            ti : time index
            k  : wave number
            Qb : breaking possibility
        """
        T = self.m_ancillary['T'][ti]
        U10 = self.m_ancillary['U10'][ti]
        indice = self.m_h > 0
        Sbrk = np.zeros(self.m_ncell)
        Hmax = 0.78 * self.m_h[indice]
        sigma = 2.0*np.pi/T
        alpha = 80.0*(roua**2)*sigma*(Cd**2)*(U10**2)/(roul**2)/ \
            (constants.g**2)/(k**2)
        Sbrk[indice] = 2*alpha/T*Qb*((Hmax/self.m_Hw[indice])**2)*self.m_E[indice]
        return Sbrk
    
    # solve differential equations
    def solve_dynamic_equations(self):
        """update the aboveground biomass and surface roughness for water flow
        Arguments:
            
        """
        # for wave energy, it assumes that the maximum wave energy allowed is
        # roul*g*(fr*h)**2/8
        
        # the vegetation canopy on the marsh surface attenuate wind waves
        # Hred(%) = 3*Bag/Bmax*Latt (Hred is the % reduction of m_Hw)
        
        # at depths smaller than 0.1 m the flow velocity and sediment transport
        # are taken equal to zero
    
    # 
    
    