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

class TAIHydroMOD(object):
    """Class of the 1D hydrodynamics model for coastal wetland.

    Primitive state variables:
        hydro_states[:,:] : [h (m), hU (m^2/s), Hw (m), css (kg/m^3), 
                             sal (PSU or g/kg)] over the domain
        
    Other state variables:    
        x_arr[:]    : transect cell coordinate (m)
        zh_arr[:]   : transect cell elevation relative to MSL (m)
        Cf_arr[:]   : flow conductance ()
        pft_arr[:]  : plant function type
        U_arr[:]    : water flow velocity (m/s)
        Ew_arr[:]   : wave energy density per unit area (J/m^2)
        tau_arr[:]  : bottom shear stress (Pa)
        
    Misc variables:
        site_id     : ID
        site_lat    : latitude (rad)
        site_lon    : longitude (rad)
        
    Parameters:
        model_params: a set of free parameters
        
    """
    
    # class member variables
    site_id = None
    site_lat = None
    site_lon = None
    x_arr = None
    zh_arr = None
    Cf_arr = None
    pft_arr = None
    U_arr = None
    Ew_arr = None
    tau_arr = None
    hydro_states = None
    forcings = None
    
    # model parameters (free + empirical)
    model_params = {}
    
    # constructor
    def __init__(self, nid, lat, lon, x, zh, pft):
        self.site_id = nid
        self.site_lat = lat
        self.site_lon = lon
        self.x_arr = x
        self.zh_arr = zh
        self.pft_arr = pft
        self.Ew_arr = np.zeros_like(x, dtype=np.float64, order='F')
        self.U_arr = np.zeros_like(x, dtype=np.float64, order='F')
        self.tau_arr = np.zeros_like(x, dtype=np.float64, order='F')
        # initialize the hydrodynamic state variables
        self.hydro_states[0] = -self.zh_arr
        self.hydro_states[0][self.hydro_states[0]<0] = 0.0
        self.hydro_states[1] = np.zeros_like(x, dtype=np.float64, order='F')
        self.hydro_states[2] = np.zeros_like(x, dtype=np.float64, order='F')
        self.hydro_states[3] = np.zeros_like(x, dtype=np.float64, order='F')
        self.hydro_states[4] = np.zeros_like(x, dtype=np.float64, order='F')
        
#    def load_ancillary_data(self, arrU0, arrTR, arrh0, arrU10, arrT, 
#                            arrTa):
#        """load ancillary data, including hydrodynamic forcings
#        Arguments:
#            arrU0 : tide speed at the seaward boundary (m/s)
#            arrTR : tidal range (m)
#            arrh0 : sea water level at the seaward boundary (m)
#            arrU10: wind speed at 10-m height over the domain (m/s)
#            arrT  : wave period (s)
#            arrTa : air temperature (celsius)
#        """
#        self.m_ancillary = {}   # remove previous loadings
#        
#        self.m_ancillary['U0'] = arrU0
#        self.m_ancillary['TR'] = arrTR
#        self.m_ancillary['h0'] = arrh0
#        self.m_ancillary['U10'] = arrU10
#        self.m_ancillary['T'] = arrT
#        self.m_ancillary['Ta'] = arrTa
        
    def set_model_parameters(self, d50=20e-6):
        """set model free parameters
        Arguments:
            
        """
        # sediment particle median diameter (m)
        self.model_params['d50'] = d50
        # flow conductance when no vegetation (van Rijn's (1984) formula)
        self.model_params['C0'] = 1.0
        
    def update_ground_roughness(self, Bag):
        """update plant-induced surface roughness for water flow
        Arguments:
            Bag : aboveground biomass (kg/m^2)
        """
        cb = self.model_params['cb']
        C0 = self.model_params['C0']
        for ii, h in enumerate(self.hydro_states[0]):
            if h<0.1:
                self.Cf_arr[ii] = 1e-20
                continue
            asb = alpha_a*(Bag[ii]**beta_a)
            dsb = alpha_d*(Bag[ii]**beta_d)
            cD = cD0 + ScD*Bag[ii]
            self.Cf_arr[ii] = C0*np.sqrt(2.0/(cD*(cb**2)*asb*h+2.0*(1-asb*dsb)))
        
    def update_shear_stress(self, T):
        """update bottom shear stress
        Arguments:
            T : wave period
        """
        d50 = self.model_params['d50']
        indice = self.hydro_states[0] > 0
        self.tau_arr[self.hydro_states[0]<=0] = 0.0
        # bottom shear stress induced by currents
        fcurr = 0.24 / (np.log(12*self.hydro_states[0][indice]/2.5/d50))**2
        tau_curr = 0.125*roul*fcurr*(self.U_arr[indice])**2
        # bottom shear stress induced by wind
        Um = np.pi*self.hydro_states[2][indice]/T/np.sinh(karman*self.hydro_states[0][indice])
        fwave = 1.39*(6.0*Um*T/np.pi/d50)**(-0.52)
        tau_wave = 0.5*fwave*roul*Um**2
        self.tau_arr[indice] = tau_curr*(1+1.2*(tau_wave/(tau_curr+tau_wave))**3.2)

    def wave_number(self, T):
        """calculate wave number
        Arguments:
            T : wave period
        """
        indice = self.hydro_states[0] > 0
        sigma = 2.0*np.pi/T
        k = np.ones_like(self.hydro_states[0][indice])
        for ii, h in enumerate(self.hydro_states[0][indice]):
            a = sigma**2/constants.g
            b = sigma / np.sqrt(constants.g*h)
            k[ii] = optimize.root_scalar(lambda x: constants.g*x*np.tanh(x*h)-sigma**2, \
                 bracket=[a,b], method='brentq').root
        return k

    def breaking_possibility(self):
        """calculate wave breaking possibility
        
        Arguments:

        """
        indice = self.hydro_states[0] > 0
        Hmax = 0.78 * self.hydro_states[0][indice]
        Qb = np.ones_like(Hmax)
        for ii, h in enumerate(self.hydro_states[0][indice]):
            Hrms = self.hydro_states[2][indice][ii]
            Hmax = 0.78 * h
            Qb[ii] = optimize.root_scalar(lambda x: 1-x+np.log(x)*(Hrms/Hmax)**2, \
                  bracket=[0,1], method='brentq').root
        return Qb
        
    def wave_generation(self, T, k, U10):
        """calculate wave generation by wind (J/m^2)
        
        Arguments:
            T  : wave period
            k  : wave number
            U10: 10-m wind speed
        """
        indice = self.hydro_states[0] > 0
        sigma = 2.0*np.pi/T
        alpha = 80.0*(roua**2)*sigma*(Cd**2)*(U10**2)/(roul**2)/ \
            (constants.g**2)/(k**2)
        beta = 5.0*roua/roul/T*(U10*k/sigma-0.9)
        Swg = np.zeros_like(self.x_arr)
        Swg[indice] = alpha + beta * self.Ew_arr[indice]
        return Swg
    
    def wave_bottom_friction(self, T, k, Qb):
        """calculate wave dissipation through bottom friction (J/m^2)
        
        Arguments:
            T  : wave period
            k  : wave number
            Qb : breaking possibility
        """
        indice = self.hydro_states[0] > 0
        Sbf = np.zeros_like(self.x_arr)
        Cf = 2*0.015*np.pi*self.hydro_states[2][indice]/T/ \
            np.sinh(k*self.hydro_states[0][indice])
        Sbf[indice] = (1-Qb)*2.0*Cf*k*self.Ew_arr[indice]/ \
            np.sinh(2*k*self.hydro_states[0][indice])
        return Sbf
        
    def wave_white_capping(self, T):
        """calculate wave dissipation through white capping (J/m^2)
        
        Arguments:
            T : wave period
        """
        indice = self.hydro_states[0] > 0
        Swc = np.zeros_like(self.x_arr)
        sigma = 2.0*np.pi/T
        gamma = self.Ew_arr[indice]*sigma**4/constants.g**2
        gamma_PM = 4.57e-3
        Swc[indice] = 1.0/3.0*1e-5*sigma*((gamma/gamma_PM)**2)*self.Ew_arr[indice]
        return Swc
        
    def wave_depth_breaking(self, T, k, U10, Qb):
        """calculate wave dissipation through depth-induced breaking (J/m^2)
        
        Arguments:
            T  : wave period
            k  : wave number
            U10: 10-m wind speed
            Qb : breaking possibility
        """
        indice = self.hydro_states[0] > 0
        Sbrk = np.zeros_like(self.x_arr)
        Hmax = 0.78 * self.hydro_states[0][indice]
        sigma = 2.0*np.pi/T
        alpha = 80.0*(roua**2)*sigma*(Cd**2)*(U10**2)/(roul**2)/ \
            (constants.g**2)/(k**2)
        Sbrk[indice] = 2*alpha/T*Qb*((Hmax/self.hydro_states[2][indice])**2)* \
            self.Ew_arr[indice]
        return Sbrk        
    
    # finite volume spatial discretization
    def odeFunc(self, uhydro_in, duhydro):
        """update hydrodynamic state variables
        
        Arguments:
            uhydro_in : the current state of hydrodynamics
            duhydro   : hydrodynamics time derivative
        """
        # for wave energy, it assumes that the maximum wave energy allowed is
        # roul*g*(fr*h)**2/8
        
        # the vegetation canopy on the marsh surface attenuate wind waves
        # Hred(%) = 3*Bag/Bmax*Latt (Hred is the % reduction of m_Hw)
        
        # at depths smaller than 0.1 m the flow velocity and sediment transport
        # are taken toward to zero (linearly interpolated to zero)
        """Define ODE function on the RHS
        
        Arguments:
            Swave   : the generation of wave energy (J/m^2)
            Ero     : erosion rate (kg/m^2)
            Dep     : deposition rate (kg/m^2)
        """
            
        N = np.shape(uhydro)[1]
        for i in range(1,N):
            phi_i = self.slope_limiter(uhydro_in,i)
            phi_im = self.slope_limiter(uhydro_in,i-1)
            phi_ip = self.slope_limiter(uhydro_in,i+1)
            uleft_minus = uhydro[:,i-1] + 0.5*phi_im*(uhydro[:,i]-uhydro[:,i-1])
            uright_minus = uhydro[:,i] - 0.5*phi_i*(uhydro[:,i+1]-uhydro[:,i])
            uleft_plus = uhydro[:,i] + 0.5*phi_i*(uhydro[:,i+1]-uhydro[:,i])
            uright_plus = uhydro[:,i+1] - 0.5*phi_ip*(uhydro[:,i+2]-uhydro[:,i+1]) 
            
    
    
class FVSKT_4th_RK_solver(object):
    """Class of the 1D Finite Volume Semidiscrete Kurganov and Tadmor (KT) 
    central scheme.
        
    """
    
    uhydro = None
    odeFunc = None
    MAXITER = 100
    
    def __init__(self, uhydro, odeFunc):
        self.uhydro = uhydro
        self.odeFunc = odeFunc
        
    def slope_limiter(self, i):
        """Define slope limiter
        
        Arguments:
            i : cell index
        """
        # Superbee scheme (https://en.wikipedia.org/wiki/Flux_limiter)
        ri = (self.uhydro[:,i]-self.uhydro[:,i-1])/(self.uhydro[:,i+1]-self.uhydro[:,i])
        return np.array([max(0,min(2*rij,1),min(rij,2)) for rij in ri])
        
    def fluxρ(self):
        return
    
    def flux(self, uhydro):
        
        return
    

    