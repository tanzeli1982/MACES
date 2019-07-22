#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  4 14:39:22 2019

1D transect-based hydrodynamics model for coastal TAI ecogeomorphology 
algorithm comparison

https://en.wikipedia.org/wiki/Airy_wave_theory

@author: Zeli Tan
"""

import numpy as np
from scipy import constants
from scipy import optimize
from numpy import linalg

roul = 1e3          # density of water (kg/m^3)
roua = 1.225        # density of air (kg/m^3)
karman = 0.4        # von Karman's constant
Cd = 1.3e-3         # drag coefficient for wind at 10-m height
G = constants.g     # gravitational constant (m/s^2)

class TAIHydroMOD(object):
    """Class of the 1D hydrodynamics model for coastal wetland.

    Primitive state variables:
        hydro_states[:,:] : [h (m), hU (m^2/s), N (m), css (kg/m^3), 
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
    dx_arr = None
    zh_arr = None
    Cf_arr = None
    pft_arr = None
    U_arr = None
    Ew_arr = None
    Hw_arr = None
    tau_arr = None
    hydro_states = None
    forcings = {'T': None, 'Bag': None, 'U10': None, 'U0': None, 'h0': None, 
                'Tair': None}
    k = None
    Qb = None
    Swg = None
    Sbf = None
    Swc = None
    Sbrk = None
    
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
        self.Hw_arr = np.zeros_like(x, dtype=np.float64, order='F')
        self.U_arr = np.zeros_like(x, dtype=np.float64, order='F')
        self.tau_arr = np.zeros_like(x, dtype=np.float64, order='F')
        # initialize the hydrodynamic state variables
        self.hydro_states[0] = -self.zh_arr
        self.hydro_states[0][self.hydro_states[0]<0] = 0.0
        self.hydro_states[1] = np.zeros_like(x, dtype=np.float64, order='F')
        self.hydro_states[2] = np.zeros_like(x, dtype=np.float64, order='F')
        self.hydro_states[3] = np.zeros_like(x, dtype=np.float64, order='F')
        self.hydro_states[4] = np.zeros_like(x, dtype=np.float64, order='F')
        self.dx_arr = np.zeros_like(x)
        for ii, val in enumerate(x):
            if ii==0:
                self.dx_arr[ii] = 0.5*(x[ii]+x[ii+1])
            elif ii==len(x)-1:
                self.dx_arr[ii] = 0.5*(x[ii-1]+x[ii])
            else:
                self.dx_arr[ii] = 0.5*(x[ii+1]-x[ii-1])
        # internal varibles
        self.k = np.zeros_like(x, dtype=np.float64, order='F')
        self.Qb = np.ones_like(x, dtype=np.float64, order='F')
        self.Swg = np.zeros_like(x, dtype=np.float64, order='F')
        self.Sbf = np.zeros_like(x, dtype=np.float64, order='F')
        self.Swc = np.zeros_like(x, dtype=np.float64, order='F')
        self.Sbrk = np.zeros_like(x, dtype=np.float64, order='F')
        
    def set_model_parameters(self, d50=20e-6):
        """set model free parameters
        Arguments:
            d50 : sediment particle median diameter (m)
        """
        self.model_params['d50'] = d50
        # flow conductance when no vegetation (van Rijn's (1984) formula)
        self.model_params['C0'] = 1.0
        # sediment dispersoin coefficient (m2/s-1) (Maan et al., 2005)
        self.model_params['K'] = 1e2
        # wave dissipation coefficient
        self.model_params['cbc'] = 0.015 # should cbc and cb or 1/C correlated
        # the fraction of water depth allowed for maximum wave height (0.5-0.8)
        self.model_params['fr'] = 0.78
        
    def update_ground_roughness(self, Bag):
        """update plant-induced surface roughness for water flow
        Arguments:
            Bag : aboveground biomass (kg/m^2)
        """
        cb = self.model_params['cb']
        C0 = self.model_params['C0']
        alpha_a = self.model_params['alpha_a']
        alpha_d = self.model_params['alpha_d']
        beta_a = self.model_params['beta_a']
        beta_d = self.model_params['beta_d']
        cD0 = self.model_params['cD0']
        ScD = self.model_params['ScD']
        for ii, h in enumerate(self.hydro_states[0]):
            if h<0.1:
                self.Cf_arr[ii] = 1e-20
                continue
            asb = alpha_a*(Bag[ii]**beta_a)
            dsb = alpha_d*(Bag[ii]**beta_d)
            cD = cD0 + ScD*Bag[ii]
            self.Cf_arr[ii] = C0*np.sqrt(2.0/(cD*(cb**2)*asb*h+2.0*(1-asb*dsb)))
        
    def update_shear_stress(self):
        """update bottom shear stress
        """
        d50 = self.model_params['d50']
        T = self.forcings['T']
        indice = self.hydro_states[0] > 0
        # bottom shear stress induced by currents
        fcurr = 0.24 / (np.log(12*self.hydro_states[0][indice]/2.5/d50))**2
        tau_curr = 0.125*roul*fcurr*(self.U_arr[indice])**2
        # bottom shear stress induced by wind
        Um = np.pi*self.Hw_arr[indice]/T/np.sinh(karman*self.hydro_states[0][indice])
        fwave = 1.39*(6.0*Um*T/np.pi/d50)**(-0.52)
        tau_wave = 0.5*fwave*roul*Um**2
        self.tau_arr = 0
        self.tau_arr[indice] = tau_curr*(1+1.2*(tau_wave/(tau_curr+tau_wave))**3.2)

    def wave_number(self):
        """calculate wave number
        
        """
        T = self.forcings['T']
        indice = self.hydro_states[0] > 0
        sigma = 2.0*np.pi/T
        self.k = 0
        for ii, h in enumerate(self.hydro_states[0][indice]):
            a = sigma**2/G
            b = sigma / np.sqrt(G*h)
            self.k[ii] = optimize.root_scalar(lambda x: G*x*np.tanh(x*h)-sigma**2, \
                 bracket=[a,b], method='brentq').root

    def breaking_possibility(self):
        """calculate wave breaking possibility
        
        Arguments:

        """
        indice = self.hydro_states[0] > 0
        Hmax = 0.78 * self.hydro_states[0][indice]
        self.Qb = 0
        for ii, h in enumerate(self.hydro_states[0][indice]):
            Hrms = self.Hw_arr[indice][ii]
            Hmax = 0.78 * h
            self.Qb[ii] = optimize.root_scalar(lambda x: 1-x+np.log(x)*(Hrms/Hmax)**2, \
                  bracket=[0,1], method='brentq').root
        
    def wave_generation(self):
        """calculate wave generation by wind (J/m^2)
        
        Arguments:
            
        """
        T = self.forcings['T']
        U10 = self.forcings['U10']
        indice = self.hydro_states[0] > 0
        sigma = 2.0*np.pi/T
        alpha = 80.0*(roua**2)*sigma*(Cd**2)*(U10**2)/(roul**2)/(G**2)/ \
            (self.k[indice]**2)
        beta = 5.0*roua/roul/T*(U10*self.k[indice]/sigma-0.9)
        self.Swg = 0
        self.Swg[indice] = alpha + beta * self.Ew_arr[indice]
    
    def wave_bottom_friction(self):
        """calculate wave dissipation through bottom friction (J/m^2)
        
        As the vegetation canopy on the marsh surface attenuate wind waves
        Hred(%) = 3*Bag/Bmax*Latt (Hred is the % reduction of Hw), cbc should
        be related to veg cover.
        
        Arguments:

        """
        T = self.forcings['T']
        cbc = self.model_params['cbc']
        indice = self.hydro_states[0] > 0
        Cf = 2*cbc*np.pi*self.Hw_arr[indice]/T/ \
            np.sinh(self.k[indice]*self.hydro_states[0][indice])
        self.Sbf = 0
        self.Sbf[indice] = (1-self.Qb[indice])*2.0*Cf*self.k[indice]* \
            self.Ew_arr[indice]/np.sinh(2*self.k[indice]*self.hydro_states[0][indice])
        
    def wave_white_capping(self):
        """calculate wave dissipation through white capping (J/m^2)
        
        Arguments
        """
        T = self.forcings['T']
        indice = self.hydro_states[0] > 0
        sigma = 2.0*np.pi/T
        gamma = self.Ew_arr[indice]*sigma**4/G**2
        gamma_PM = 4.57e-3
        self.Swc = 0
        self.Swc[indice] = 1.0/3.0*1e-5*sigma*((gamma/gamma_PM)**2)*self.Ew_arr[indice]
        
    def wave_depth_breaking(self):
        """calculate wave dissipation through depth-induced breaking (J/m^2)
        
        Arguments:

        """
        T = self.forcings['T']
        U10 = self.forcings['U10']
        indice = self.hydro_states[0] > 0
        Hmax = 0.78 * self.hydro_states[0][indice]
        sigma = 2.0*np.pi/T
        alpha = 80.0*(roua**2)*sigma*(Cd**2)*(U10**2)/(roul**2)/(G**2)/ \
            (self.k[indice]**2)
        self.Sbrk = 0
        self.Sbrk[indice] = 2*alpha/T*self.Qb[indice]*((Hmax/ \
            self.Hw_arr[indice])**2)*self.Ew_arr[indice]

    def setup(self, T, U0, h0, U10, Tair, Bag, zh, Ero, Dep):
        """update the status of the non-status variables
        
        Arguments:
            T   : wave period (s)
            U0  : tide speed at the seaward boundary (m/s)
            h0  : sea water level at the seaward boundary (m)
            U10 : wind speed at 10-m height over the domain (m/s)
            Tair: air temperature (celsius)
            Bag : aboveground biomass (kg/m^2)
            zh  : platform elevation (m)
            Ero : sediment suspending rate (kg/m^2/s)
            Dep : sediment deposition rate (kg/m^2/s)
        """
        self.forcings['T'] = T
        self.forcings['U0'] = U0
        self.forcings['h0'] = h0
        self.forcings['U10'] = U10
        self.forcings['Tair'] = Tair
        self.forcings['Bag'] = Bag
        self.forcings['Ero'] = Ero
        self.forcings['Dep'] = Dep
        self.zh_arr = zh
        self.update_ground_roughness(Bag)   # update ground roughness
        self.wave_number()                  # update wave number
        self.breaking_possibility()         # update wave breaking possibility
        self.wave_generation()              # update wave energy source
        self.wave_bottom_friction()         # update wave energy sink by friciton
        self.wave_white_capping()           # update wave energy sink by white capping
        self.wave_depth_breaking()          # update wave energy sink by depth breaking   
        
    def callback(self, uhydro):
        """update the status of the state variables

        at depths smaller than 0.1 m the flow velocity and sediment transport
        are taken toward to zero (linearly interpolated to zero)
        
        Arguments:
            uhydro : the updated hydrodynamics
            dt     : time step (s)
        """
        self.hydro_states = uhydro
        # update current velocity
        for ii, h in enumerate(self.hydro_states[0]):
            self.U_arr[ii] = self.hydro_states[1,ii] / max(0.1,self.hydro_states[0,ii])
        # update wave energy and significant wave height
        T = self.forcings['T']
        sigma = 2*np.pi/T
        self.Ew_arr = sigma * uhydro[2,:]
        self.Hw_arr = np.sqrt(8*self.Ew_arr/G/roul)
        
    def dF_1st(self, uhydro):
        """calculate the 1st order space derivative of the state variables
        
        Arguments:
            uhydro  : instant hydrodynamic states
        """
        T = self.forcings['T']
        Cg = np.zeros_like(self.x_arr, dtype=np.float64, order='F')
        indice = uhydro[0] > 0
        sigma = 2*np.pi/T
        Cg[indice] = 0.5*sigma/self.k[indice]*(1+2*self.k[indice]*uhydro[0][indice]/ \
          np.sinh(2*self.k[indice]*uhydro[0][indice]))
        dF = np.zeros_like(uhydro, dtype=np.float64, order='F')
        dF[0,:] = uhydro[1,:]
        dF[1,:] = uhydro[1,:]**2/uhydro[0,:] + 0.5*G*uhydro[0,:]**2
        dF[2,:] = Cg * uhydro[2,:]
        dF[3,:] = uhydro[1,:]*uhydro[3,:]
        dF[4,:] = uhydro[1,:]*uhydro[4,:]
        return dF
    
    def dF_jacobian(self, uhydro):
        """return the maximum eigenvalue of the Jacobian of the 1st order 
        space derivative of the state variables
        
        Ignore the non-diagonal terms of Jacobi for N, Css and Cj
        
        Arguments:
            uhydro  : instant hydrodynamic states
        """
        T = self.forcings['T']
        sigma = 2*np.pi/T
        Cg = np.zeros_like(self.x_arr, dtype=np.float64, order='F')
        indice = uhydro[0] > 0
        Cg[indice] = 0.5*sigma/self.k[indice]*(1+2*self.k[indice]* \
            uhydro[0][indice]/np.sinh(2*self.k[indice]*uhydro[0][indice]))
        lamda = np.zeros_like(self.x_arr, dtype=np.float64, order='F')
        for ii, x in enumerate(self.x_arr):
            h = uhydro[0,ii]
            Uh = uhydro[1,ii]
            U = Uh / max(h,0.1)
            #Css = uhydro[3,ii]
            #Cj = uhydro[4,ii]
            #jacobi = np.array([[U,1,0,0,0],[U**2+G*h,2*U,0,0,0],
            #                   [0.5*(G/max(h,0.1))**0.5,0,Cg[ii],0,0],
            #                   [U*Css,Css,0,Uh,0],[U*Cj,Cj,0,0,Uh]])
            jacobi = np.array([[U,1,0,0,0],[U**2+G*h,2*U,0,0,0],
                               [0,0,Cg[ii],0,0],[0,0,0,Uh,0],[0,0,0,0,Uh]])
            lamda[ii] = np.max(np.abs(linalg.eigvals(jacobi)))
        return lamda
    
    def dP_1st(self, uhydro):
        """calculate the 1st order space derivative of the diffusive term
        
        Arguments:
            uhydro  : instant hydrodynamic states
        """
        K = self.model_params['K']
        dP = np.zeros_like(uhydro, dtype=np.float64, order='F')
        N = np.size(self.x_arr)
        for ii in range(1,N-1):
            dP[3,ii] = 0.5*K*(uhydro[0,ii]+uhydro[0,ii+1])* \
                (uhydro[3,ii+1]-uhydro[3,ii])/self.dx_arr[ii]
            dP[4,ii] = 0.5*K*(uhydro[0,ii]+uhydro[0,ii+1])* \
                (uhydro[4,ii+1]-uhydro[4,ii])/self.dx_arr[ii]
        return dP
    
    def Src_0rd(self, uhydro):
        """calculate the source term
        
        for wave energy, it assumes that the maximum wave energy allowed is
        roul*g*(fr*h)**2/8
        
        Arguments:
            uhydro  : instant hydrodynamic states
        """
        T = self.forcings['T']
        Ero = self.forcings['Ero']
        Dep = self.forcings['Dep']
        fr = self.model_params['fr']
        sigma = 2*np.pi/T
        Dep_inst = np.copy(Dep)
        indice = Dep_inst > 0
        Dep_inst[indice] = Dep_inst[indice]*uhydro[3,indice]/self.hydro_states[3,indice]
        Src = np.zeros_like(uhydro, dtype=np.float64, order='F')
        U = np.zeros_like(uhydro[1,:], dtype=np.float64, order='F')
        dB = np.zeros_like(uhydro[1,:], dtype=np.float64, order='F')
        for ii, h in enumerate(uhydro[0]):
            U[ii] = uhydro[1,ii] / max(0.1,uhydro[0,ii])
            if ii==0:
                dB[ii] = 0.5*(self.zh_arr[ii+1]-self.zh_arr[ii]) / self.dx_arr[ii]
            elif ii==len(self.x_arr)-1:
                dB[ii] = 0.5*(self.zh_arr[ii]-self.zh_arr[ii-1]) / self.dx_arr[ii]
            else:
                dB[ii] = 0.5*(self.zh_arr[ii+1]-self.zh_arr[ii-1]) / self.dx_arr[ii]
        Src[1,:] = -U*np.abs(U)/(self.Cf_arr)**2 - G*uhydro[0,:]*dB
        Src[2,:] = (self.Swg - self.Sbf - self.Swc - self.Sbrk) / sigma
        for ii, h in enumerate(uhydro[0]):
            Nmax = 0.125*roul*G*(fr*h)**2/sigma
            kx = self.k[ii]
            dx = self.dx_arr[ii]
            if uhydro[2,ii]>Nmax:
                Cg = 0.5*sigma/kx*(1+2*kx*h/np.sinh(2*kx*h))
                Src[2,ii] = Src[2,ii] - Cg*(Nmax-uhydro[2,ii])/dx
        Src[3,:] = Ero - Dep_inst
    
    # finite volume spatial discretization
    def odeFunc(self, uhydro, duhydro):
        """Define ODE function on the RHS
        
        Arguments:
            uhydro    : the state of hydrodynamics
            duhydro   : hydrodynamics time derivative
        """  
        N = np.size(self.x_arr)
        # calculate slope limiter
        uhydro_l = np.copy(uhydro)
        uhydro_r = np.copy(uhydro)
        for ii in range(N):
            phi = FVSKT.slope_limiter(uhydro,ii)
            if ii==0 or ii==N-1:
                uhydro_l[:,ii] = uhydro[:,ii]
                uhydro_r[:,ii] = uhydro[:,ii]
            else:
                uhydro_l[:,ii] = uhydro[:,ii] - 0.5*phi*(uhydro[:,ii+1]-uhydro[:,ii])
                uhydro_r[:,ii] = uhydro[:,ii] + 0.5*phi*(uhydro[:,ii+1]-uhydro[:,ii])
        F_l = self.dF_1st(uhydro_l)
        F_r = self.dF_1st(uhydro_r)
        Src = self.Src_0rd(uhydro)
        P = self.dP_1st(uhydro)
        F_minus = np.zeros_like(uhydro, dtype=np.float64, order='F')
        F_plus = np.zeros_like(uhydro, dtype=np.float64, order='F')
        aL = self.dF_jacobian(uhydro_l)
        aR = self.dF_jacobian(uhydro_r)
        for ii in range(1,N-1):
            ap = max(0.0, aR[ii], aL[ii+1])
            am = max(0.0, aR[ii-1], aL[ii])
            F_minus[:,ii] = 0.5*(F_r[:,ii-1]+F_l[:,ii]) - \
                0.5*am*(uhydro_l[:,ii]-uhydro_r[:,ii-1])
            F_plus[:,ii] = 0.5*(F_r[:,ii]+F_l[:,ii+1]) - \
                0.5*ap*(uhydro_l[:,ii+1]-uhydro_r[:,ii])
        for ii, dx in enumerate(self.dx_arr):
            if ii==0:
                # seaward boundary conditions
                duhydro[:,ii] = 0.0
            elif ii==N-1:
                # landward boundary conditions
                duhydro[:,ii] = 0.0
            else:
                duhydro[:,ii] = -(F_plus[:,ii] - F_minus[:,ii]) / dx + \
                    (P[:,ii] - P[:,ii-1]) / dx + Src[:,ii]  
    
class FVSKT(object):
    """Class of the 1D Finite Volume Semidiscrete Kurganov and Tadmor (KT) 
    central scheme.
        
    """
        
    @staticmethod
    def slope_limiter(uhydro, idx):
        """Define slope limiter
        
        Arguments:
            uhydro : the state of hydrodynamics
            idx    : cell index
        """
        # Superbee scheme (https://en.wikipedia.org/wiki/Flux_limiter)
        ri = (uhydro[:,idx]-uhydro[:,idx-1])/(uhydro[:,idx+1]-uhydro[:,idx])
        return np.array([max(0,min(2*rij,1),min(rij,2)) for rij in ri])
    

    