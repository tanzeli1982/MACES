#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  5 21:53:49 2019

Base class for coastal wetland models (geomorphology + biogeochemistry)

@author: Zeli Tan
"""

from abc import ABCMeta, abstractmethod

class TAIMODSuper(object):
    """Abstract base class for TAI models.

    Attributes:
        sediment_arr[:]: sediment height
        biomass_arr[:]: plant biomass

    """
    
    sediment_arr = None
    biomass_arr = None
    
    __metaclass__ = ABCMeta
    
    @abstractmethod 
    def mineral_accretion(self, sed, shear):
        """"Calculate mineral accretion rate.
        Parameters:
            sed : suspended sediment concentration (g m-3)
            shear : water flow shear stress (Pa)
        Returns: mineral accretion rate (g m-2 s-1)
        """
        pass
    
    @abstractmethod 
    def organic_accretion(self, tide):
        """"Calculate organic matter accretion rate.
        Parameters:
            tide : inundation depth (m)
        Returns: organic matter accretion rate (g m-2 s-1)
        """
        pass

    @abstractmethod 
    def wind_erosion(self, wave):
        """"Calculate storm surge erosion rate.
        Parameters:
            wave : wave height
        Returns: storm surge erosion rate (m s-1)
        """
        pass  
    
    @abstractmethod 
    def landward_migration(self):
        """"Calculate coastal wetland landward migration.
        Parameters:
            
        Returns: coastal wetland landward migration
        """
        pass      