# -*- coding: utf-8 -*-
"""
Created on Sat Nov  7 21:05:18 2020

@author: Michael Thrippleton

PURPOSE: Module containin pharmacokinetic models.
    Includes functions to generate total and compartmental CA concentrations.
    Includes dictionary of all model functions, to allow switching between models.
"""

from abc import ABC, abstractmethod

import numpy as np
from scipy.optimize import LinearConstraint

from dce import aifs


class pk_model(ABC):
    # PK model containing time points, aif and hct
    # used to return concentrations for specified PK parameters
    
    PARAMETERS = None
    DEFAULT_TYPICAL_PARS = None
    DEFAULT_CONSTRAINTS = None

    def __init__(self, t, dt_interp_request, aif, hct):
        self.t = t        
        self.aif = aif
        self.hct = hct
        
        #interpolate time points and AIF
        self.dt_interp, self.t_interp = interpolate_time_series(dt_interp_request, t)
        self.c_ap_interp = aif.c_ap(self.t_interp)
        self.n_interp = self.t_interp.size
        self.n = self.t.size
        
        self.typical_pars = type(self).DEFAULT_TYPICAL_PARS
        self.constraints = type(self).DEFAULT_CONSTRAINTS

    def conc(self, *pk_pars):        

        h_cp, h_e = self.irf(*pk_pars)
    
        # Do the convolutions, taking only results in the required range    
        C_cp_interp = self.dt_interp * np.convolve(self.c_ap_interp, h_cp , mode='full')[:self.n_interp]
        C_e_interp  = self.dt_interp * np.convolve(self.c_ap_interp, h_e  , mode='full')[:self.n_interp]
        
        # Resample concentrations at the measured time points
        C_cp = np.interp(self.t, self.t_interp, C_cp_interp)
        C_e  = np.interp(self.t, self.t_interp, C_e_interp)
        
        c_t = C_cp + C_e
       
        return  c_t, C_cp, C_e

    @abstractmethod
    def irf(self):
        pass  
        
    def pkp_array(self, pkp_dict):
        return np.array([ pkp_dict[p] for p in type(self).PARAMETER_NAMES ])
    
    def pkp_dict(self, pkp_array):
        return dict(zip(type(self).PARAMETER_NAMES, pkp_array))
        pass



    
class patlak(pk_model):
    #Patlak model subclass
    
    PARAMETER_NAMES = ['vp', 'ps']
    DEFAULT_TYPICAL_PARS = np.array([ 0.1, 1.e-3])
    DEFAULT_CONSTRAINTS = ()
    
    def irf(self, vp, ps):

        #calculate irf for capillary plasma (delta function centred at time zero)
        irf_cp = np.zeros(self.n_interp, dtype=float)
        irf_cp[0] = vp / self.dt_interp

        #calculate irf for the EES (constant term)
        irf_e = np.ones(self.n_interp, dtype=float) * (1./60.) * ps
        irf_e[0] = irf_e[0]/2.
        
        return irf_cp, irf_e



    
def interpolate_time_series(dt_required, t):
    #interpolate time series t to evenly spaced values from dt/2 to max(t)
    max_t = np.max(t)
    n_interp = np.round(max_t/dt_required + 0.5).astype(int)
    dt_actual = max_t/(n_interp-0.5)
    t_interp = np.linspace(0.5*dt_actual, max_t, num=n_interp)
    return dt_actual, t_interp
