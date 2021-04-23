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

    def __init__(self, t, dt_interp_request, aif, hct):
        self.t = t        
        self.aif = aif
        self.hct = hct
        
        #interpolate time points and AIF
        self.dt_interp, self.t_interp = interpolate_time_series(dt_interp_request, t)
        self.c_ap_interp = aif.c_ap(self.t_interp)
        self.n_interp = self.t_interp.size
        self.n = self.t.size
     
    @abstractmethod
    def var_pars(self, pk_pars):
        pass
    
    @abstractmethod
    def all_pars(self, all_pars):
        pass
    
    @abstractmethod
    def typicalx(self):
        pass
    
    @abstractmethod
    def irf(self, pk_pars):
        pass    
    
    @abstractmethod
    def constraints(self):
        pass
        
    def conc(self, pk_pars):        
      
        h_cp, h_e = self.irf(pk_pars)
    
        # Do the convolutions, taking only results in the required range    
        c_cp_interp = self.dt_interp * np.convolve(self.c_ap_interp, h_cp , mode='full')[:self.n_interp]
        c_e_interp  = self.dt_interp * np.convolve(self.c_ap_interp, h_e  , mode='full')[:self.n_interp]
        
        # Resample concentrations at the measured time points
        c_cp = np.interp(self.t, self.t_interp, c_cp_interp)
        c_e  = np.interp(self.t, self.t_interp, c_e_interp)
        
        c_t = pk_pars['vp'] * c_cp + pk_pars['ve'] * c_e
        c = {
            'b': c_cp * (1-self.hct),
            'e': c_e,
            'i': np.zeros(self.n)
            }
        
        return  c, c_t


    
class patlak(pk_model):
    #Patlak model subclass
    
    def var_pars(self, pk_pars):
        # take a dict of parameters and return variable parameters as a list
        var_pars = [ pk_pars['vp'], pk_pars['ps'] ]
        return var_pars
    
    def all_pars(self, vp, ps):
        # take variable parameters and all parameters as a dict, making sure
        # they are consistent with the model
        pk_pars = {'vp': vp,
                   'ps': ps,
                   've': 1. - vp/(1.-self.hct),
                   'vi': 0.,
                   'vb': vp/(1.-self.hct)
                   }
        return pk_pars
    
    def typicalx(self):
        pk_typical = { 'vp': 0.1,
                       'ps': 1e-3 }
        x_typical = self.var_pars(pk_typical)
        return x_typical
    
    def irf(self, pk_pars):

        #calculate h_cp, i.e. delta function at time zero
        h_cp = np.zeros(self.n_interp, dtype=float)
        h_cp[0] = 1./self.dt_interp

        #calculate h_e
        h_e = np.ones(self.n_interp, dtype=float) * (1./60.) * ( pk_pars['ps'] / pk_pars['ve'] )
        h_e[0] = h_e[0]/2.
        
        return h_cp, h_e
    
    def constraints(self):
        return ()


    
def interpolate_time_series(dt_required, t):
    #interpolate time series t to evenly spaced values from dt/2 to max(t)
    max_t = np.max(t)
    n_interp = np.round(max_t/dt_required + 0.5).astype(int)
    dt_actual = max_t/(n_interp-0.5)
    t_interp = np.linspace(0.5*dt_actual, max_t, num=n_interp)
    return dt_actual, t_interp
