# -*- coding: utf-8 -*-
"""
Created on Sat Nov  7 21:05:18 2020

@author: Michael Thrippleton

PURPOSE: Module containin pharmacokinetic models.
    Includes functions to generate total and compartmental CA concentrations.
    Includes dictionary of all model functions, to allow switching between models.
"""

import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import LinearConstraint


# Dictionary of parameters associated with each model (WIP)
req_pars = {'patlak': ('vp', 'ps', 've'), # (v_e required to determine local EES concentration)         
          }

# Dictionary of constraints for each model (WIP)
pkp_constraints = {'patlak': LinearConstraint([[1., 0., 0., 0.],[1., 0., 1., 0.]], [0., 0.], [1., 1.])     
                 }


class pk_model():

    def __init__(self, t, dt_interp_request, c_p_aif, hct):
        self.t = t        
        self.c_p_aif = c_p_aif
        self.hct = hct
        
        #interpolate time points and AIF
        self.dt_interp, self.t_interp = interpolate_time_series(dt_interp_request, t)
        #interpolate AIF
        f=interp1d(t, c_p_aif, kind='quadratic', bounds_error=False, fill_value = c_p_aif[0])
        self.c_p_aif_interp = f(self.t_interp) 
        self.n_interp = self.t_interp.size()
        self.n = self.t.size()
        

        
    def irf(self, pk_pars):
        pass    
        
    def conc(self, pk_pars):        
   
        h_cp, h_e = self.irf(pk_pars)
    
        # Do the convolutions, taking only results in the required range    
        c_cp_interp = self.dt_interp * np.convolve(self.c_p_aif_interp, h_cp , mode='full')[:self.n_interp]
        c_e_interp  = self.dt_interp * np.convolve(self.c_p_aif_interp, h_e  , mode='full')[:self.n_interp]
        
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
    
    def split_pk_pars(self, pk_pars):
        var_pars = np.array( pk_pars['vp'], pk_pars['ps'] )
        fixed_pars = { 've': pk_pars['ve'] }
        return var_pars, fixed_pars    
    
    def join_pk_pars(self, var_pars, fixed_pars):
        pk_pars = {**fixed_pars,
                   **{'vp': var_pars[0],
                      've': var_pars[1] }
                   }
        return pk_pars
    
    def irf(self, pk_pars):

        #calculate h_cp, i.e. delta function at time zero
        h_cp = np.zeros(self.n_interp, dtype=float)
        h_cp[0] = 1./self.dt_interp

        #calculate h_e
        h_e = np.ones(self.n_interp, dtype=float) * (1./60.) * ( pk_pars['ps'] / pk_pars['ve'] )
        h_e[0] = h_e[0]/2.
        
        return h_cp, h_e

    
def interpolate_time_series(dt_required, t):
    #interpolate time series t to evenly spaced values from dt/2 to max(t)
    max_t = np.max(t)
    n_interp = np.round(max_t/dt_required + 0.5).astype(int)
    dt_actual = max_t/(n_interp-0.5)
    t_interp = np.linspace(0.5*dt_actual, max_t, num=n_interp)
    return dt_actual, t_interp
