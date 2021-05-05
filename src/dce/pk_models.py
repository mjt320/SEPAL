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

from dce import aifs


class pk_model(ABC):
    # PK model containing time points, aif
    # used to return concentrations for specified PK parameters
    
    PARAMETERS = None
    DEFAULT_TYPICAL_PARS = None
    DEFAULT_CONSTRAINTS = None

    def __init__(self, t, dt_interp_request, aif):
        self.t = t        
        self.aif = aif
        
        #interpolate time points and AIF
        self.dt_interp, self.t_interp = interpolate_time_series(dt_interp_request, t)
        self.c_ap_interp = aif.c_ap(self.t_interp)
        self.n_interp = self.t_interp.size
        self.n = self.t.size
        
        self.typical_pars = type(self).DEFAULT_TYPICAL_PARS
        self.constraints = type(self).DEFAULT_CONSTRAINTS

    def conc(self, *pk_pars, **pk_pars_kw):        

        h_cp, h_e = self.irf(*pk_pars, **pk_pars_kw)
    
        # Do the convolutions, taking only results in the required range    
        C_cp_interp = self.dt_interp * np.convolve(self.c_ap_interp, h_cp , mode='full')[:self.n_interp]
        C_e_interp  = self.dt_interp * np.convolve(self.c_ap_interp, h_e  , mode='full')[:self.n_interp]
        
        # Resample concentrations at the measured time points
        C_cp = np.interp(self.t, self.t_interp, C_cp_interp)
        C_e  = np.interp(self.t, self.t_interp, C_e_interp)
        
        C_t = C_cp + C_e
       
        return  C_t, C_cp, C_e

    @abstractmethod
    def irf(self):
        pass  
        
    def pkp_array(self, pkp_dict):
        return np.array([ pkp_dict[p] for p in type(self).PARAMETER_NAMES ])
    
    def pkp_dict(self, pkp_array):
        return dict(zip(type(self).PARAMETER_NAMES, pkp_array))
        pass


class steady_state_vp(pk_model):
    
    PARAMETER_NAMES = ('vp')
    DEFAULT_TYPICAL_PARS = np.array([ 0.1 ])
    DEFAULT_CONSTRAINTS = ()
    
    def irf(self, vp, **kwargs):

        #calculate irf for capillary plasma (delta function centred at time zero)
        irf_cp = np.zeros(self.n_interp, dtype=float)
        irf_cp[0] = vp / self.dt_interp

        #calculate irf for the EES (constant term)
        irf_e = np.zeros(self.n_interp, dtype=float)
        
        return irf_cp, irf_e    
    
class patlak(pk_model):
    #Patlak model subclass
    
    PARAMETER_NAMES = ('vp', 'ps')
    DEFAULT_TYPICAL_PARS = np.array([ 0.1, 1.e-3])
    DEFAULT_CONSTRAINTS = ()

       
    def irf(self, vp, ps, **kwargs):

        #calculate irf for capillary plasma (delta function centred at time zero)
        irf_cp = np.zeros(self.n_interp, dtype=float)
        irf_cp[0] = vp / self.dt_interp

        #calculate irf for the EES (constant term)
        irf_e = np.ones(self.n_interp, dtype=float) * (1./60.) * ps
        irf_e[0] = irf_e[0]/2.
        
        return irf_cp, irf_e

class extended_tofts(pk_model):
    
    PARAMETER_NAMES = ('vp', 'ps', 've')
    DEFAULT_TYPICAL_PARS = np.array([0.1, 1e-3, 0.2])
    DEFAULT_CONSTRAINTS = ()
    
    def irf(self, vp, ps, ve, **kwargs):
        
        #calculate irf for capillary plasma (delta function centred at time zero)
        irf_cp = np.zeros(self.n_interp, dtype=float)
        irf_cp[0] = vp / self.dt_interp

        #calculate irf for the EES
        irf_e = (1./60.) * ps * np.exp(-(self.t_interp * ps)/(60. * ve) )
        irf_e[0] = irf_e[0]/2.
        
        return irf_cp, irf_e

class tcxm(pk_model):
    
    PARAMETER_NAMES = ('vp', 'ps', 've', 'fp')
    DEFAULT_TYPICAL_PARS = np.array([0.1, 1e-3, 0.2, 20.])
    DEFAULT_CONSTRAINTS = ()
    
    def irf(self, vp, ps, ve, fp, **kwargs):

        fp_per_s = fp / (60. * 100.)
        ps_per_s = ps / 60.
        v = ve + vp
        T = v / fp_per_s
        tc = vp / fp_per_s
        te = ve / ps_per_s
        sig_p = ( (T + te) + np.sqrt((T + te)**2 - (4 * tc * te)) ) / (2 * tc * te)
        sig_n = ( (T + te) - np.sqrt((T + te)**2 - (4 * tc * te)) ) / (2 * tc * te)
                
        #calculate irf for capillary plasma
        irf_cp = vp * sig_p * sig_n * \
            ( (1 - te*sig_n) * np.exp(-self.t_interp*sig_n) + (te*sig_p - 1.) * np.exp(-self.t_interp*sig_p) ) \
                / ( sig_p - sig_n )
        irf_cp[0] /= 2.

        #calculate irf for the EES (constant term)
        irf_e = ve * sig_p * sig_n * ( np.exp(-self.t_interp*sig_n) - np.exp(-self.t_interp*sig_p) ) \
            / ( sig_p - sig_n )
        irf_e[0] /= 2.
        
        return irf_cp, irf_e
    
class tofts(pk_model):
    
    PARAMETER_NAMES = ('ktrans', 've')
    DEFAULT_TYPICAL_PARS = np.array([1e-2, 0.2])
    DEFAULT_CONSTRAINTS = ()
    
    def irf(self, ktrans, ve, **kwargs):
        
        #calculate irf for capillary plasma (zeros)
        irf_cp = np.zeros(self.n_interp, dtype=float)
        
        #calculate irf for the EES
        irf_e = ktrans * np.exp(-self.t_interp * ktrans/ve)
        
        return irf_cp, irf_e
    
    
def interpolate_time_series(dt_required, t):
    #interpolate time series t to evenly spaced values from dt/2 to max(t)
    max_t = np.max(t)
    n_interp = np.round(max_t/dt_required + 0.5).astype(int)
    dt_actual = max_t/(n_interp-0.5)
    t_interp = np.linspace(0.5*dt_actual, max_t, num=n_interp)
    return dt_actual, t_interp
