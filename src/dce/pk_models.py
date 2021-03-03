# -*- coding: utf-8 -*-
"""
Created on Sat Nov  7 21:05:18 2020

@author: Michael Thrippleton

PURPOSE: Module containin pharmacokinetic models.
    Includes functions to generate total and compartmental CA concentrations.
    Includes dictionary of all model functions, to allow switching between models.
"""

import numpy as np
from scipy.optimize import LinearConstraint


# Dictionary of parameters associated with each model (WIP)
req_pars = {'patlak': ('vp', 'ps', 've'), # (v_e required to determine local EES concentration)         
            }


# Dictionary of constraints for each model (WIP)
pkp_constraints = {'patlak': LinearConstraint([[1., 0., 0., 0.],[1., 0., 1., 0.]], [0., 0.], [1., 1.])     
                   }

def irf_patlak(dt, n, pk_pars):
    """
    Get impulse response functions for the Patlak model.
    
    Parameters
    ----------
    dt : separation between time points (s)
    n : number of time points
    pk_pars : dict containing pharmacokinetic parameters
        {'ps': <PS (min^-1),
         'vp': <vP>,
         've': <ve> }

    Returns
    -------
    h_cp : 1D numpy array
        IRF for capillary plasma compartment.
    h_e : 1D numpy array
        IRF for EES compartment.

    """
    
    #calculate h_cp, i.e. delta function at time zero
    h_cp = np.zeros(n, dtype=float)
    h_cp[0] = 1./dt    
    #calculate h_e
    h_e = np.ones(n, dtype=float) * (1./60.) * ( pk_pars['ps'] / pk_pars['ve'] )
    h_e[0] = h_e[0]/2.
    
    return h_cp, h_e


# Dictinoary of IRF functions
irfs = {
    'patlak': irf_patlak
}


def pkp_to_c(t, t_interp, c_p_aif_interp, pk_pars, hct, irf_model, options=None):
    """
    Calculate compartment and tissue concentrations.
    
    Calculates concentrations by convolving the AIF with the IRF. Both are interpolated
    to maximise accuracy but concentrations are returned at the measured time points.    

    Parameters
    ----------
    t : 1D numpy array
        Measured time points (s).
    t_interp : 1D numpy array
        Time points following interpolation.
    c_p_aif_interp : 1D numpy array
        AIF plasma concentration at interpolated time points.
    pk_pars : dict
        Dict containing pharmacokinetic parameters.
    hct : float
        Haematocrit.
    irf_model : string
        Pharmacokinetic model.
    options : dict, optional
        Not implemented.

    Returns
    -------
    c : dict
        Dict of compartment concentrations at measured time points.
        Each value is a 1D numpy array:
            {'b': <capillary blood concentration>,
             'e': <EES concentration>,
             'i': <intracellular concentration> }
    c_t: 1D numpy array
        Tissue concentration.

    """    
    n = t.size
    n_interp = t_interp.size
    dt_interp = t_interp[1]-t_interp[0]
    
    h_cp, h_e = irfs[irf_model](dt_interp, n_interp, pk_pars)

    # Do the convolutions, taking only results in the required range    
    c_cp_interp = dt_interp * np.convolve(c_p_aif_interp, h_cp, mode='full')[:n_interp]
    c_e_interp  = dt_interp * np.convolve(c_p_aif_interp, h_e  , mode='full')[:n_interp]
    
    # Resample concentrations at the measured time points
    c_cp = np.interp(t, t_interp, c_cp_interp)
    c_e = np.interp(t, t_interp, c_e_interp)
    
    c_b = c_cp*(1-hct)
    c_t = pk_pars['vp'] * c_cp + pk_pars['ve'] * c_e
    c = {'b': c_b, 'e': c_e, 'i': np.zeros(n)}
    
    return  c, c_t

    

class pk_model():
    #model base class
    #TODO: change so that constructor interpolates and generates interped aif
    # pk_pars become a parameter of conc; only need this once
    # constructur takes dt
    def __init__(self, t, t_interp, c_p_aif_interp, pk_pars, hct):
        self.t = t
        self.t_interp = t_interp
        self.c_p_aif_interp = c_p_aif_interp
        self.pk_pars = pk_pars
        self.hct = hct
        
    def irf(self, dt, n):
        pass    
        
    def conc(self):
        
        n = self.t.size
        n_interp = self.t_interp.size
        dt_interp = self.t_interp[1]-self.t_interp[0]
        
        h_cp, h_e = self.irf(dt_interp, n_interp)
    
        # Do the convolutions, taking only results in the required range    
        c_cp_interp = dt_interp * np.convolve(self.c_p_aif_interp, h_cp, mode='full')[:n_interp]
        c_e_interp  = dt_interp * np.convolve(self.c_p_aif_interp, h_e  , mode='full')[:n_interp]
        
        # Resample concentrations at the measured time points
        c_cp = np.interp(self.t, self.t_interp, c_cp_interp)
        c_e = np.interp(self.t, self.t_interp, c_e_interp)
        
        c_b = c_cp*(1-self.hct)
        c_t = self.pk_pars['vp'] * c_cp + self.pk_pars['ve'] * c_e
        c = {'b': c_b, 'e': c_e, 'i': np.zeros(n)}
        
        return  c, c_t

    
class patlak(pk_model):
    #Patlak model subclass
        
    def irf(self, dt, n):

        #calculate h_cp, i.e. delta function at time zero
        h_cp = np.zeros(n, dtype=float)
        h_cp[0] = 1./dt    

        #calculate h_e
        h_e = np.ones(n, dtype=float) * (1./60.) * ( self.pk_pars['ps'] / self.pk_pars['ve'] )
        h_e[0] = h_e[0]/2.
        
        return h_cp, h_e

    
    
    