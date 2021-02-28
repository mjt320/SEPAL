# -*- coding: utf-8 -*-
"""
Created on Sat Nov  7 21:05:18 2020

@author: mthripp1
"""

import numpy as np
from scipy.optimize import LinearConstraint

req_pars = {'patlak': ('vp', 'ps', 've'), # (v_e required to determine local EES concentration)
            }

pkp_constraints = {'patlak': LinearConstraint([[1., 0., 0., 0.],[1., 0., 1., 0.]], [0., 0.], [1., 1.])}

def irf_patlak(t, pk_pars):
    dt = t[1]-t[0] #to do: generalise this
    n = t.size
    h_cp = np.zeros(n, dtype=float)
    h_cp[0] = 1./dt
    h_e = np.ones(n, dtype=float) * (1./60.) * ( pk_pars['ps'] / pk_pars['ve'] )
    h_e[0] = h_e[0]/2.
    return h_cp, h_e

irfs = {
    'patlak': irf_patlak
}

def pkp_to_c(t, c_p_aif, pk_pars, hct, irf_model, options=None):
#version without interpolation
    dt = t[1]-t[0] #need to generalise, add interpolation
    n = t.size
    h_cp, h_e = irfs[irf_model](t, pk_pars)
    c_cp = dt * np.convolve(c_p_aif, h_cp, mode='full')[:n]
    c_e  = dt * np.convolve(c_p_aif, h_e  , mode='full')[:n]
    c_b = c_cp*(1-hct)
    c_t = pk_pars['vp'] * c_cp + pk_pars['ve'] * c_e
    c = {'b': c_b, 'e': c_e, 'i': np.zeros(n)}
    return  c, c_t

def pkp_to_c_2(t, t_interp, c_p_aif_interp, pk_pars, hct, irf_model, options=None):
#version with interpolation and conv
    dt = t_interp[1]-t_interp[0]
    n = t_interp.size
    h_cp, h_e = irfs[irf_model](t_interp, pk_pars)
    c_cp_interp = dt * np.convolve(c_p_aif_interp, h_cp, mode='full')[:n]
    c_e_interp  = dt * np.convolve(c_p_aif_interp, h_e  , mode='full')[:n]
    
    c_cp = np.interp(t, t_interp, c_cp_interp)
    c_e = np.interp(t, t_interp, c_e_interp)
    c_b = c_cp*(1-hct)
    c_t = pk_pars['vp'] * c_cp + pk_pars['ve'] * c_e
    c = {'b': c_b, 'e': c_e, 'i': np.zeros(n)}
    return  c, c_t
    
def pkp_to_c_3(t, t_interp, c_p_aif_interp, pk_pars, hct, irf_model, options=None):
#version with interpolation but no conv
#unless time series is highly interpolated this is slower than using np.convolve
#there is a loss of precision due to using the nearest interpolated time point
    dt_interp = t_interp[1]-t_interp[0]
    n_interp = t_interp.size
    n = t.size
    h_cp, h_e = irfs[irf_model](t_interp, pk_pars)
    c_cp = np.zeros(n)
    c_e= np.zeros(n)
    for i in range(n): #loop through time points
        i_interp = np.around(t[i]/dt_interp - 0.5).astype(int) #index of nearest interpolated time point - this is an approximation!
        c_cp[i] = dt_interp * np.sum( c_p_aif_interp[:i_interp+1] * h_cp[i_interp::-1] )
        c_e[i]  = dt_interp * np.sum( c_p_aif_interp[:i_interp+1] * h_e[i_interp::-1] )
    c_b = c_cp*(1-hct)
    c_t = pk_pars['vp'] * c_cp + pk_pars['ve'] * c_e
    c = {'b': c_b, 'e': c_e, 'i': np.zeros(n)}
    return  c, c_t 
    
def interpolate_time_series(dt_required, t):
    #interpolate time series t to evenly spaced values from dt/2 to max(t)
    max_t = np.max(t)
    n_interp = np.round(max_t/dt_required + 0.5).astype(int)
    dt_actual = max_t/(n_interp-0.5)
    t_interp = np.linspace(0.5*dt_actual, max_t, num=n_interp)
    return t_interp