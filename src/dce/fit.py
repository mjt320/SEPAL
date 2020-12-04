# -*- coding: utf-8 -*-
"""
Created on Sat Oct 31 09:03:05 2020

@author: mthripp1
"""

import numpy as np
from scipy.optimize import root, minimize

from . import models


def c_to_r_linear(c, r_0, rlxy_pars):
    r = {
        'r1': r_0['r1'] + rlxy_pars['r1']*c,
        'r2s': r_0['r2s'] + rlxy_pars['r2s']*c
    }
    return r

# dictionary containing references to c_to_r functions
c_to_r = {
    'linear': c_to_r_linear
}

def r_to_s_spgr(acq_pars, s0, r):    
    s = s0 * (((1.0-np.exp(-acq_pars['tr']*r['r1']))*np.sin(acq_pars['fa_rad'])) /
              (1.0-np.exp(-acq_pars['tr']*r['r1'])*np.cos(acq_pars['fa_rad'])) ) \
        * np.exp(-acq_pars['te']*r['r2s'])    
    return s

#dictionary containing references to r_to_s functions
r_to_s = {
    'spgr': r_to_s_spgr
}

def s_to_e(s, base_idx):
    s_pre = np.mean(s[base_idx])
    e = 100.*((s - s_pre)/s_pre)
    return e

def e_to_c(e, r_0, rlxy_pars, acq_pars, c_to_r_model, r_to_s_model):
    res = root(lambda c: e - c_to_e(c, r_0, rlxy_pars, acq_pars, c_to_r_model, r_to_s_model),
               x0= 0., method='hybr', options={'maxfev': 1000, 'xtol': 1e-7})
    assert res.success, 'Root finding failed.'
    return min(res.x)

def c_to_e(c, r_0, rlxy_pars, acq_pars, c_to_r_model, r_to_s_model):
    r = c_to_r[c_to_r_model](c, r_0, rlxy_pars)
    s_pre = r_to_s[r_to_s_model](acq_pars, 1., r_0)
    s = r_to_s[r_to_s_model](acq_pars, 1., r)
    e = 100. * ((s - s_pre) / s_pre)
    return e

def pkp_to_c(t, c_p_aif, pk_pars, irf_model):
    dt = t[1]-t[0] #need to generalise, add interpolation
    n = t.size
    h_cp, h_e = models.irfs[irf_model](t, pk_pars)
    c_cp = dt * np.convolve(c_p_aif, h_cp, mode='full')[:n]
    c_e  = dt * np.convolve(c_p_aif, h_e  , mode='full')[:n]
    c_t = pk_pars['vp'] * c_cp + pk_pars['ve'] * c_e
    return c_cp, c_e, c_t

def c_to_r(c_b, c_e, r_0_cp, r_0_e, c_to_r_model, rlxy_pars):
    return r

def r_to_comps(r, p, water_ex_model):
    return r_comps, p_comps

def comps_to_s(r_comps, p_comps, acq_pars, signal_model):
    return s

def pkp_to_s(t, c_p_aif, pk_pars, s0, r_0, rlxy_pars, acq_pars, irf_model, c_to_r_model, r_to_s_model):
    c_cp, c_e, c_t = pkp_to_c(t, c_p_aif, pk_pars, irf_model)
    #r = c_to_r[c_to_r_model](c_t, r_0, rlxy_pars)
    #s = r_to_s[r_to_s_model](acq_pars, s0, r)
    r = c_to_r(c_b, c_e, r_0_cp, r_0_e, c_to_r_model, rlxy_pars)
    r_comps, p_comps = r_to_comps(r, p, water_ex_model)
    s = comps_to_s(r_comps, p_comps, acq_pars, signal_model)
    return s

def pkp_to_x(pk_pars, s0, irf_model):
    # create vector x listing values in pk_pars dict (in standard order), followed by s0
    x = np.asarray([pk_pars[par] for par in models.req_pars[irf_model]] + [s0], dtype=np.float32)
    return x

def x_to_pkp(x, irf_model):
    # convert vector x of parameter values to pk_pars dict and s0 value
    parnames_in_x = [parname for parname in models.req_pars[irf_model]] #list pk parameter names corresponding to elements of x
    pk_pars = {parname:x[n] for n, parname in enumerate(parnames_in_x)} # dictionary containing pk parameters from x
    s0 = x[-1]
    return pk_pars, s0

def obj_fun(x, x_sf, t_mask, t, s, c_p_aif, r_0, rlxy_pars, acq_pars, irf_model, c_to_r_model, r_to_s_model):
    pk_pars_try, s0_try = x_to_pkp(x * x_sf, irf_model)
    s_try = pkp_to_s(t, c_p_aif, pk_pars_try, s0_try, r_0, rlxy_pars, acq_pars, irf_model, c_to_r_model, r_to_s_model)
    ssq = np.sum(t_mask*((s_try - s)**2)) #add masking later
    return ssq

# TO DO:
# INTERPOLATE BEFORE CONVOLUTION
# CHECK SPEED IF PARAMETERS FIXED VIA BOUNDS
# ADD WATER EXCHANGE
# ADD LINEAR CONSTRAINTS
# ADD MULTISTART
def s_to_pkp(t, s, c_p_aif, r_0, rlxy_pars, acq_pars, irf_model, c_to_r_model, r_to_s_model, fit_opts):
    s0_0 = s[0] / r_to_s[r_to_s_model](acq_pars, 1., r_0)
    print(s0_0)
    x_0 = pkp_to_x(fit_opts['pk_pars_0'], s0_0, irf_model) 
    x_sf= x_0 # use initial values to scale variables    
    x_lb_norm = pkp_to_x(fit_opts['pk_pars_lb'], s0_0*0.5, irf_model) / x_sf
    x_ub_norm = pkp_to_x(fit_opts['pk_pars_ub'], s0_0*1.5, irf_model) / x_sf
    x_0_norm = x_0/x_sf
    
    bounds = list(zip(x_lb_norm, x_ub_norm))    
    res = minimize(obj_fun, x_0_norm, args=(x_sf, fit_opts['t_mask'], t, s, c_p_aif, r_0, rlxy_pars, acq_pars, irf_model, c_to_r_model, r_to_s_model),
             method='trust-constr', bounds=bounds)#method='trust-constr', bounds=bounds, constraints=models.pkp_constraints[irf_model])
    pk_pars_opt, s0_opt = x_to_pkp(res.x * x_sf, irf_model)
    s_fit = pkp_to_s(t, c_p_aif, pk_pars_opt, s0_opt, r_0, rlxy_pars, acq_pars, irf_model, c_to_r_model, r_to_s_model)
    return pk_pars_opt, s0_opt, s_fit


