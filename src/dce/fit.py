# -*- coding: utf-8 -*-
"""
Created on Sat Oct 31 09:03:05 2020

@author: mthripp1
"""

import numpy as np
from scipy.optimize import root, minimize

from . import models
from .conc_to_r import c_to_r, c_to_r_single
from .r_to_s import r_to_s, r_compartments_to_components, r_components_to_s


def s_to_e(s, base_idx):
    s_pre = np.mean(s[base_idx])
    e = 100.*((s - s_pre)/s_pre)
    return e

def e_to_c_fxl(e, r_0, rlxy_pars, acq_pars, c_to_r_model, r_to_s_model):
    res = root(lambda c: e - c_to_e_fxl(c, r_0, rlxy_pars, acq_pars, c_to_r_model, r_to_s_model),
               x0= 0., method='hybr', options={'maxfev': 1000, 'xtol': 1e-7})
    assert res.success, 'Root finding failed.'
    return min(res.x)

def c_to_e_fxl(c_t, r_0, rlxy_pars, acq_pars, c_to_r_model, r_to_s_model):
    #convert concentration to enhancement assuming fxl
    r = c_to_r_single[c_to_r_model](c_t, r_0, rlxy_pars)
    s_pre = r_to_s[r_to_s_model](acq_pars, 1., r_0)
    s_post = r_to_s[r_to_s_model](acq_pars, 1., r)
    e = 100. * ((s_post - s_pre) / s_pre)
    return e

def pkp_to_s(t, c_p_aif, pk_pars, s0, r_0_tissue, r_0_blood, rlxy_pars, acq_pars, hct, irf_model, c_to_r_model, r_to_s_model, water_model):    
    p_compartments = { 'b': pk_pars['vp']/(1-hct), 'e': pk_pars['ve'], 'i': 1-pk_pars['vp']/(1-hct)-pk_pars['ve'] } # DEFINE VB/USE hct HERE        
    r_0_extravasc = { 'r1': (r_0_tissue['r1']-p_compartments['b']*r_0_blood['r1'])/(1-p_compartments['b']),
                     'r2s': r_0_tissue['r2s'] }
    r_0_compartments = { 'b': r_0_blood,
                         'e': r_0_extravasc,
                         'i': r_0_extravasc }
    c_compartments, c_t = models.pkp_to_c(t, c_p_aif, pk_pars, hct, irf_model)   #dict of lists    
    r_compartments = c_to_r(c_compartments, r_0_compartments, c_to_r_model, rlxy_pars)
    r_components, p_components = r_compartments_to_components[water_model](r_compartments, p_compartments)    
    s = r_components_to_s(s0, r_components, p_components, acq_pars, r_to_s_model)
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

def obj_fun(x, x_sf, t_mask, t, s, c_p_aif, r_0_tissue, r_0_blood, rlxy_pars, acq_pars, hct, irf_model, c_to_r_model, r_to_s_model, water_model):
    pk_pars_try, s0_try = x_to_pkp(x * x_sf, irf_model)
    s_try = pkp_to_s(t, c_p_aif, pk_pars_try, s0_try, r_0_tissue, r_0_blood, rlxy_pars, acq_pars, hct, irf_model, c_to_r_model, r_to_s_model, water_model)
    ssq = np.sum(t_mask*((s_try - s)**2)) #add masking later
    return ssq

# TO DO:
# SIMPLIFY T2 IMPLMENTATION
# ADD DATA POINT MASKING
# CHECK SPEED IF PARAMETERS FIXED VIA BOUNDS
# INTERPOLATE BEFORE CONVOLUTION
# ADD LINEAR CONSTRAINTS
# CONSTRAINT AND MODEL FLEXIBILITY
# ADD MULTISTART
def s_to_pkp(t, s, c_p_aif, r_0_tissue, r_0_blood, rlxy_pars, acq_pars, hct, irf_model, c_to_r_model, r_to_s_model, water_model='fxl', fit_opts=None):
    s0_0 = s[0] / r_to_s[r_to_s_model](acq_pars, 1., r_0_tissue)
    print(s0_0)
    x_0 = pkp_to_x(fit_opts['pk_pars_0'], s0_0, irf_model) 
    x_sf= x_0 # use initial values to scale variables    
    x_lb_norm = pkp_to_x(fit_opts['pk_pars_lb'], s0_0*0.5, irf_model) / x_sf
    x_ub_norm = pkp_to_x(fit_opts['pk_pars_ub'], s0_0*1.5, irf_model) / x_sf
    x_0_norm = x_0/x_sf    
    bounds = list(zip(x_lb_norm, x_ub_norm))    
    res = minimize(obj_fun, x_0_norm, args=(x_sf, fit_opts['t_mask'], t, s, c_p_aif, r_0_tissue, r_0_blood, rlxy_pars, acq_pars, hct, irf_model, c_to_r_model, r_to_s_model, water_model),
             method='trust-constr', bounds=bounds)#method='trust-constr', bounds=bounds, constraints=models.pkp_constraints[irf_model])
    pk_pars_opt, s0_opt = x_to_pkp(res.x * x_sf, irf_model)
    s_fit = pkp_to_s(t, c_p_aif, pk_pars_opt, s0_opt, r_0_tissue, r_0_blood, rlxy_pars, acq_pars, hct, irf_model, c_to_r_model, r_to_s_model, water_model)
    return pk_pars_opt, s0_opt, s_fit


