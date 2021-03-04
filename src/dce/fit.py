# -*- coding: utf-8 -*-
"""
Created on Sat Oct 31 09:03:05 2020

@author: mthripp1
"""

import numpy as np
from scipy.optimize import root, minimize
from scipy.interpolate import interp1d

from . import pk_models
from .conc_to_r import c_to_r, c_to_r_single
from .r_to_s import r_to_s, r_compartments_to_components, r_components_to_s

# TODO
# TIDY, USE CLASSES
# ADD MODELS
# TESTING
# CONSTRAINT AND MODEL FLEXIBILITY
# DEAL WITH FIT WARNINGS
# EXCEPTIONS
# ADD MULTISTART

class expt:
    def __init__(self, t, acq_pars, r_to_s_model):
        self.t = t
        self.acq_pars = acq_pars
        self.r_to_s_model = r_to_s_model

        
class proc:
    def __init__(self, irf_model, c_to_r_model, rlxy_pars, water_model, fit_opts=None, delta_t_interp=1.0):
        self.irf_model = irf_model
        self.c_to_r_model = c_to_r_model
        self.rlxy_pars = rlxy_pars
        self.water_model = water_model
        self.fit_opts = fit_opts
        self.delta_t_interp = delta_t_interp
       
        
class data:
    def __init__(self, s, c_p_aif, r_0_tissue, r_0_blood, hct):
        self.s = s
        self.c_p_aif = c_p_aif
        self.r_0_tissue = r_0_tissue
        self.r_0_blood = r_0_blood
        self.hct = hct


def s_to_pkp(s, acq_pars, proc):


    # estimate s0, scale parameters
    s0_0 = data.s[0] / r_to_s[proc.r_to_s_model](acq_pars, 1., data.r_0_tissue)
    
    x_0 = np.array( [ s0_0, *proc.pk_model.split_pk_pars(proc.pk_pars[0]) ] )

    # use initial values to scale variables    
    x_sf= x_0 
    x_0_norm = x_0/x_sf    
    #TODO: implement bounds and constraints
    #x_lb_norm = pkp_to_x(proc.fit_opts['pk_pars_lb'], s0_0*0.5, proc.irf_model) / x_sf
    #x_ub_norm = pkp_to_x(proc.fit_opts['pk_pars_ub'], s0_0*1.5, proc.irf_model) / x_sf
    #bounds = list(zip(x_lb_norm, x_ub_norm))   
    
    #perform fitting
    #TODO: implement multistart
    result = minimize(obj_fun, x_0_norm, args=(x_sf, s, proc),
             method='trust-constr', bounds=None)#method='trust-constr', bounds=bounds, constraints=models.pkp_constraints[irf_model])

    s0_opt = result.x[0] * x_sf[0]
    
    pk_pars_opt = proc.pk_model.join_pk_pars(proc.pk_pars_0, result.x[2:] * x_sf[2:])

    s_fit = pkp_to_s(pk_pars_opt, s0_opt, r_0_tissue, r_0_blood, proc)
    s_fit[np.logical_not(proc.fit_opts['t_mask'])]=np.nan
    
    return pk_pars_opt, s0_opt, s_fit

def obj_fun(x, x_sf, s, proc):
    
    s0_try = x[0] * x_sf[0]
    
    pk_pars_try = proc.pk_model.join_pk_pars(proc.pk_pars_0, x[2:] * x_sf[2:])
    
    s_try = pkp_to_s(pk_pars_try, s0_try, data, proc)
    
    ssq = np.sum(proc.fit_opts['t_mask']*((s_try - s)**2))
    
    return ssq


def pkp_to_s(pk_pars, s0, r_0_tissue, r_0_blood, pk_model, rlxy_model, water_ex_model, signal_model):   
    
    p_compa = {
        'b': pk_pars['vp']/(1-pk_model.hct),
        'e': pk_pars['ve'],
        'i': 1-pk_pars['vp']/(1-pk_model.hct)-pk_pars['ve']
        } # DEFINE VB/USE hct HERE        
    
    r_0_extravasc = { 'r1': (r_0_tissue['r1']-p_compa ['b']*r_0_blood['r1'])/(1-p_compa['b']),
                     'r2s': r_0_tissue['r2s'] }
    
    r_0_compa = { 'b': {'r1': r_0_blood['r1'], 'r2s': r_0_tissue['r2s']},
                         'e': r_0_extravasc,
                         'i': r_0_extravasc }
    
    c_compa, c_t = pk_model.conc(pk_pars)
    
    r_compa = rlxy_model.r_compa(r_0_compa, c_compa)
        
    s = water_ex_model.r_compa_to_s(s0, p_compa, r_compa, signal_model)
    
    return s


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


