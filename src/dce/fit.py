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

def interpolate_time_series(dt_required, t):
    #interpolate time series t to evenly spaced values from dt/2 to max(t)
    max_t = np.max(t)
    n_interp = np.round(max_t/dt_required + 0.5).astype(int)
    dt_actual = max_t/(n_interp-0.5)
    t_interp = np.linspace(0.5*dt_actual, max_t, num=n_interp)
    return t_interp

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


def s_to_pkp(expt, data, proc):

    #interpolate time points and AIF
    t_interp=interpolate_time_series(proc.delta_t_interp, expt.t)

    #interpolate AIF
    f=interp1d(expt.t, data.c_p_aif, kind='quadratic', bounds_error=False, fill_value=data.c_p_aif[0])
    c_p_aif_interp=f(t_interp)

    # estimate s0, scale parameters
    s0_0 = data.s[0] / r_to_s[expt.r_to_s_model](expt.acq_pars, 1., data.r_0_tissue)
    x_0 = pkp_to_x(proc.fit_opts['pk_pars_0'], s0_0, proc.irf_model) 
    x_sf= x_0 # use initial values to scale variables    
    x_lb_norm = pkp_to_x(proc.fit_opts['pk_pars_lb'], s0_0*0.5, proc.irf_model) / x_sf
    x_ub_norm = pkp_to_x(proc.fit_opts['pk_pars_ub'], s0_0*1.5, proc.irf_model) / x_sf
    x_0_norm = x_0/x_sf    
    bounds = list(zip(x_lb_norm, x_ub_norm))   
    
    #perform fitting
    res = minimize(obj_fun, x_0_norm, args=(x_sf, expt, data, proc, t_interp),
             method='trust-constr', bounds=bounds)#method='trust-constr', bounds=bounds, constraints=models.pkp_constraints[irf_model])
    pk_pars_opt, s0_opt = x_to_pkp(res.x * x_sf, proc.irf_model)
    s_fit = pkp_to_s(t_interp, expt, data, proc, pk_pars_opt, s0_opt)
    s_fit[np.logical_not(proc.fit_opts['t_mask'])]=np.nan
    
    return pk_pars_opt, s0_opt, s_fit

def obj_fun(x, x_sf, t_interp, expt, data, proc, c_p_aif_interp):
    pk_pars_try, s0_try = x_to_pkp(x * x_sf, proc.irf_model)
    s_try = pkp_to_s(t_interp, c_p_aif_interp, pk_pars_try, s0_try, expt, data, proc)
    ssq = np.sum(proc.fit_opts['t_mask']*((s_try - data.s)**2))
    return ssq

def pkp_to_s(t_interp, c_p_aif_interp, pk_pars, s0_try, expt, data, proc):    
    p_compartments = { 'b': pk_pars['vp']/(1-data.hct), 'e': pk_pars['ve'], 'i': 1-pk_pars['vp']/(1-data.hct)-pk_pars['ve'] } # DEFINE VB/USE hct HERE        
    r_0_extravasc = { 'r1': (data.r_0_tissue['r1']-p_compartments['b']*data.r_0_blood['r1'])/(1-p_compartments['b']),
                     'r2s': data.r_0_tissue['r2s'] }
    r_0_compartments = { 'b': {'r1': data.r_0_blood['r1'], 'r2s': data.r_0_tissue['r2s']},
                         'e': r_0_extravasc,
                         'i': r_0_extravasc }
    
    #c_compartments, c_t = models.pkp_to_c(t, t_interp, c_p_aif_interp, pk_pars, hct, irf_model)   # pk model -> c
    pk_model = pk_model_type(expt.t, t_interp, c_p_aif_interp, pk_pars, proc.hct)
    c_compartments, c_t = pk_model.conc()
        
    #r_compartments = c_to_r(c_compartments, r_0_compartments, c_to_r_model, rlxy_pars) # cr model, r0_compartments
    
    #r_components, p_components = r_compartments_to_components[water_model](r_compartments, p_compartments) # water model, r_compartments
    
    #s = r_components_to_s(s0, r_components, p_components, acq_pars, r_to_s_model) #signal model, r_components, s0 -> s
    
    #return s

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


def interpolate_time_series(dt_required, t):
    #interpolate time series t to evenly spaced values from dt/2 to max(t)
    max_t = np.max(t)
    n_interp = np.round(max_t/dt_required + 0.5).astype(int)
    dt_actual = max_t/(n_interp-0.5)
    t_interp = np.linspace(0.5*dt_actual, max_t, num=n_interp)
    return t_interp

