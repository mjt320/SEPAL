# -*- coding: utf-8 -*-
"""
Created on Sat Oct 31 09:03:05 2020

@author: mthripp1
"""

from dataclasses import dataclass

import numpy as np
from scipy.optimize import root, minimize
from scipy.interpolate import interp1d

from dce import relax

# TODO
# get conc calculation working
# fit concentration
# test with constraints
# ADD more MODELS
# ADD T1 FITTING
# TESTING - more jupyter notebook tests
# DOCUMENTATION
# DEAL WITH FIT WARNINGS
# EXCEPTIONS
# ADD MULTISTART

def s_to_e(s, base_idx):
    s_pre = np.mean(s[base_idx])
    e = 100.*((s - s_pre)/s_pre)
    return e

def e_to_c(e, k, r0, c_to_r_model, signal_model):    
    res = root(lambda c: e - c_to_e(c, k, r0, c_to_r_model, signal_model),
                x0= 0., method='hybr', options={'maxfev': 1000, 'xtol': 1e-7})
    assert res.success, 'Enh-to-conc root finding failed.'
    return min(res.x)

def c_to_e(c, k, r0, c_to_r_model, signal_model):
    #convert concentration to enhancement assuming fxl
    r = c_to_r_model.r_single(r0, c)
    s_pre = signal_model.signal(1., r0, k)
    s_post = signal_model.signal(1., r, k)
    e = 100. * ((s_post - s_pre) / s_pre)
    return e


def s_to_pkp(s, k, r0_tissue, r0_blood, pk_model, c_to_r_model, water_ex_model, signal_model, fit_opts=None):

    if 't_mask' not in fit_opts:
        fit_opts['t_mask'] = np.ones(s.shape)

    # estimate s0, scale parameters
    s0_0 = s[0] / signal_model.signal(1., r0_tissue)
    
    # list of variable parameters = s0 + variable PK parameters
    x_0 = np.array( [ s0_0, *pk_model.pars_asvect(fit_opts['pk_pars_0']) ] )

    x_scalefactor = [ s0_0, *pk_model.typicalx() ]
    
    x_0_norm = x_0 / x_scalefactor    
    #TODO: implement bounds and constraints
    #x_lb_norm = pkp_to_x(proc.fit_opts['pk_pars_lb'], s0_0*0.5, proc.irf_model) / x_sf
    #x_ub_norm = pkp_to_x(proc.fit_opts['pk_pars_ub'], s0_0*1.5, proc.irf_model) / x_sf
    #bounds = list(zip(x_lb_norm, x_ub_norm))   
    
    #define function to minimise
    def obj_fun(x, *args):
        s0_try = x[0] * x_scalefactor[0]
        pk_pars_try = pk_model.pars_asdict(fit_opts['pk_pars_0'], x[1:] * x_scalefactor[1:])
        s_try = pkp_to_s(pk_pars_try, s0_try, k, r0_tissue, r0_blood, pk_model, c_to_r_model, water_ex_model, signal_model)
        ssq = np.sum(fit_opts['t_mask']*((s_try - s)**2))    
        return ssq
    
    #perform fitting
    result = minimize(obj_fun, x_0_norm, args=None,
             method='trust-constr', bounds=None, constraints=pk_model.constraints())#method='trust-constr', bounds=bounds, constraints=models.pkp_constraints[irf_model])

    s0_opt = result.x[0] * x_scalefactor[0]
    pk_pars_opt = pk_model.pars_asdict(fit_opts['pk_pars_0'], result.x[1:] * x_scalefactor[1:])

    s_fit = pkp_to_s(pk_pars_opt, s0_opt, k, r0_tissue, r0_blood, pk_model, c_to_r_model, water_ex_model, signal_model)
    s_fit[np.logical_not(fit_opts['t_mask'])]=np.nan
    
    return pk_pars_opt, s0_opt, s_fit




def pkp_to_s(pk_pars, s0, k, r0_tissue, r0_blood, pk_model, c_to_r_model, water_ex_model, signal_model):   
    
    p_compa = {
        'b': pk_pars['vp']/(1-pk_model.hct),
        'e': pk_pars['ve'],
        'i': 1. - pk_pars['vp']/(1-pk_model.hct) - pk_pars['ve']
        } # DEFINE VB/USE hct HERE        
       
    r0_extravasc = relax.relaxation(
        r_1 = (r0_tissue.r_1-p_compa['b']*r0_blood.r_1)/(1-p_compa['b']),
        r_2s = r0_tissue.r_2s )
    
    r0_compa = { 'b': r0_blood,
                 'e': r0_extravasc,
                 'i': r0_extravasc }
   
    c_compa, c_t = pk_model.conc(pk_pars)    
    r_compa = c_to_r_model.r_compa(r0_compa, c_compa)        
    s = water_ex_model.r_compa_to_s(s0, p_compa, r_compa, signal_model, k)
    
    return s
