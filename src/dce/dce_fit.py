# -*- coding: utf-8 -*-
"""
Created on Sat Oct 31 09:03:05 2020

@author: mthripp1
"""

import numpy as np
from scipy.optimize import root, minimize, basinhopping, shgo

from dce import relax

def sig_to_enh(s, base_idx):
    s_pre = np.mean(s[base_idx])
    e = 100.*((s - s_pre)/s_pre)
    return e

def enh_to_conc(e, k, r0, c_to_r_model, signal_model):
    # convert enh to conc for one time point
    def enh_to_conc_single(enh):
        res = root(lambda c: enh - conc_to_enh(c, k, r0, c_to_r_model, signal_model),
                x0= 0., method='hybr', options={'maxfev': 1000, 'xtol': 1e-7})
        assert res.success, 'Enh-to-conc root finding failed.'
        return min(res.x)
    # apply to all time points
    c = np.asarray( [enh_to_conc_single(enh) for enh in e] )
    return c

def conc_to_enh(c, k, r0, c_to_r_model, signal_model):
    #convert concentration to enhancement assuming fxl
    r = c_to_r_model.r_single(r0, c)
    s_pre = signal_model.signal(1., r0, k)
    s_post = signal_model.signal(1., r, k)
    e = 100. * ((s_post - s_pre) / s_pre)
    return e

def conc_to_pkp(C_t, pk_model, pk_pars_0 = None, weights = None):
    if pk_pars_0 is None:
        pk_pars_0 = [pk_model.pkp_dict(pk_model.typical_vals)]
    if weights is None:
        weights = np.ones(C_t.shape)
        
    x_0_all = [pk_model.pkp_array(pars)  # get starting values as array
           for pars in pk_pars_0]
    x_scalefactor = pk_model.typical_vals
    x_0_norm_all = [x_0 / x_scalefactor for x_0 in x_0_all]
    
    #define function to minimise
    def cost(x_norm, *args):
        x = x_norm * x_scalefactor
        C_t_try, _C_cp, _C_e = pk_model.conc(*x)
        ssq = np.sum(weights * ((C_t_try - C_t)**2))
        return ssq
    
    #perform fitting
    result = minimize_global(cost, x_0_norm_all, args=None,
              bounds=None, constraints=pk_model.constraints, method='trust-constr')

    x_opt = result.x * x_scalefactor
    pk_pars_opt = pk_model.pkp_dict(x_opt)
    C_fit, _C_cp, _C_e = pk_model.conc(*x_opt)
    C_fit[weights == 0] = np.nan
    
    return pk_pars_opt, C_fit

def enh_to_pkp(enh, hct, k, r0_tissue, r0_blood, pk_model, c_to_r_model, water_ex_model, signal_model, pk_pars_0 = None, weights = None):
    
    if pk_pars_0 is None:
        pk_pars_0 = [pk_model.pkp_dict(pk_model.typical_vals)]
    if weights is None:
        weights = np.ones(enh.shape)
        
    x_0_all = [pk_model.pkp_array(pars)  # get starting values as array
           for pars in pk_pars_0]
    x_scalefactor = pk_model.typical_vals
    x_0_norm_all = [x_0 / x_scalefactor for x_0 in x_0_all]
    
    #define function to minimise
    def cost(x_norm, *args):
        x = x_norm * x_scalefactor
        pk_pars_try = pk_model.pkp_dict(x)
        enh_try = pkp_to_enh(pk_pars_try, hct, k, r0_tissue, r0_blood, pk_model, c_to_r_model, water_ex_model, signal_model)
        ssq = np.sum(weights * ((enh_try - enh)**2))    
        return ssq
    
    #perform fitting
    result = minimize_global(cost, x_0_norm_all, args=None,
             bounds=None, constraints=pk_model.constraints, method='trust-constr')

    x_opt = result.x * x_scalefactor
    pk_pars_opt = pk_model.pkp_dict(x_opt)
    enh_fit = pkp_to_enh(pk_pars_opt, hct, k, r0_tissue, r0_blood, pk_model, c_to_r_model, water_ex_model, signal_model)
    enh_fit[weights == 0]=np.nan
    
    return pk_pars_opt, enh_fit



def volume_fractions(pk_pars, hct):
    # if vp exists, calculate vb, otherwise set vb to zero
    if 'vp' in pk_pars:
        vb = pk_pars['vp'] / (1 - hct)
    else:
        vb = 0
    
    # if ve exists define vi as remaining volume, otherwise set to vi zero
    if 've' in pk_pars:
        ve = pk_pars['ve']
        vi = 1 - vb - ve
    else:
        ve = 1 - vb
        vi = 0    
    
    v = {'b': vb, 'e': ve, 'i': vi}
    return v
    


def pkp_to_enh(pk_pars, hct, k, r0_tissue, r0_blood, pk_model, c_to_r_model, water_ex_model, signal_model):   
   
    C_t, C_cp, C_e = pk_model.conc(**pk_pars)     
    v = volume_fractions(pk_pars, hct)    
    c = { 'b': C_cp / v['b'],
         'e': C_e / v['e'],
         'i': np.zeros(C_e.shape) }
    p = v
       
    r0_extravasc = relax.relaxation(
        r_1 = (r0_tissue.r_1-p['b']*r0_blood.r_1)/(1-p['b']),
        r_2s = r0_tissue.r_2s )
    
    r0 = { 'b': r0_blood,
           'e': r0_extravasc,
           'i': r0_extravasc }
   
    r = c_to_r_model.r_compa(r0, c)
    
    r_compo_0, p_compo_0 = water_ex_model.r_components(p, r0)
    r_compo, p_compo = water_ex_model.r_components(p, r)
    
    s_eq = 100.
    s0 = signal_model.r_to_s(s_eq, p_compo_0, r_compo_0, k)
    s = signal_model.r_to_s(s_eq, p_compo, r_compo, k)
    enh = 100. * (s - s0) / s0
    
    return enh





# HELPERS
def minimize_global(cost, x_0_all, **kwargs):
    results = [minimize(cost, x_0, **kwargs) for x_0 in x_0_all]
    costs = [result.fun for result in results]
    cost = min(costs)
    idx = costs.index(cost)
    result = results[idx]
    return result