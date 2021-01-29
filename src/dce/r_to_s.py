# -*- coding: utf-8 -*-
"""
Created on Fri Jan 29 10:09:02 2021

@author: mthripp1
"""

import numpy as np

# Start by defining r_to_s relationships for a single component/compartment

def r_to_s_spgr(acq_pars, s0, r):    
    s = s0 * (((1.0-np.exp(-acq_pars['tr']*r['r1']))*np.sin(acq_pars['fa_rad'])) /
              (1.0-np.exp(-acq_pars['tr']*r['r1'])*np.cos(acq_pars['fa_rad'])) ) \
        * np.exp(-acq_pars['te']*r['r2s'])    
    return s

r_to_s = { #dictionary containing references to r_to_s functions
    'spgr': r_to_s_spgr
}

## Define relationship between compartment relaxation and relaxation components, based on water exchange properties
def r_compartments_to_components_fxl(r, p):
    #convert compartment r and p to single values for FXL
    #in the FXL, p=1 and r is a p-weighted average of compartment r values
    #r is a dict (key=compartment) of dicts (key=r1/r2)
    #make average relaxation parameter dict
    r_comps = [
        {'r1': np.sum( [ p[comp]*r[comp]['r1'] for comp in r.keys() ] ,0),
         'r2s': np.sum( [ p[comp]*r[comp]['r1'] for comp in r.keys() ] ,0) }
        ]
    p_comps = [ np.ones(len(r_comps[0]['r1'])) ]
    return r_comps, p_comps

def r_compartments_to_components_nxl(r, p):
    #convert compartment r and p to single values for FXL
    #in the FXL, p=1 and r is a p-weighted average of compartment r values
    #r is a dict (key=compartment) of dicts (key=r1/r2)
    #make average relaxation parameter dict
    r_comps = [ r['b'], r['e'], r['i'] ]
    p_comps = [ p['b'], p['e'], p['i'] ]
    return r_comps, p_comps

r_compartments_to_components = { #dictionary containing references to r_to_comps functions
    'fxl': r_compartments_to_components_fxl,
    'nxl': r_compartments_to_components_nxl
    }

def r_components_to_s(s0, r_comps, p_comps, acq_pars, r_to_s_model):
    #Calculate the total signal based on relaxation components
    #this is thel
    #just p-weighted average signal over each relaxation component    
    s = np.sum([ p * r_to_s[r_to_s_model](acq_pars, s0, r_comps[idx]) for idx, p in enumerate(p_comps) ], 0)
    return s

