# -*- coding: utf-8 -*-
"""
Created on Fri Jan 29 09:49:37 2021

@author: mthripp1
"""

def c_to_r_linear(c, r_0, rlxy_pars):
    #convert concentration to relaxation for a single compartment
    r = {
        'r1': r_0['r1'] + rlxy_pars['r1']*c,
        'r2s': r_0['r2s'] + rlxy_pars['r2s']*c
    }
    return r

# dict of references to functions that convert concentration to relaxation (for a single compartment)
c_to_r_single = {
    'linear': c_to_r_linear
}

def c_to_r(c, r_0, c_to_r_model, rlxy_pars):
    #wrapper function to convert concentration to relaxation for all compartments
    #(assumes equal relaxivity parameters for all copmartments)
    r = { compartment: c_to_r_single[c_to_r_model](conc, r_0[compartment], rlxy_pars)
         for (compartment, conc) in c.items() }
    return r

