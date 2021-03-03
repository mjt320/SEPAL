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

class c_to_r_model:
    
    def __init__(self, r_0, rlxy):
        self.r_0 = r_0
        self.rlxy = rlxy
        
    def r(self, c):
        r = { compartment: self.r_single(c, self.r_0[compartment], self.rlxy_pars)
            for (compartment, c) in c.items() }
    
    def r_single(self, c):
        pass


class c_to_r_linear(c_to_r_model):
    
    def r(self, c):
        r = {
            'r1': self.r_0['r1'] + self.rlxy['r1'] * c,
            'r2s': self.r_0['r2s'] + self.rlxy['r2s'] * c
            }        
        
        return r