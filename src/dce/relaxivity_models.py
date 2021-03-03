# -*- coding: utf-8 -*-
"""
Created on Fri Jan 29 09:49:37 2021

@author: mthripp1
"""

#TODO: define relaxivity per compartment

class relaxivity_model:
    
    def __init__(self, rlxy):
        self.rlxy = rlxy
        
    def r_compa(self, r_0, c_compa):
        r_compa = {
            compartment: self.r_single(r_0[compartment], c)
            for (compartment, c) in c_compa.items()
            }
        return r_compa
    
    def r_single(self, r_0, c):
        pass


class c_to_r_linear(relaxivity_model):
    
    def r_single(self, r_0, c):
        r = {
            'r1' : r_0['r1']  + self.rlxy['r1']  * c,
            'r2s': r_0['r2s'] + self.rlxy['r2s'] * c
            }        
        
        return r