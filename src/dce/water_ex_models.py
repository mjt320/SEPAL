# -*- coding: utf-8 -*-
"""
Created on Wed Mar  3 09:41:24 2021

@author: mthripp1
"""

import numpy as np

from. import signal_models

class water_ex_model:
    
    def __init__(self):
        pass
        
    def r_compa_to_compo(self):
        pass
    
    def r_compa_to_s(self, s0, p_compa, r_compa, signal_model):
    #Calculate the total signal based on relaxation components
    #this is the p-weighted average signal over each relaxation component    
    
        r_compo, p_compo = self.r_compa_to_compo(r_compa, p_compa)
        
        s = np.sum([ p * signal_model.signal(s0, r_compo[idx]) for idx, p in enumerate(p_compo) ], 0)
    
        return s
    
    
class fxl(water_ex_model):
    
    def r_compa_to_compo(self, r_compa, p_compo):
        #r_compa is a dict with an entry for each of 3 compartments.
        # Each entry is a dict containing an entry for r1 and r2*
        # Each value is an array of relaxation rate values
        #
        #r_compo is a list of dicts, 1 per componet
        #each dict contains r1 and r2* entries
        #each entry is an array of relaxation rate values
        r_compo = [
            {
            'r1':  np.sum( [ self.p_compa[compartment]*self.r_compa[compartment]['r1']  for compartment in self.r.keys() ] ,0),
            'r2s': np.sum( [ self.p_compa[compartment]*self.r_compa[compartment]['r2s'] for compartment in self.r.keys() ] ,0)
            }
        ]
        p_compo = [ np.ones(len(r_compo[0]['r1'])) ]
        return r_compo, p_compo
    
class nxl(water_ex_model):
    pass


def r_compa_to_compo_nxl(r, p):
    #convert compartment r and p to single values for NXL
    r2s_mean = np.sum([p[compartment]*r[compartment]['r2s'] for compartment in r.keys()] ,0)
    r_components = [ {'r1': r['b']['r1'], 'r2s': r2s_mean},
                     {'r1': r['e']['r1'], 'r2s': r2s_mean},
                     {'r1': r['i']['r1'], 'r2s': r2s_mean} ]
    p_components = [ p['b'], p['e'], p['i'] ]
    return r_components, p_components

