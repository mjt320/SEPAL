# -*- coding: utf-8 -*-
"""
Created on Wed Mar  3 09:41:24 2021

@author: mthripp1
"""

from abc import ABC, abstractmethod
import numpy as np

from dce import signal_models, relax

class water_ex_model(ABC):
    # water exchange model.
    # used to convert compartment pops and relaxation rates to exponential relaxation pops and values

    @abstractmethod        
    def r_components(self, r_compa, p_compa):
        pass
    
    def r_compa_to_s(self, s0, p_compa, r_compa, signal_model, k=1.):
    #Calculate the total signal based on relaxation components
    #this is the p-weighted average signal over each relaxation component    
    
        r_compo, p_compo = self.r_components(r_compa, p_compa)
        
        s = np.sum([ p * signal_model.signal(s0, r_compo[idx], k) for idx, p in enumerate(p_compo) ], 0)
    
        return s
    
    
class fxl(water_ex_model):
    
    def r_components(self, r_compa, p_compa):
        #r_compa is a dict of relaxations (1 per tissue compartment)
        #p_compa is a dict of populations (1 per tissue compartment)
        #r_compo is a list of relaxations, (1 per relaxation component)
        #p_compo is a list of relaxations, (1 per relaxation component)
        r_compo = [ relax.relaxation(
            r_1 = np.sum( [ p_compa[compartment]*r_compa[compartment].r_1  for compartment in r_compa.keys() ] ,0) ,
            r_2s = np.sum( [ p_compa[compartment]*r_compa[compartment].r_2s for compartment in r_compa.keys() ] ,0)
            ) ]

        p_compo = [ 1. ]
        return r_compo, p_compo
    
class nxl(water_ex_model):
    
    def r_components(self, r_compa, p_compa):
        r2s_mean = np.sum([p_compa[compartment]*r_compa[compartment].r_2s for compartment in r_compa.keys()] ,0)
        r_compo = [ relax.relaxation(r_1 = r_compa['b'].r_1, r_2s = r2s_mean) ,
                    relax.relaxation(r_1 = r_compa['e'].r_1, r_2s = r2s_mean) ,
                    relax.relaxation(r_1 = r_compa['i'].r_1, r_2s = r2s_mean) ]
        p_compo = [ p_compa['b'], p_compa['e'], p_compa['i'] ]
        return r_compo, p_compo
    
class ntexl(water_ex_model):
    
    def r_components(self, r_compa, p_compa):
        r2s_mean = np.sum([p_compa[compartment]*r_compa[compartment].r_2s for compartment in r_compa.keys()] ,0)

        r1_ev = p_compa['e']*r_compa['e'].r_1 + p_compa['i']*r_compa['i'].r_1
        p_ev = p_compa['e'] + p_compa['i']
        
        r_compo = [ relax.relaxation(r_compa['b'].r_1, r2s_mean) ,
                    relax.relaxation(r1_ev, r2s_mean) ]
        p_compo = [ p_compa['b'], p_ev ]
        return r_compo, p_compo
