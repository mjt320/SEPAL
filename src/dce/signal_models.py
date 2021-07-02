# -*- coding: utf-8 -*-
"""
Created on Wed Mar  3 19:14:33 2021

@author: mthripp1
"""

from abc import ABC, abstractmethod
import numpy as np


class signal_model(ABC):
    # signal model. defined by acq parameters
    # used to return signal for specified s0 and relaxation parameters
    # should be vectorised
    
    def r_to_s(self, s0, p_compo, r_compo, k=1.):
    #Calculate the total signal based on relaxation components
    #this is the p-weighted average signal over each relaxation component    
        
        s = np.sum([ p * self.signal(s0, r_compo[idx], k) for idx, p in enumerate(p_compo) ], 0)
    
        return s    
    
    @abstractmethod
    def signal(self, s0, r, k):
        pass


class spgr(signal_model):
    
    def __init__(self, tr, fa_rad, te):
        self.tr = tr
        self.fa_rad = fa_rad
        self.te = te
    
    def signal(self, s0, r, k=1.):
        fa = k * self.fa_rad
        s = s0 * (((1.0-np.exp(-self.tr*r.r_1))*np.sin(fa)) /
              (1.0-np.exp(-self.tr*r.r_1)*np.cos(fa)) ) \
        * np.exp(-self.te*r.r_2s)    
        
        return s
