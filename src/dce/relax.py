# -*- coding: utf-8 -*-
"""
Created on Fri Jan 29 09:49:37 2021

@author: mthripp1
"""

#TODO: define relaxivity per compartment

from abc import ABC, abstractmethod
from dataclasses import dataclass


import numpy as np

@dataclass
class relaxation:
    # each entry is a 1D np array of relaxation rates
    r_1: np.ndarray
    r_2s: np.ndarray
    

@dataclass
class relaxivity:
    # each entry is a relaxivity
    r_1: float
    r_2s: float
    

class c_to_r_model(ABC):
    # relaxivity model. defined by relaxivity parameters
    # r_single returns relaxation rate for a given initial relaxation rate and concentration
    # r_compa returns relaxation rate for each compartment
    def __init__(self, rlxy):
        self.rlxy = rlxy
        
    def r_compa(self, r0_compa, c_compa):
        # input: r_0 = dict of relaxatations, c_compa = dict of concentrations
        # output: r_compa = dict of relaxations
        r_compa = {
            compartment: self.r_single(r0_compa[compartment], c)
            for (compartment, c) in c_compa.items()
            }
        return r_compa
    
    @abstractmethod
    def r_single(self, r0, c):
        # input: r_0 = initial relaxation, c = 1D array of concentrations
        # output: r = relaxation
        pass


class c_to_r_linear(c_to_r_model):
    
    def r_single(self, r0, c):
        r = relaxation(r_1 = r0.r_1 + self.rlxy.r_1 * c,
                       r_2s = r0.r_2s + self.rlxy.r_2s * c)
        return r