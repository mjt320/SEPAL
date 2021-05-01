# -*- coding: utf-8 -*-
"""
Created on Tue Apr  6 10:40:31 2021

@author: mthripp1
"""

#TODO: time delay/shift

from abc import ABC, abstractmethod

import numpy as np
from scipy.interpolate import interp1d

class aif(ABC):
    
    @abstractmethod
    def c_ap(self):
        pass
 
    
class patient_specific(aif):
    
    def __init__(self, t_data, c_ap_data):
        self.t_data = t_data
        self.c_ap_data = c_ap_data
        self.c_ap_func=interp1d(t_data, c_ap_data, kind='quadratic', bounds_error=False, fill_value = c_ap_data[0])        
        
    def c_ap(self, t):        
        c_ap = self.c_ap_func(t) 
        return c_ap

    
class parker_like(aif):
    
    def __init__(self, hct, a1 = 0.809, a2 = 0.330,
                 t1 = 0.17046, t2 = 0.365,
                 sigma1 = 0.0563, sigma2 = 0.132,
                 s = 38.078, tau = 0.483,
                 alpha = 0, beta = 0, alpha2 = 1.050, beta2 = 0.1685):
        self.a1, self.a2, self.t1, self.t2 = a1, a2, t1, t2
        self.sigma1, self.sigma2 = sigma1, sigma2
        self.s, self.tau = s, tau
        self.alpha, self.alpha2, self.beta, self.beta2 = alpha, alpha2, beta, beta2
        self.hct = hct
    
    def c_ap(self, t):
        
        t_mins = t/60.
        
        c_ab = (self.a1 / (self.sigma1*np.sqrt(2.*np.pi))) * np.exp(-((t_mins-self.t1)**2)/(2.*self.sigma1**2)) + \
            (self.a2 / (self.sigma2*np.sqrt(2.*np.pi))) * np.exp(-((t_mins-self.t2)**2)/(2.*self.sigma2**2)) + \
            (self.alpha*np.exp(-self.beta*t_mins) + self.alpha2*np.exp(-self.beta2*t_mins)) / (1+np.exp(-self.s*(t_mins-self.tau)))
        
        c_ap = c_ab / (1 - self.hct)
        
        return c_ap
    

class parker(parker_like):    
    def __init__(self, hct):
        super().__init__(hct)
