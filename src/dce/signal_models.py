# -*- coding: utf-8 -*-
"""
Created on Wed Mar  3 19:14:33 2021

@author: mthripp1
"""

import numpy as np

class signal_model:
    def __init__(self, acq_pars):
        self.acq_pars = acq_pars
        
    def signal(self, s0, r):
        pass



class spgr(signal_model):
    
    def signal(self, s0, r):
        s = s0 * (((1.0-np.exp(-self.acq_pars['tr']*r['r1']))*np.sin(self.acq_pars['fa_rad'])) /
              (1.0-np.exp(-self.acq_pars['tr']*r['r1'])*np.cos(self.acq_pars['fa_rad'])) ) \
        * np.exp(-self.acq_pars['te']*r['r2s'])    
        
        return s
