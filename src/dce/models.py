# -*- coding: utf-8 -*-
"""
Created on Sat Nov  7 21:05:18 2020

@author: mthripp1
"""

import numpy as np

req_pars = {'patlak': ('vp', 'ps', 've'),
            }

def irf_patlak(t, pk_pars):
    n = t.size
    h_cp = np.zeros(n, dtype=float)
    h_cp[0] = 1. * pk_pars['vp']    
    h_e = np.ones(n, dtype=float) * ( pk_pars['ps'] / pk_pars['ve'] )
    h_e[0] = h_e[0]/2.
    return h_cp, h_e

irfs = {
    'patlak': irf_patlak
}