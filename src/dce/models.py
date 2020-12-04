# -*- coding: utf-8 -*-
"""
Created on Sat Nov  7 21:05:18 2020

@author: mthripp1
"""

import numpy as np
from scipy.optimize import LinearConstraint

req_pars = {'patlak': ('vp', 'ps', 've'), # (v_e required to determine local EES concentration)
            }

pkp_constraints = {'patlak': LinearConstraint([[1., 0., 0., 0.],[1., 0., 1., 0.]], [0., 0.], [1., 1.])}

def irf_patlak(t, pk_pars):
    dt = t[1]-t[0] #to do: generalise this
    n = t.size
    h_cp = np.zeros(n, dtype=float)
    h_cp[0] = 1./dt
    h_e = np.ones(n, dtype=float) * (1./60.) * ( pk_pars['ps'] / pk_pars['ve'] )
    h_e[0] = h_e[0]/2.
    return h_cp, h_e

irfs = {
    'patlak': irf_patlak
}