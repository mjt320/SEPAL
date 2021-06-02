# -*- coding: utf-8 -*-
"""
Created on Wed Jul  8 17:50:09 2020.

@author: Michael Thrippleton

Functions to fit MRI signal to obtain T1:
    fit_vfa_2_point: obtain T1 using analytical formula based on two images
    fit_vfa_linear: obtain T1 using linear regression
    fit_vfa_nonlinear: obtain T1 using non-linear least squares fit

"""

import numpy as np
from scipy.optimize import curve_fit

from dce import signal_models

def fit_vfa_2_point(s, fa_rad, tr):
    """Return T1 based on 2 SPGR signals with different FA but same TR and TE.
    """
    
    with np.errstate(divide='ignore', invalid='ignore'):
        sr=s[0]/s[1]
        t1 = tr / np.log( (sr*np.sin(fa_rad[1])*np.cos(fa_rad[0]) - np.sin(fa_rad[0])*np.cos(fa_rad[1])) / (sr*np.sin(fa_rad[1]) - np.sin(fa_rad[0])) )
        s0 = s[0] * ( (1-np.exp(-tr/t1)*np.cos(fa_rad[0])) / ((1-np.exp(-tr/t1))*np.sin(fa_rad[0])) )
    
    t1 = np.nan if ~np.isreal(t1) | (t1<=0) | np.isinf(t1) else t1
    s0 = np.nan if (s0<=0) | np.isinf(s0) else s0
   
    return s0, t1


def fit_vfa_linear(s, fa_rad, tr):
    """Return T1 based on VFA signals using linear regression.

    Parameters
    ----------
        s: np array of signals. Can be >=1 dimensional but last dimensions should
            correspond to acquisition with particular flip angle
        fa_rad: 1-D np array of flip angles
        tr: TR

    Returns
    -------
        t1: T1
        s0: signal corresponding to equilibrium magnetisation excited with 90
            degree pulse            
    """
    y = s / np.sin(fa_rad)
    x = s / np.tan(fa_rad)
    A = np.stack([x, np.ones(x.shape)], axis=1)
    slope, intercept = np.linalg.lstsq(A, y, rcond=None)[0]

    is_valid = (intercept > 0) and (0. < slope < 1.)
    t1, s0 = (-tr/np.log(slope), intercept/(1-slope)) if is_valid else (np.nan, np.nan)
    
    return s0, t1


def fit_vfa_nonlinear(s,fa_rad,tr):
    """Return T1 based on VFA signals using NLLS fitting.
    
    Parameters
    ----------
        s: np array of signals. Can be >=1 dimensional but last dimensions should
            correspond to acquisition with particular flip angle
        fa_rad: 1-D np array of flip angles
        tr: TR

    Returns
    -------
        t1: T1
        s0: signal corresponding to equilibrium magnetisation excited with 90
            degree pulse  
    """
    p_linear = np.array(fit_vfa_linear(s, fa_rad, tr)) #linear fit to obtain initial guess
    
    # calculate typical values for scaling
    p_typ = np.array([ s[0] / spgr_signal(1., 1., tr, fa_rad[0]), 1. ])

    # if linear result is valid, use this otherwise use typical values
    p0 = p_linear if (~np.isnan(p_linear[0]) & ~np.isnan(p_linear[1])) else p_typ
    p0_scaled = p0 / p_typ

    popt_scaled, pcov = curve_fit(lambda x_fa_rad, p_s0, p_t1: spgr_signal(p_s0 * p_typ[0], p_t1 * p_typ[1], tr, x_fa_rad), fa_rad, s,
                           p0=p0_scaled, sigma=None, absolute_sigma=False, check_finite=True, bounds=(1e-5, np.inf), method='trf', jac=None) 
    popt = p_typ * popt_scaled
    s0 = popt[0]
    t1 = popt[1]
    return s0, t1
        

def spgr_signal(s0, t1, tr, fa_rad):      
    """Return signal for SPGR sequence.

    Parameters 
    ----------
        s0 (float or float array): equilibrium signal
        t1 (float or float array): T1 value
        tr (float or float array): TR value
        fa_rad (float or float array): flip angle

    Returns
    -------
        float or float array: steady-state signal

    """
    s=s0 * (((1.0-np.exp(-tr/t1))*np.sin(fa_rad)) / (1.0-np.exp(-tr/t1)*np.cos(fa_rad)) )
    
    return s

