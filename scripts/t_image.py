# -*- coding: utf-8 -*-
"""
Created on Mon Aug 17 15:43:41 2020

@author: Michael Thrippleton

Test script for mrifit.t1 (work in progress)

"""


import sys

import numpy as np
import nibabel as nib

sys.path.append("..")
from mrifit import t1, signal

#create 2D signal data with s0 and T1 varying in dimension 1 and 2 respectively
nx, ny, nz = 3, 4, 5
s0_map_true=np.tile( np.linspace(1,5,nx)[:,np.newaxis,np.newaxis] ,(1,ny,nz))
t1_map_true=np.tile( np.linspace(0.5,3,ny)[np.newaxis,:,np.newaxis] , (nx,1,nz))
fa1_rad, fa2_rad, fa3_rad = (2/360.0)*2*np.pi , (10/360.0)*2*np.pi , (18/360.0)*2*np.pi
fa_rad = np.array([fa1_rad, fa2_rad, fa3_rad])
tr_s=10e-3
s1=signal.spgr(s0_map_true, t1_map_true, tr_s, fa1_rad)
s2=signal.spgr(s0_map_true, t1_map_true, tr_s, fa2_rad)
s3=signal.spgr(s0_map_true, t1_map_true, tr_s, fa3_rad)
    
s_4d = np.stack((s1,s2,s3), axis=3)
s_1d=np.copy(s_4d[0,0,0,:])
s_2d=np.copy(s_4d[0,:,0,:])
s_3d=np.copy(s_4d[0,:,:,:])

s=s_1d

mask = np.atleast_1d(np.ones(s.shape[:-1]))


# s0_obs, t1_s_obs = t1.fit_vfa_2_point(s,fa_rad,tr_s,mask=mask)
# print(t1_s_obs)

#s0_map_obs, t1_map_obs = t1.vfit_vfa_2_point(s1,s2,fa1_rad,fa2_rad,tr_s)


#s0_map_obs_lin4d, t1_map_obs_lin4d = t1.fit_vfa_nonlinear(s_2d,fa_rad,tr_s)

s0_map_obs_lin4d, t1_map_obs_lin4d = t1.fit_vfa_linear(s_4d,fa_rad,tr_s)
print(t1_map_obs_lin4d)