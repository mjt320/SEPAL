# -*- coding: utf-8 -*-
"""
Created on Thu Sep 24 22:18:12 2020.

@author: Michael Thrippleton

Test script to process VFA nii data using mrifit.t1 functions
"""

import os
import sys

import numpy as np
import nibabel as nib

sys.path.append("..")
from mrifit import t1

#acquisition and processing parameters
tr_s = 4.e-3
fa_rad = (np.pi/180.) * np.array((2.,5.,12.))
threshold = 50.

#load data from nii files
data_path = "./test_data/phantom_vfa"
file_names = ("series19.nii", "series20.nii", "series21.nii")
img_1 = nib.load(os.path.join(data_path, 'series19.nii'))
images = [nib.load(os.path.join(data_path,file_name)) for file_name in file_names]
signals = np.stack( [img.get_fdata() for img in images], axis=3 ) #4D signal array

#zero all slices except one to make things faster
signals[:,0:120,:,:]=0.
signals[:,122:,:,:]=0.

#create a mask based on threshold intensity
mask = (signals[...,0] > threshold).astype(int)

# s0, t1 = t1.fit_vfa_2_point(signals, fa_rad, tr_s, idx=[0,-1], mask=mask) #2-point method (fast)
s0, t1 = t1.fit_vfa_linear(signals, fa_rad, tr_s, mask=mask) #linear method (medium)
# s0, t1 = t1.fit_vfa_nonlinear(signals, fa_rad, tr_s, mask=mask) #non-linear method (slow)

#save output as nii images
new_hdr = images[0].header.copy()
new_hdr.set_data_dtype(np.float32)
t1_img = nib.nifti1.Nifti1Image(t1, None, header = new_hdr)
s0_img = nib.nifti1.Nifti1Image(s0, None, header = new_hdr)
nib.save(t1_img, "./results/t1_map.nii")
nib.save(s0_img, "./results/s0_map.nii")