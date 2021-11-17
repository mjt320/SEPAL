"""Utilities.

Created 6 October 2021
@authors: Michael Thrippleton
@email: m.j.thrippleton@ed.ac.uk
@institution: University of Edinburgh, UK

Functions:
    minimize_global
"""

import nibabel as nib
import numpy as np
from scipy.optimize import minimize, least_squares


def read_images(images):
    images = [images] if type(images) is not list else images
    if all([type(i) is str for i in images]):
            data = np.stack([nib.load(i).get_fdata() for i in images], axis=-1)
            header = nib.load(images[0]).header
    elif all([type(i) is np.ndarray for i in images]):
            data = np.stack(images, axis=-1)
            header = None
    else:
        raise TypeError('Argument images should contain all strings or all ndarrays.')
    if data.shape[-1] == 1:
        data = data.squeeze(axis=-1)
    return data, header

def roi_measure(image, mask_image):
    data, _hdr = read_images(image)
    mask, _hdr = read_images(mask_image)
    
    if mask.ndim == data.ndim:
        data = np.expand_dims(data, axis=0)
    
    # flatten spatial dimensions
    data_2d = data.reshape(-1, data.shape[-1])
    mask_1d = mask.reshape(-1)

    if not all((mask_1d==0) | (mask_1d==1)):
        raise ValueError('Mask contains values that are not 0 or 1.')

    masked_voxels = data_2d[mask_1d==1, :]
    mean = np.squeeze([np.nanmean(m_d) for m_d in masked_voxels.transpose()])
    median = np.squeeze([np.nanmedian(m_d) for m_d in masked_voxels.transpose()])
    sd = np.squeeze([np.nanstd(m_d) for m_d in masked_voxels.transpose()])
    
    return {'mean': mean, 'median': median, 'sd': sd}
    
