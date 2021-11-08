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
    
    # THIS WONT WORK. SORT OUT DIMENSIONS AND LOOPING
    if mask.ndim == data.ndim:
        data = np.expand_dims(data, axis=0)
    
    if not all((mask==0) | (mask==1)):
        raise ValueError('Mask contains values that are not 0 or 1.')
    
    # flatten spatial dimensions
    data = data.reshape(data.shape[0], -1)
    mask = mask.reshape(mask.shape[0], -1)
    masked_data = data[:, mask==1]

    mean = np.squeeze([np.nanmean(m_d) for m_d in masked_data])
    median = np.squeeze([np.nanmedian(m_d) for m_d in masked_data])
    sd = np.squeeze([np.nanstd(m_d) for m_d in masked_data])
    
    return {'mean': mean, 'median': median, 'sd': sd}
    
