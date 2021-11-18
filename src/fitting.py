"""Fitting.
Created on Thu Oct 21 15:50:47 2021
@authors: Michael Thrippleton
@email: m.j.thrippleton@ed.ac.uk
@institution: University of Edinburgh, UK

Classes: 
Functions: 
"""

from abc import ABC, abstractmethod
import nibabel as nib
import numpy as np
from utils.imaging import read_images


class calculator(ABC):
    # interface for classes that process data, e.g. fit DCE curve, estimate T1
    @abstractmethod
    def proc(self):
        # method to process a single data series (e.g. fit time series for one voxel)
        # should be overridden in subclasses
        # returns dict(values), series
        pass

    @abstractmethod
    def output_info(self):
        pass

    def proc_image(self, input_images, arg_images=None, mask=None, threshold=-np.inf, write_output=True, prefix="",
                   suffix="", filters=None, template=None):
        # method to process every voxel in an image (e.g every DCE time series in 4D image)
        # works by essentially looping over subclass proc() method

        # read image(s) containing input data to fit, then reshape to 2D (1 series of values per voxel)
        data, input_header = read_images(input_images)
        data_2d = data.reshape(-1, data.shape[-1])  # N voxels x N datapoints per voxel
        n_voxels, n_points = data_2d.shape

        # read argument images and reshape to 1D
        if arg_images is not None:
            # args = tuple; each element contains all voxels for an argument
            args, _hdrs = zip(
                *[read_images(a) if type(a) is not float else (np.tile(a, data.shape[:-1]), None) for a in arg_images])
            # args_1d = tuple; each element is a tuple of all arguments for a voxel
            args_1d = tuple(zip(*[a.reshape(-1) for a in args]))
        else:
            args_1d = [()] * n_voxels

        # read mask and reshape to 1D
        if mask is not None:
            mask_1d = nib.load(mask).get_fdata().reshape(-1) > 0
        else:
            mask_1d = np.empty(n_voxels, dtype=bool)
            mask_1d[:] = True

        # Prepare dict of pre-allocated output arrays
        outputs = {}
        for name, is1d in self.output_info().items():
            n_values = n_points if is1d else 1
            outputs[name] = np.empty((n_voxels, n_values), dtype=np.float32)
            outputs[name][:] = np.nan

        # process each voxel
        for i, voxel_data in enumerate(data_2d):
            if max(voxel_data) >= threshold and mask_1d[i]:
                try:
                    voxel_output = self.proc(voxel_data, *args_1d[i])
                    for name, values in voxel_output.items():
                        outputs[name][i, :] = values
                except (ValueError, ArithmeticError):
                    pass

        # filter outputs
        if filters is not None:
            for name, limits in filters.items():
                outputs[name][~(limits[0] <= outputs[name] <= limits[1])] = np.nan

        # reshape arrays to match image
        for name, values in outputs.items():
            outputs[name] = np.squeeze(outputs[name].reshape((*data.shape[:-1], values.shape[-1])))
        del data, data_2d

        # write output images
        if write_output:
            if template is not None:
                hdr = nib.load(template).header
            elif input_header is not None:
                hdr = input_header
            else:
                raise ValueError("Need input nifti files or template nifti file to write output images.")
            hdr.set_data_dtype(np.float32)
            for name, values in outputs.items():
                # MODIFY ACCORDING TO DIMENSIONALITY
                img = nib.nifti1.Nifti1Image(outputs[name], None, header=hdr)
                filename = f"{prefix}{name}{suffix}.nii"
                nib.save(img, filename)

        return outputs
