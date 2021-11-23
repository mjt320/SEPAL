"""Fitting.
Created on Thu Oct 21 15:50:47 2021
@authors: Michael Thrippleton
@email: m.j.thrippleton@ed.ac.uk
@institution: University of Edinburgh, UK

Classes:
Functions:
"""

from abc import ABC, abstractmethod
import os

import nibabel as nib
import numpy as np
from joblib import Parallel, delayed

from utils.imaging import read_images


class calculator(ABC):
    """Abstract base class for fitting algorithms.

    Subclasses must implement the proc method, which process a single data
    series, e.g. DCE time series for one voxel/ROI, set of values with
    different flip angles for T1 measurement. The proc_image method is
    provided for processing entire images by calling the proc method on each
    voxel.
    """

    @abstractmethod
    def proc(self, *args):
        """Abstract method processing a single data series.

        For example, estimating pharmacokinetic parameters from a DCE
        concentration-time series, or estimating T1 from a series of
        signals acquired at different flip angles. Overridden by subclass.

        Args:
            *args: First argument is the data, followed by any other arguments.

        Returns:
            dict: dictionary containing output parameters, which should be
                scalar (e.g. KTrans) or 1D ndarrays (e.g. fitted concentration
                series).
        """
        pass

    @abstractmethod
    def output_info(self):
        """Abstract method returning output names (dict keys) and dimension.

        Returns:
             dict: key=parameter name, value=True if 1D, False if scalar
        """
        pass

    def proc_image(self, input_images, arg_images=None, mask=None,
                   threshold=-np.inf, write_output=False, out_dir=".",
                   prefix="",
                   suffix="", filters=None, template=None, n_procs=1):
        """Process image using subclass proc method.

        Args:
            input_images (list): One or more input images to be processed. List
                can contain nifti filenames (str) or arrays (ndarray). If there
                is one image in list, the last dimension is assumed to be the
                series dimension (e.g. time, flip angle). If list contains >1
                images, they are concatenated along a new series dimension.
            arg_images (tuple): Tuple containing one image (str or ndarray,
                as above) for each argument required by the proc method -
                refer to subclass proc docstring for more info. Defaults to None
            mask (str or ndarray): Mask image (str or ndarray, as above). Must
                contain 1 or 0 only. 1 indicates voxels to be processed.
                Defaults to None (process all voxels).
            threshold (float): Voxel is processed if max input value in
                series (e.g. flip angle or time series) is >= threshold.
                Defaults to -np.inf
            write_output (bool): Outputs written to nifti files if True.
                Defaults to False.
            out_dir (str): Directory for output images.
            prefix (str): filename prefix for output images. Defaults to "".
            suffix (str): filename suffix for output images. Defaults to "".
            filters(dict): Dict of 2-tuples: key=parameter name, value=(lower
                limit, upper limit). Output values outside the range are set
                to nan.
            template (str): Nifti filename. Uses the header of this image to
                write output images. Defaults to None, in which case the header
                of the first input image will be used, otherwise an exception is
                raised.
            n_procs (int): Number of processes for parallel computation.
            n_chunks:

        Returns:
            dict: key=output parameter name, value=ndarray of values
        """

        # move to end if this works
        out_dir_path = os.path.join(os.getcwd(), out_dir)
        if write_output and not os.path.isdir(out_dir_path):
            os.mkdir(out_dir_path)

        # read source images, e.g. signal-time images
        data, input_header = read_images(input_images)
        # reshape data to 2D array n_voxels x n_points (length of series)
        data_2d = data.reshape(-1, data.shape[-1])
        n_voxels, n_points = data_2d.shape
        data_shape = data.shape

        # read argument images, e.g. flip angle correction, T10
        if arg_images is None:
            args_1d = [()] * n_voxels
        else:
            # get list of N-D arrays for each argument
            arg_arrays, _hdrs = zip(*[
                read_images(a) if type(a) is not float else
                (np.tile(a, data_shape[:-1]), None) for a in arg_images])
            # convert to list of (list of arguments) for each voxel
            args_1d = list(zip(*[a.reshape(-1) for a in arg_arrays]))
        del arg_images

        # read mask image
        if mask is None:
            mask_1d = np.empty(n_voxels, dtype=bool)
            mask_1d[:] = True
        else:
            mask_1d = nib.load(mask).get_fdata().reshape(-1) > 0
            if any((mask_1d != 0) & (mask_1d != 1)):
                raise ValueError('Mask contains elements that are not 0 or 1.')
            mask_1d = mask_1d.astype(bool)

        # divide data into 1+ "chunks" of voxels for parallel processing
        n_chunks = min(5 * n_procs, n_voxels)
        chunks_start_idx = np.int32(
            n_voxels * (np.array(range(n_chunks)) / n_chunks))

        # function to process a chunk of voxels, to be called by joblib
        def _proc_chunk(i_chunk):
            # work out voxel indices corresponding to the chunk
            start_voxel = chunks_start_idx[i_chunk]
            stop_voxel = chunks_start_idx[i_chunk + 1] if (
                    i_chunk != n_chunks - 1) else n_voxels
            n_chunk_voxels = stop_voxel - start_voxel
            # preallocate output arrays
            chunk_output = {}
            for name_, is1d in self.output_info().items():
                n_values = n_points if is1d else 1
                chunk_output[name_] = np.empty((n_chunk_voxels, n_values),
                                               dtype=np.float32)
                chunk_output[name_][:] = np.nan
            # process all voxels in the chunk
            for i_vox_chunk, i_vox in enumerate(np.arange(start_voxel,
                                                          stop_voxel, 1)):
                voxel_data = data_2d[i_vox, :]
                if max(voxel_data) >= threshold and mask_1d[i_vox]:
                    try:
                        voxel_output = self.proc(voxel_data, *args_1d[i_vox])
                        for name_, values_ in voxel_output.items():
                            chunk_output[name_][i_vox_chunk, :] = values_
                    except (ValueError, ArithmeticError):
                        pass
            return chunk_output

        # run the processing
        chunk_outputs = Parallel(n_jobs=n_procs)(
            delayed(_proc_chunk)(i_chunk) for i_chunk in range(n_chunks))
        del data, data_2d, args_1d, mask_1d

        # Combine chunks into single output dict
        outputs = {key: np.concatenate(
            [co[key] for co in chunk_outputs], axis=0
        )
            for key in self.output_info().keys()}
        del chunk_outputs

        # filter outputs
        if filters is not None:
            for name, limits in filters.items():
                outputs[name][
                    ~(limits[0] <= outputs[name] <= limits[1])] = np.nan

        # reshape arrays to match image dimensions
        for name, values in outputs.items():
            outputs[name] = np.squeeze(outputs[name].reshape(
                (*data_shape[:-1], values.shape[-1])))

        # write output images if required
        if write_output:
            if template is not None:
                hdr = nib.load(template).header
            elif input_header is not None:
                hdr = input_header
            else:
                raise ValueError("Need input nifti files or template nifti "
                                 "file to write output images.")
            hdr.set_data_dtype(np.float32)
            for name, values in outputs.items():
                img = nib.nifti1.Nifti1Image(outputs[name], None, header=hdr)
                filename = f"{prefix}{name}{suffix}.nii"
                nib.save(img, filename)

        return outputs
