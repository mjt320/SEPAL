"""Functions to fit MRI SPGR signal to obtain T2* or T2.

Created 17 October 2023
@authors: Michael Thrippleton
@email: m.j.thrippleton@ed.ac.uk
@institution: University of Edinburgh, UK

Classes:
    MultiEchoT2sLinear
    MultiEchoT2sNonLinear

Functions:
    multiecho_signal: get signal
"""

import numpy as np
from scipy.optimize import least_squares
from sepal.fitting import Fitter


class MultiEchoT2sLinear(Fitter):
    """Linear multi-echo T2* estimation.

    Subclass of Fitter.
    """

    def __init__(self, te, min_signal=0.):
        """

        Args:
            te (ndarray): 1D array containing the echo times (s)
            min_signal (float): exclude echoes where signal is below this
                value. Defaults to zero.

        """
        self.te = te
        self.min_signal = min_signal

    def output_info(self):
        """Get output info. Overrides superclass method.
        """
        return ('s0', False), ('t2s', False)

    def proc(self, s):
        """Estimate T2*. Overrides superclass method.

        Args:
            s (ndarray): 1D array containing the signals

        Returns:
            tuple: (s0, t1)
                s0 (float): signal without T2* decay
                t2s (float): T2* (s)

        """
        if any(np.isnan(s)):
            raise ArithmeticError(
                f'Unable to calculate T2s: nan signal received.')

        te_mask = s > self.min_signal
        if np.sum(te_mask) < 2:
            raise ArithmeticError(
                f'Unable to calculate T2s: fewer than 2 included TE values.')

        y = np.log(s[te_mask])
        x = self.te[te_mask]
        A = np.stack([x, np.ones(x.shape)], axis=1)
        slope, intercept = np.linalg.lstsq(A, y, rcond=None)[0]

        if (intercept < 0) or (0. < slope):
            raise ArithmeticError('T2s estimation failed.')

        t2s, s0 = -1. / slope, np.exp(intercept)

        return s0, t2s


class MultiEchoT2sNonLinear(Fitter):
    """Non-linear multi-echo T2* estimation.

    Subclass of Fitter.
    """

    def __init__(self, te, min_signal=0.):
        """

        Args:
            te (ndarray): 1D array containing the echo times (s)
            min_signal (float): exclude echoes where signal is below this
                value. Defaults to zero.

        """
        self.te = te
        self.min_signal = min_signal
        self.linear_fitter = MultiEchoT2sLinear(te, min_signal)

    def output_info(self):
        """Get output info. Overrides superclass method.
        """
        return ('s0', False), ('t2s', False)

    def proc(self, s):
        """Estimate T2*. Overrides superclass method.

        Args:
            s (ndarray): 1D array containing the signals

        Returns:
            tuple: (s0, t1)
                s0 (float): signal without T2* decay
                t2s (float): T2* (s)

        """
        if any(np.isnan(s)):
            raise ArithmeticError(
                f'Unable to calculate T2s: nan signal received.')

        te_mask = s > self.min_signal
        if np.sum(te_mask) < 2:
            raise ArithmeticError(
                f'Unable to calculate T2s: fewer than 2 included TE values.')

        # use linear fit to obtain initial guess, otherwise set arbitrary values
        try:
            x0 = np.array(self.linear_fitter.proc(s))
        except ArithmeticError:
            x0 = np.array([s[0], 10e-3])

        s = s[te_mask]
        te = self.te[te_mask]

        result = least_squares(self.__residuals, x0, args=(s, te), bounds=(
                (1e-8, 1e-8), (np.inf, np.inf)), method='trf', x_scale=x0)
        if result.success is False:
            raise ArithmeticError(f'Unable to fit multiecho data:'
                                  f' {result.message}')

        s0, t2s = result.x
        return s0, t2s

    @staticmethod
    def __residuals(x, s, te):
        s0, t2s = x
        s_est = multiecho_signal(s0, t2s, te)
        return s - s_est


def multiecho_signal(s0, t2s, te):
    """Return signal for multiple echoes / echo times.

    Parameters
    ----------
        s0 : float
             Signal without T2/T2* decay.
        t2s : float
             T2*/T2 value (s).
        te : float
             TE value (s).

    Returns
    -------
        s : float
            Signal.
    """
    s = s0 * np.exp(-te/t2s)

    return s
