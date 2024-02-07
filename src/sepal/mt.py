"""Functions to fit MRI SPGR signal to obtain MT parameters.

Created 13 November 2023
@authors: Michael Thrippleton
@email: m.j.thrippleton@ed.ac.uk
@institution: University of Edinburgh, UK

Classes:
    MTR
    MTSat

"""

import numpy as np
from sepal.fitting import Fitter


class MTR(Fitter):
    """Magnetisation transfer ratio calculation.

    Calculates MTR based on MT ON and MT OFF acquisitions.
    Subclass of Fitter.
    """

    def __init__(self):
        """

        Args:
        """
        pass

    def output_info(self):
        """Get output info. Overrides superclass method.
        """
        return ('mtr', False),

    def proc(self, s):
        """Calculate MTR. Overrides superclass method.

        Args:
            s (ndarray): 1D np array containing the MT ON, MT OFF signals

        Returns:
            tuple: (mtr)
                mtr (float): magnetisation transfer ratio (percent)
        """
        if any(np.isnan(s)):
            raise ArithmeticError(
                f'Unable to calculate MTR: nan signal received.')

        if s.shape != (2,):
            raise ValueError(
                f'Argument s should have shape (2,)'
            )

        return 100*(s[1]-s[0])/s[1]


class MTSat(Fitter):
    """Magnetisation transfer saturation calculation.

    Calculates MTSat, MTR and T1 based on three acquisitions with MT, PD and
    T1 weightings.
    Subclass of Fitter.
    """

    def __init__(self, fa_mt, fa_pd, fa_t1, tr_pd, tr_t1):
        """

        Args:
            fa_mt (float): excitation flip angle for MTw acquisition (deg)
            fa_pd (float): excitation flip angle for PDw acquisition (deg)
            fa_t1 (float): excitation flip angle for T1w acquisition (deg)
            tr_pd (float): TR for PD and MT weighted acquisitions (s)
            tr_t1 (float): TR for T1 weighted acquisition (s)
        """
        self.fa_mt = fa_mt * np.pi / 180
        self.fa_pd = fa_pd * np.pi / 180
        self.fa_t1 = fa_t1 * np.pi / 180
        self.tr_pd = tr_pd
        self.tr_t1 = tr_t1

    def output_info(self):
        """Get output info. Overrides superclass method.
        """
        return ('mtsat', False), ('mtr', False), ('t1', False), ('r1', False), \
               ('a', False)

    def proc(self, s):
        """Calculate MTSat, MTR and T1. Overrides superclass method.
        Uses equations from:
            Helms et al., Magnetic Resonance in Medicine 60:1396–1407 (2008)
            https://doi.org/10.1002/mrm.21732
            and erratum https://doi.org/10.1002/mrm.22607

        Args:
            s (ndarray): 1D np array containing the MT, PD and T1
            weighted signals in that sequence

        Returns:
            tuple: (mtsat, mtr and t1)
                mtsat (float): magnetisation transfer saturation
                mtr (float): magnetisation transfer ratio (percent)
                t1 (float): T1 estimation (s)
                r1 (float): R1 estimation (s^-1)
                a (float): amplitude estimation
        """
        if any(np.isnan(s)):
            raise ArithmeticError(
                f'Unable to calculate MTSat: nan signal received.')

        if s.shape != (3,):
            raise ValueError(
                f'Argument s should have shape (3,)'
            )
        s_mt, s_pd, s_t1 = s
        mtr = 100*(s_pd-s_mt)/s_pd
        a = (s_t1*(s_pd*self.tr_t1*self.fa_pd**2 -
                   s_pd*self.tr_pd*self.fa_t1**2)) / (
            self.fa_pd*self.fa_t1*(s_pd*self.tr_t1*self.fa_pd -
                                 s_t1*self.tr_pd*self.fa_t1))
        r1 = -(self.fa_pd*self.fa_t1*(s_pd*self.tr_t1*self.fa_pd -
                                    s_t1*self.tr_pd*self.fa_t1)) / (
                2*self.tr_t1*self.tr_pd*(s_pd*self.fa_t1 - s_t1*self.fa_pd))
        mtsat = -self.fa_pd**2/2 + (a*r1*self.tr_pd*self.fa_pd)/s_mt - \
                r1*self.tr_pd
        t1 = 1/r1

        return mtsat, mtr, t1, r1, a
