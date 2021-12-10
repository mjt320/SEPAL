"""Functions to fit MRI SPGR signal to obtain T1.

Created 28 September 2020
@authors: Michael Thrippleton
@email: m.j.thrippleton@ed.ac.uk
@institution: University of Edinburgh, UK

Functions:
    fit_vfa_2_point: obtain T1 using analytical formula based on two images
    fit_vfa_linear: obtain T1 using linear regression
    fit_vfa_nonlinear: obtain T1 using non-linear least squares fit
    fit_hifi: obtain T1 by fitting a combination of SPGR and IR-SPGR scans
    spgr_signal: get SPGR signal
    irspgr_signal: get IR-SPGR signal
"""

import numpy as np
from scipy.optimize import curve_fit, least_squares
from fitting import Fitter


class vfa_2points(Fitter):
    def __init__(self, fa, tr):
        self.fa = np.asarray(fa)
        self.tr = tr
        self.fa_rad = np.pi*self.fa/180

    def proc(self, s, k_fa=1):
        with np.errstate(divide='ignore', invalid='ignore'):
            fa_true = k_fa * self.fa_rad
            sr = s[0] / s[1]
            t1 = self.tr / np.log(
                (sr*np.sin(fa_true[1])*np.cos(fa_true[0]) -
                 np.sin(fa_true[0])*np.cos(fa_true[1])) /
                (sr*np.sin(fa_true[1]) - np.sin(fa_true[0])))
            s0 = s[0] * ((1-np.exp(-self.tr/t1)*np.cos(fa_true[0])) /
                         ((1-np.exp(-self.tr/t1))*np.sin(fa_true[0])))

        t1 = np.nan if ~np.isreal(t1) | (t1 <= 0) | np.isinf(t1) else t1
        s0 = np.nan if (s0 <= 0) | np.isinf(s0) else s0

        return {'s0': s0, 't1': t1}


class vfa_linear(Fitter):
    def __init__(self, fa, tr):
        self.fa = np.asarray(fa)
        self.tr = tr
        self.fa_rad = np.pi*self.fa/180

    def proc(self, s, k_fa=1):
        fa_true = k_fa * self.fa_rad
        y = s / np.sin(fa_true)
        x = s / np.tan(fa_true)
        A = np.stack([x, np.ones(x.shape)], axis=1)
        slope, intercept = np.linalg.lstsq(A, y, rcond=None)[0]

        is_valid = (intercept > 0) and (0. < slope < 1.)
        t1, s0 = (-self.tr/np.log(slope),
                  intercept/(1-slope)) if is_valid else (np.nan, np.nan)

        return {'s0': s0, 't1': t1}


class vfa_nonlinear(Fitter):
    def __init__(self, fa, tr):
        self.fa = np.asarray(fa)
        self.tr = tr
        self.fa_rad = np.pi*self.fa/180
        self.linear_fitter = vfa_linear(fa, tr)

    def proc(self, s, k_fa=1):
        # use linear fit to obtain initial guess
        result_linear = self.linear_fitter.proc(s, k_fa=k_fa)
        x_linear = np.array((result_linear['s0'], result_linear['t1']))
        if (~np.isnan(x_linear[0]) & ~np.isnan(x_linear[1])):
            x0 = x_linear
        else:
            x0 = np.array([s[0] / spgr_signal(1., 1., self.tr, k_fa*self.fa[0]), 1.])

        result = least_squares(self.__residuals, x0, args=(s, k_fa), bounds=((1e-8,1e-8),(np.inf,np.inf)), method='trf',
                           x_scale=x0
                           )
        if result.success is False:
            raise ArithmeticError(f'Unable to fit VFA data'
                                  f': {result.message}')

        s0, t1 = result.x
        return {'s0': s0, 't1': t1}

    def __residuals(self, x, s, k_fa):
        s0, t1 = x
        s_est = spgr_signal(s0, t1, self.tr, k_fa*self.fa)
        return s - s_est


class hifi(Fitter):
    def __init__(self, esp, ti, n, b, td, centre):
        self.esp = esp
        self.ti = ti
        self.n = n
        self.b = b
        self.td = td
        self.centre = centre
        # get information about the scans
        self.n_scans = len(esp)
        self.is_ir = ~np.isnan(ti)
        self.is_spgr = ~self.is_ir
        self.idx_spgr = np.where(self.is_spgr)[0]
        self.n_spgr = self.idx_spgr.size
        self.get_linear_estimate = self.n_spgr > 1 and np.all(
            np.isclose(esp[self.idx_spgr], esp[self.idx_spgr[0]]))
        self.linear_fitter = vfa_linear( b[self.is_spgr], esp[self.idx_spgr[0]])

    def proc(self, s, k_fa_fixed=None):
        # First get a quick linear T1 estimate
        if self.get_linear_estimate:  # If >1 SPGR, use linear VFA fit
            result_lin = self.linear_fitter.proc(s[self.is_spgr])
            if ~np.isnan(result_lin['s0']) and ~np.isnan(result_lin['t1']):
                s0_init, t1_init = result_lin['s0'], result_lin['t1']
            else:  # if result invalid, assume T1=1
                t1_init = 1
                s0_init = s[self.idx_spgr[0]] / spgr_signal(1, t1_init,
                                                   self.esp[self.idx_spgr[0]],
                                                   self.b[self.idx_spgr[0]])
        elif self.n_spgr == 1: # If 1 SPGR, assume T1=1 and estimate s0 based on this scan
            t1_init = 1
            s0_init = s[self.idx_spgr[0]] / spgr_signal(1, t1_init,
                                               self.esp[self.idx_spgr[0]],
                                               self.b[self.idx_spgr[0]])
        else: # If 0 SPGR, assume T1=1 and estimate s0 based on 1st scan
            t1_init = 1
            s0_init = s[0] / irspgr_signal(1, t1_init, self.esp[0], self.ti[0], self.n[0], self.b[0],
                                       180, self.td[0], self.centre[0])

        # Non-linear fit
        if k_fa_fixed is None:
            k_init = 1
            bounds = ([0, 0, 0], [np.inf, np.inf, np.inf])
        else:
            k_init = k_fa_fixed
            bounds = ([0, 0, 1], [np.inf, np.inf, 1])
        x_0 = np.array([t1_init, s0_init, k_init])
        result = least_squares(self.__residuals, x_0, args=(s,), bounds=bounds, method='trf',
                           x_scale=(t1_init, s0_init, k_init)
                           )
        x_opt = result.x if result.success else (np.nan, np.nan, np.nan)
        t1_opt, s0_opt, k_fa_opt = x_opt
        s_opt = self.__signal(x_opt)
        return {'t1': t1_opt, 's0': s0_opt, 'k_fa': k_fa_opt, 's_opt': s_opt}

    def __residuals(self, x, s):
        return s - self.__signal(x)

    def __signal(self, x):
        t1, s0, k_fa = x
        s = np.zeros(self.n_scans)
        s[self.is_ir] = irspgr_signal(s0, t1, self.esp[self.is_ir], self.ti[self.is_ir],
                                     self.n[self.is_ir], k_fa*self.b[self.is_ir], self.td[self.is_ir],
                                     self.centre[self.is_ir])
        s[self.is_spgr] = spgr_signal(s0, t1, self.esp[self.is_spgr],
                                     k_fa*self.b[self.is_spgr])
        return s


def spgr_signal(s0, t1, tr, fa):
    """Return signal for SPGR sequence.

    Parameters
    ----------
        s0 : float
             Equilibrium signal.
        t1 : float
             T1 value (s).
        tr : float
             TR value (s).
        fa : float
                 Flip angle (deg).

    Returns
    -------
        s : float
            Steady-state SPGR signal.
    """
    fa_rad = np.pi*fa/180

    e = np.exp(-tr/t1)
    s = s0 * (((1-e)*np.sin(fa_rad)) /
              (1-e*np.cos(fa_rad)))

    return s


def irspgr_signal(s0, t1, esp, ti, n, b, td=0, centre=0.5):
    """Return signal for IR-SPGR sequence.

    Uses formula by Deichmann et al. (2000) to account for modified
    apparent relaxation rate during the pulse train. Note inversion is assumed
    to be ideal.

    Parameters
    ----------
        s0 : float
             Equilibrium signal.
        t1 : float
             T1 value (s).
        esp : float
             Echo spacing (s). For SPGR, this is the TR.
        ti : float
             Inversion time (s). Note this is the actual time delay between the
             inversion pulse and the start of the echo train. The effective TI
             may be different, e.g for linear phase encoding of the echo train.
        n : int
            Number of excitation pulses per inversion pulse
        b : float
            Readout pulse flip angle (deg)
        td : float
             Delay between end of readout train and the next inversion (s).
        centre : float
                 Time in readout train when centre of k-space is acquired,
                 expressed as a fraction of the readout duration. e.g. = 0 for
                 centric phase encoding, = 0.5 for linear phase encoding.

    Returns
    -------
        s : float
            Steady-state IR-SPGR signal.
    """
    b_rad = np.pi*b/180
    tau = esp * n
    t1_star = (1/t1 - 1/esp*np.log(np.cos(b_rad)))**-1
    m0_star = s0 * ((1-np.exp(-esp/t1)) / (1-np.exp(-esp/t1_star)))

    r1 = -tau/t1_star
    e1 = np.exp(r1)
    e2 = np.exp(-td/t1)
    e3 = np.exp(-ti/t1)

    a1 = m0_star * (1-e1)
    a2 = s0 * (1 - e2)
    a3 = s0 * (1 - e3)

    a = a3 - a2*e3 - a1*e2*e3
    b = -e1*e2*e3

    m1 = a/(1-b)

    s = np.abs((
        m0_star + (m1-m0_star)*np.exp(centre*r1))*np.sin(b_rad))

    return s