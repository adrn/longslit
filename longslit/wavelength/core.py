# Third-party
import numpy as np
from scipy.optimize import minimize

# Project
from ..utils import gaussian_constant

__all__ = ['fit_emission_line']

def _errfunc(p, pix, flux, flux_ivar):
    if p[0] < 0 or p[2] < 0:
        return np.inf

    return np.sum((gaussian_constant(pix, *p) - flux)**2)

def fit_emission_line(pix_grid, flux, flux_ivar=None,
                      centroid0=None, sigma0=None, amp0=None, offset0=None):
    """
    TODO:

    Parameters
    ----------
    pix_grid : array_like
        Must be the same shape as ``flux``.
    flux : array_like
        Must be the same shape as ``pix_grid``.
    centroid0 : numeric (optional)
        Initial guess for line centroid.
    amp0 : numeric (optional)
        Initial guess for line amplitude.
    offset0 : numeric (optional)
        Initial guess for line offset above continuum.
    """

    if centroid0 is None: # then estimate the initial guess for the centroid
        centroid0 = np.argmax(flux)

    int_ctrd0 = int(round(centroid0))
    if amp0 is None: # then estimate the initial guess for amplitude
        amp0 = flux[int_ctrd0] # flux at initial guess

    if offset0 is None: # then estimate the initial guess for offset
        # TODO / MAGIC NUMBER: buffer hard-coded to 16??
        offset0 = flux[int_ctrd0-16:int_ctrd0+16].min()

    if sigma0 is None:
        sigma0 = 4. # MAGIC NUMBER

    if flux_ivar is None:
        flux_ivar = 1.

    p0 = (amp0, centroid0, sigma0, offset0)
    res = minimize(_errfunc, x0=p0, args=(pix_grid, flux, flux_ivar))
    p = res.x

    fail_msg = "Fitting spectral line in comp lamp spectrum failed. {msg}"
    if p[1] < min(pix_grid) or p[1] > max(pix_grid):
        raise ValueError(fail_msg.format(msg="Unphysical peak centroid: {:.3f}".format(p[0])))

    elif p[2] < 0.1 or p[2] > 10.:
        raise ValueError(fail_msg.format(msg="Unphysical line width: {:.3f}".format(p[2])))

    return dict(amp=p[0], centroid=p[1], stddev=p[2], const=p[3])
