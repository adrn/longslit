# Standard library
import os
import sys

# Third-party
import astropy.units as u
import numpy as np
from scipy.optimize import curve_fit

# Project
from ..log import logger
from ..utils import gaussian_constant

__all__ = ['get_emission_line_centroid']

def get_emission_line_centroid(pix_grid, flux, centroid0=None, amp0=None, offset0=None):
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

    p0 = (amp0, centroid0, 1., offset0)
    p,_ = curve_fit(gaussian_constant, centroid0, flux, p0=p0)
    peak_pix = p[1]

    fail_msg = "Fitting spectral line in comp lamp spectrum failed. {msg}"
    if peak_pix < 0. or peak_pix > len(flux):
        raise ValueError(fail_msg.format(msg="Unphysical peak centroid: {:.3f}".format(p[0])))

    elif p[2] < 1. or p[2] > 10.:
        raise ValueError(fail_msg.format(msg="Unphysical line width: {:.3f}".format(p[2])))

    return peak_pix
