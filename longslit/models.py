# Third-party
import numpy as np
from scipy.special import wofz

sqrt_2pi = np.sqrt(2*np.pi)
def gaussian(x, amp, mu, sigma):
    return amp/(sqrt_2pi*sigma) * np.exp(-0.5 * ((np.array(x) - mu)/sigma)**2)

def gaussian_polynomial(x, amp, mu, sigma, *coeff):
    """
    Normalized Gaussian plus a polynomial.

    Parameters
    ----------
    x : numeric, array_like
    amp : numeric
    mu : numeric
    sigma : numeric
    *coeff :
        Any other arguments are interpreted as coefficients for a
        polynomial in x. Follows ordering in `numpy.polyval` - decreasing
        power!
    """
    return gaussian(x, amp, mu, sigma) + np.polyval(coeff, x)

def gaussian_constant(x, amp, mu, sigma, offset):
    """
    Normalized Gaussian plus a constant offset.

    Parameters
    ----------
    x : numeric, array_like
    amp : numeric
    mu : numeric
    sigma : numeric
    offset : numeric
    """
    return gaussian_polynomial(x, amp, mu, sigma, offset)

def voigt(x, amp, x0, G_std, L_fwhm):
    """
    Voigt profile - convolution of a Gaussian and Lorentzian.

    When G_std -> 0, the profile approaches a Lorentzian. When L_fwhm=0,
    the profile is a Gaussian.

    Parameters
    ----------
    x : numeric, array_like
    amp : numeric
        Amplitude of the profile (integral).
    x0 : numeric
        Centroid.
    G_std : numeric
        Standard of deviation of the Gaussian component.
    L_fwhm : numeric
        FWHM of the Lorentzian component.
    """
    _x = x-x0
    z = (_x + 1j*L_fwhm/2.) / (np.sqrt(2.)*G_std)
    return amp * wofz(z).real / (np.sqrt(2.*np.pi)*G_std)
