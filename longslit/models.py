# Third-party
import numpy as np

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

