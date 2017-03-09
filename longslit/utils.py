# Third-party
import numpy as np

sqrt_2pi = np.sqrt(2*np.pi)

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
    return amp/(sqrt_2pi*sigma) * np.exp(-0.5 * ((np.array(x) - mu)/sigma)**2) + offset
