"""
Computing errors against an exact solution
"""

import numpy


def error_norm(q_numerical, q_exact, dx, p=2):
    """
    Compute the discrete error in q in the p norm

    Parameters
    ----------
    q_numerical : numpy vector
        The numerical solution, an array size (N,) or (N,1)
    q_exact : numpy vector
        The exact solution, whose size matches q_numerical
    dx : float
        The relevant grid spacing
    p : int or 'inf', optional
        The norm. The default is 2.

    Returns
    -------
    error_value : float
        (dx * sum((q_n - q_e)**p))**(1/p)
    """
    
    if p == 'inf':
        error_value = numpy.max(numpy.abs(q_numerical - q_exact))
    else:
        error_value = (dx * numpy.sum(numpy.abs(q_numerical - q_exact)**p))**(1/p)
    
    return error_value
