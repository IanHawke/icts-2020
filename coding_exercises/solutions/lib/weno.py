"""
Weighted Essentially Non-Oscillatory (WENO) reconstruction
"""

import numpy
from numba import jit


# Coefficients of order r=2
# On smooth solutions this should converge with order r=3
C_2 = numpy.array([ 1,  2 ]) / 3
a_2 = numpy.array([
                   [ 3, -1],
                   [ 1,  1],
                  ]) / 2
sigma_2 = numpy.array([
                        [
                          [ 1,  0],
                          [-2,  1]
                        ],
                        [
                          [ 1,  0],
                          [-2,  1]
                        ]
                      ])

# Coefficients of order r=3
# On smooth solutions this should converge with order r=5
C_3 = numpy.array([ 1,  6,  3 ]) / 10
a_3 = numpy.array([
                   [ 11,  -7,   2],
                   [  2,   5,  -1],
                   [ -1,   5,   2],
                  ]) / 6
sigma_3 = numpy.array([
                        [
                          [ 10,   0,   0],
                          [-31,  25,   0],
                          [ 11, -19,   4]
                        ],
                        [
                          [  4,   0,   0],
                          [-13,  13,   0],
                          [  5, -13,   4]
                        ],
                        [
                          [  4,   0,   0],
                          [-19,  25,   0],
                          [ 11, -31,  10]
                        ]
                      ]) / 3


@jit
def weno_kernel(q, order, C, a, sigma):
    """
    Do WENO reconstruction

    Parameters
    ----------

    q : numpy array
        Scalar data to reconstruct
    C, a, sigma : arrays of float
        Constants needed for the algorithm

    Returns
    -------

    q_{i+1/2} : float
        Reconstructed data
    """
    
    beta = numpy.zeros(order)
    w = numpy.zeros_like(beta)
    epsilon = 1e-16
    q_stencils = numpy.zeros_like(beta)
    alpha = numpy.zeros_like(beta)
    for k in range(order):
        for l in range(order):
            for m in range(l+1):
                beta[k] += sigma[k, l, m] * q[order-1+k-l] * q[order-1+k-m]
        alpha[k] = C[k] / (epsilon + beta[k]**2)
        for l in range(order):
            q_stencils[k] += a[k, l] * q[order-1+k-l]
    w = alpha / numpy.sum(alpha)
    q_weno = numpy.dot(w, q_stencils)

    return q_weno


def weno(q, order):
    """
    Do WENO reconstruction

    Parameters
    ----------

    q : numpy array
        Scalar data to reconstruct
    order : int
        Size of stencil order leads to optimal accuracy 2 order - 1

    Returns
    -------

    q_{i+1/2} : float
        Reconstructed data
    """
    if order == 2:
        q_weno = weno_kernel(q, order, C_2, a_2, sigma_2)
    elif order == 3:
        q_weno = weno_kernel(q, order, C_3, a_3, sigma_3)
    else:
        raise NotImplementedError

    return q_weno