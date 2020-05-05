"""
This file implements the class for Euler's equation in Special Relativity.

To solve the conservative -> primitive transformation we use the scipy library.
"""

import numpy
from scipy.optimize import root_scalar


def prims_given_p(p, D, S, tau, gamma):
    """
    A helper function. Given the conserved variables and a value for the
    pressure, return the primitive variables

    Parameters
    ----------
    p : float
        Pressure
    D : float
        Density
    S : float
        Momentum
    tau : float
        Reduced energy
    gamma : float
        EOS information. (In general we should pass the EOS functions through)


    Returns
    -------
    rho, v, epsilon, W, h : float
        Primitive and auxilliary variables
    """
    v2 = S**2 / (tau + p + D)**2
    W = 1.0 / numpy.sqrt(1.0 - v2)
    rho = D / W
    h = (tau + p + D) / (rho * W**2)
    epsilon = h - 1.0 - p / rho
    v = S / (rho * h * W**2)
    
    return rho, v, epsilon, W, h
    

def f_for_c2p(pbar, D, S, tau, gamma):
    """
    Nonlinear function to find pressure consistent with conserved variables.

    Parameters
    ----------
    pbar : float
        Pressure guess
    D : float
        Density
    S : float
        Momentum
    tau : float
        Reduced energy
    gamma : float
        EOS information. (In general we should pass the EOS functions through)

    Returns
    -------
    f : float
        Residual f(pbar)
    """
    rhobar, _, epsilonbar, _, _ = prims_given_p(pbar, D, S, tau, gamma)
    
    return pbar - (gamma - 1.0) * rhobar * epsilonbar

class EulerSR(object):
    
    def __init__(self, gamma):
        self.gamma = gamma
        
    def p_from_eos(self, rho, epsilon):
        return (self.gamma - 1.0) * rho * epsilon
    
    def p2c(self, rho, v, epsilon):
        q = numpy.zeros([len(rho), 3])
        W = 1 / numpy.sqrt(1 - v**2)
        p = self.p_from_eos(rho, epsilon)
        h = 1 + epsilon + p / rho
        q[:, 0] = rho * W
        q[:, 1] = rho * h * W**2 * v
        q[:, 2] = rho * h * W**2 - p - rho * W
        return q
    
    def c2p(self, q):
        rho = numpy.zeros(len(q))
        v = numpy.zeros_like(rho)
        epsilon = numpy.zeros_like(rho)
        p = numpy.zeros_like(rho)
        W = numpy.zeros_like(rho)
        h = numpy.zeros_like(rho)
        cs = numpy.zeros_like(rho)
        for i in range(len(q)):
            D, S, tau = q[i, :]
            pmin = max(1e-10, abs(S) - tau - D)
            root_result = root_scalar(f_for_c2p, args=(D, S, tau, self.gamma), method='brentq', bracket=[pmin, 1e10])
            pbar = root_result.root
            rho[i], v[i], epsilon[i], W[i], h[i] = prims_given_p(pbar, D, S, tau, self.gamma)
            p[i] = pbar
            cs[i] = numpy.sqrt(self.gamma * p[i] / (rho[i] * h[i]))
        return rho, v, epsilon, p, W, h, cs
        
    def flux(self, q):
        D = q[:, 0]
        S = q[:, 1]
        tau = q[:, 2]
        _, v, _, p, _, _, _ = self.c2p(q)
        f = numpy.zeros_like(q)
        f[:, 0] = D * v
        f[:, 1] = S * v + p
        f[:, 2] = (tau + p) * v
        return f
    
    def max_lambda(self, q):
        _, v, _, _, _, _, cs = self.c2p(q)
        lambda_plus = (v + cs) / (1.0 + v * cs)
        lambda_minus = (v + cs) / (1.0 + v * cs)
        return max(numpy.max(numpy.abs(lambda_plus)),
                   numpy.max(numpy.abs(lambda_minus)))
    
