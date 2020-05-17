"""
This file implements the class for Euler's equation in General Relativity,
for evolving the TOV equations.

To solve the conservative -> primitive transformation we use the scipy library.

Note: this is notably different to previous classes as it stores the
coordinates and the metric (which do not evolve - Cowling approximation).
That means a new instance is required to change grid resolution
"""

import numpy
from scipy.optimize import root_scalar


def prims_given_p(p, D, S, tau, a, gamma):
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
    v2 = S**2 / (tau + p + D)**2 / a**2
    W = 1.0 / numpy.sqrt(1.0 - v2)
    rho = D / W
    h = (tau + p + D) / (rho * W**2)
    epsilon = h - 1.0 - p / rho
    v = S / (rho * h * W**2) / a**2
    
    return rho, v, epsilon, W, h
    

def f_for_c2p(pbar, D, S, tau, a, gamma):
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
    rhobar, _, epsilonbar, _, _ = prims_given_p(pbar, D, S, tau, a, gamma)
    
    return pbar - (gamma - 1.0) * rhobar * epsilonbar

class EulerTOV(object):
    
    def __init__(self, gamma, r, initial_data_tov, rho_atmosphere=1e-10,
                 epsilon_atmosphere=1e-10):
        self.gamma = gamma
        self.r = r
        self.rho_atmosphere = rho_atmosphere
        self.epsilon_atmosphere = epsilon_atmosphere
        q0, alpha, a = initial_data_tov(r)
        self.q0 = q0
        self.alpha = alpha
        self.a = a
        self.m = self.r * (1.0 - 1.0 / self.a**2) / 2
        
    def p_from_eos(self, rho, epsilon):
        return (self.gamma - 1.0) * rho * epsilon
    
    def p2c(self, rho, v, epsilon):
        q = numpy.zeros([len(rho), 3])
        # Impose atmosphere based on rho
        v[rho < self.rho_atmosphere] = 0
        epsilon[rho < self.rho_atmosphere] = self.epsilon_atmosphere
        rho[rho < self.rho_atmosphere] = self.rho_atmosphere
        # Impose atmosphere based on epsilon
        v[epsilon < self.epsilon_atmosphere] = 0
        rho[epsilon < self.epsilon_atmosphere] = self.rho_atmosphere
        epsilon[epsilon < self.epsilon_atmosphere] = self.epsilon_atmosphere
        W = 1 / numpy.sqrt(1 - self.a**2 * v**2)
        p = self.p_from_eos(rho, epsilon)
        h = 1 + epsilon + p / rho
        q[:, 0] = rho * W
        q[:, 1] = rho * h * W**2 * v * self.a**2
        q[:, 2] = rho * h * W**2 - p - rho * W
        q *= self.a[:, numpy.newaxis] * self.r[:, numpy.newaxis]**2
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
            D, S, tau = q[i, :] / self.a[i] / self.r[i]**2
            # Impose atmosphere on D, tau
            if D < self.rho_atmosphere or tau < self.rho_atmosphere:
                D = self.rho_atmosphere
                S = 0
                tau = self.rho_atmosphere
            pmin = 1e-30
            root_result = root_scalar(f_for_c2p, args=(D, S, tau, self.a[i], self.gamma), method='brentq', bracket=[pmin, 1e10])
            pbar = root_result.root
            rho[i], v[i], epsilon[i], W[i], h[i] = prims_given_p(pbar, D, S, tau, self.a[i], self.gamma)
            # Impose atmosphere again
            if rho[i] < self.rho_atmosphere:
                v[i] = 0.0
                rho[i] = self.rho_atmosphere
                epsilon[i] = self.epsilon_atmosphere
                W[i] = 1.0
                pbar = self.rho_atmosphere * self.epsilon_atmosphere
                h[i] = 1.0 + 2 * epsilon[i]
            p[i] = pbar
            cs[i] = numpy.sqrt(self.gamma * p[i] / (rho[i] * h[i]))
        return rho, v, epsilon, p, W, h, cs
        
    def flux(self, q):
        D = q[:, 0] / self.a / self.r**2
        S = q[:, 1] / self.a / self.r**2
        tau = q[:, 2] / self.a / self.r**2
        _, v, _, p, _, _, _ = self.c2p(q)
        f = numpy.zeros_like(q)
        f[:, 0] = D * v
        f[:, 1] = S * v + p
        f[:, 2] = (tau + p) * v
        f *= self.alpha[:, numpy.newaxis] * self.a[:, numpy.newaxis] * self.r[:, numpy.newaxis]**2
        return f
    
    def source(self, q):
        s = numpy.zeros_like(q)
        D = q[:, 0] / self.a / self.r**2
        S = q[:, 1] / self.a / self.r**2
        tau = q[:, 2] / self.a / self.r**2
        _, v, _, p, _, _, _ = self.c2p(q)
        s[:, 1] = -self.a**2 * self.m / self.r**2 * (S * v + tau + p + D)
        s[:, 2] = -self.m / self.r**2 * S
        s *= self.alpha[:, numpy.newaxis] * self.a[:, numpy.newaxis]
        return s
    
    def max_lambda(self, q):
        """
        This could in principle be worked out, but it is easiest to use c=1.
        """
        return 1.0    
