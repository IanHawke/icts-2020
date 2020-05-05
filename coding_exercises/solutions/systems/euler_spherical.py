"""
This file implements the class for Euler's equation in spherical coordinates.
"""

import numpy


class EulerSpherical(object):
    
    def __init__(self, gamma=7/5):
        self.gamma = gamma
        
    def p_from_eos(self, rho, e):
        return (self.gamma - 1.0) * rho * e
    
    def p2c(self, rho, v, e, r):
        q = numpy.zeros([len(rho), 3])
        q[:, 0] = rho
        q[:, 1] = rho*v
        q[:, 2] = rho*(e + v**2/2)
        return r[:, numpy.newaxis]**2 * q
    
    def c2p(self, q, r):
        rho = q[:, 0] / r**2
        S = q[:, 1] / r**2
        E = q[:, 2] / r**2
        v = S / rho
        e = E / rho - v**2 / 2
        p = self.p_from_eos(rho, e)
        cs = numpy.sqrt(self.gamma * p / rho)
        return rho, v, e, p, cs
        
    def flux(self, q, r):
        # rho = q[:, 0, :]
        S = q[:, 1] / r**2
        E = q[:, 2] / r**2
        _, v, _, p, _ = self.c2p(q, r)
        f = numpy.zeros_like(q)
        f[:, 0] = S
        f[:, 1] = S * v + p
        f[:, 2] = (E + p) * v
        return r[:, numpy.newaxis]**2 * f
    
    def source(self, q, r):
        _, _, _, p, _ = self.c2p(q, r)
        s = numpy.zeros_like(q)
        s[:, 1] = 2 * p * r
        return s
    
    def max_lambda(self, q, r):
        """
        Note that the geometric factors may have an impact on the cell
        interpretation here!
        """
        _, v, _, _, cs = self.c2p(q, r)
        return numpy.max(numpy.abs(v)) + numpy.max(numpy.abs(cs))
    
