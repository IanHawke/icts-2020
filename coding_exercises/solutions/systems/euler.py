"""
This file implements the class for Euler's equation.

Note that, for the Newtonian Euler equation, this could be simplified.
In particular, often the primitive variables are never explicit, as
everything is algebraic in terms of the conserved. This structure is used
as it carries over to more complex cases.
"""

import numpy


class Euler(object):
    
    def __init__(self, gamma):
        self.gamma = gamma
        
    def p_from_eos(self, rho, e):
        return (self.gamma - 1.0) * rho * e
    
    def p2c(self, rho, v, e):
        q = numpy.zeros([len(rho), 3])
        q[:, 0] = rho
        q[:, 1] = rho*v
        q[:, 2] = rho*(e + v**2/2)
        return q
    
    def c2p(self, q):
        rho = q[:, 0]
        S = q[:, 1]
        E = q[:, 2]
        v = S / rho
        e = E / rho - v**2 / 2
        p = self.p_from_eos(rho, e)
        cs = numpy.sqrt(self.gamma * p / rho)
        return rho, v, e, p, cs
        
    def flux(self, q):
        # rho = q[:, 0, :]
        S = q[:, 1]
        E = q[:, 2]
        _, v, _, p, _ = self.c2p(q)
        f = numpy.zeros_like(q)
        f[:, 0] = S
        f[:, 1] = S * v + p
        f[:, 2] = (E + p) * v
        return f
    
    def max_lambda(self, q):
        _, v, _, _, cs = self.c2p(q)
        return numpy.max(numpy.abs(v)) + numpy.max(numpy.abs(cs))
    
