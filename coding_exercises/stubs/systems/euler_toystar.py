"""
This file implements the class for Euler's equation.

This includes a source term to model a toy star, and an atmosphere.
"""

import numpy


class EulerToyStar(object):
    
    def __init__(self, gamma=2, rho_atmosphere=1e-6, e_atmosphere=1e-6):
        assert gamma == 2
        self.gamma = gamma
        self.rho_atmosphere = rho_atmosphere
        self.e_atmosphere = e_atmosphere
        
    def p_from_eos(self, rho, e):
        return (self.gamma - 1.0) * rho * e
    
    def p2c(self, rho, v, e):
        q = numpy.zeros([len(rho), 3])
        # Impose atmosphere based on rho
        v[rho < self.rho_atmosphere] = 0
        e[rho < self.rho_atmosphere] = #! To be completed
        rho[rho < self.rho_atmosphere] = #! To be completed
        # Impose atmosphere based on e
        #! To be completed
        q[:, 0] = #! To be completed
        q[:, 1] = #! To be completed
        q[:, 2] = #! To be completed
        return q
    
    def c2p(self, q):
        rho = q[:, 0]
        S = q[:, 1]
        E = q[:, 2]
        # Impose atmosphere based on rho
        #! To be completed
        v = #! To be completed
        e = #! To be completed
        # Impose atmosphere based on e
        #! To be completed
        p = self.p_from_eos(rho, e)
        cs = #! To be completed
        return rho, v, e, p, cs
        
    def flux(self, q):
        # rho = q[:, 0, :]
        S = q[:, 1]
        E = q[:, 2]
        _, v, _, p, _ = self.c2p(q)
        f = numpy.zeros_like(q)
        f[:, 0] = #! To be completed
        f[:, 1] = #! To be completed
        f[:, 2] = #! To be completed
        return f
    
    def source(self, q, x):
        rho, v, _, _, _ = self.c2p(q)
        s = numpy.zeros_like(q)
        s[:, 1] = #! To be completed
        s[:, 2] = #! To be completed
        return s
    
    def max_lambda(self, q):
        _, v, _, _, cs = self.c2p(q)
        return #! To be completed
    
