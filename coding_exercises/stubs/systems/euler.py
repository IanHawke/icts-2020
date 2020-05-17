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
        #! To be completed
        return q
    
    def c2p(self, q):
        rho = #! To be completed
        S = #! To be completed
        E = #! To be completed
        v = #! To be completed
        e = #! To be completed
        p = #! To be completed
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
    
    def max_lambda(self, q):
        _, v, _, _, cs = self.c2p(q)
        return #! To be completed
    
