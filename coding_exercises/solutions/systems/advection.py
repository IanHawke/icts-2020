"""
A system class gives the key functions for implementing conservation laws.

The conservation law is
$$
  \partial_t q + \partial_x f = 0.
$$
We need a function to return the flux f, and a function to return the maximum
characteristic speed (in order to find the safe timestep).

This file implements the class for the advection equation.
"""

class Advection(object):
    
    def __init__(self, v=1):
        self.v = v
        
    def flux(self, q):
        return self.v * q
    
    def max_lambda(self, q):
        return self.v