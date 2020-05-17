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

import numpy


class Advection(object):
    
    def __init__(self, v=1):
        self.v = v
        
    def flux(self, q):
        return #! To be completed
    
    def max_lambda(self, q):
        return #! To be completed
    

class Advection2d(object):
    
    def __init__(self, v=(1, 1)):
        self.v = numpy.array(v)
        
    def flux(self, q, dirn):
        return #! To be completed
    
    def max_lambda(self, q):
        return #! To be completed
