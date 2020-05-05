"""
This file implements the class for Burgers equation.
"""

import numpy


class Burgers(object):
    
    def __init__(self):
        pass
        
    def flux(self, q):
        return q**2 / 2
    
    def max_lambda(self, q):
        return numpy.max(numpy.abs(q))
