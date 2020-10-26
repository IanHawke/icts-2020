"""
Constructing any grids that might be required
"""

import numpy


class Grid(object):
    """
    A finite volume style 1d grid
    """
    def __init__(self, interval, nx, ngz=1, cfl_factor=0.9):
        """
        Set up a finite-volume style 1d grid
    
        Parameters
        ----------
        interval : [float, float]
            The problem domain [x_L, x_R]
        nx : int
            The number of (interior) grid cells.
        ngz : int
            The number of ghost zones for the boundaries
        cfl_factor : float
            dt / dx
    
        Returns
        -------
        x : array of float
            Cell centre coordinates
        dx : float
            Cell widths
        """
        
        x_L, x_R = interval
        self.x_L = x_L
        self.x_R = x_R
        self.nx = nx
        self.ngz = ngz
        self.dx = (x_R - x_L) / nx
        self.dt = cfl_factor * self.dx
        self.x = numpy.arange(x_L - (ngz - 0.5) * self.dx, 
                              x_R + ngz * self.dx, self.dx)


class GridSpherical(object):
    """
    A finite volume style 1d grid
    
    This is identical to a standard grid, just with `x`->`r`, to be consistent
    within the higher code. We could (should) do smarter things.
    """
    def __init__(self, interval, nr, ngz=1, cfl_factor=0.9):
        """
        Set up a finite-volume style 1d grid
    
        Parameters
        ----------
        interval : [float, float]
            The problem domain [r_L, r_R]
        nr : int
            The number of (interior) grid cells.
        ngz : int
            The number of ghost zones for the boundaries
        cfl_factor : float
            dt / dr
    
        Returns
        -------
        r : array of float
            Cell centre coordinates
        dr : float
            Cell widths
        """
        
        r_L, r_R = interval
        self.r_L = r_L
        self.r_R = r_R
        self.nr = nr
        self.ngz = ngz
        self.dr = (r_R - r_L) / nr
        self.dt = cfl_factor * self.dr
        self.r = numpy.arange(r_L - (ngz - 0.5) * self.dr, 
                              r_R + ngz * self.dr, self.dr)
        
        
class Grid2d(object):
    """
    A finite volume style 2d grid
    """
    def __init__(self, interval, nx, ngz=1, cfl_factor=0.9):
        """
        Set up a finite-volume style 2d grid
    
        Parameters
        ----------
        interval : [[float, float], [float, float]]
            The problem domain [x_L, x_R]
        nx : [int, int]
            The number of (interior) grid cells.
        ngz : int
            The number of ghost zones for the boundaries
        cfl_factor : float
            dt / dx
    
        Returns
        -------
        x : array of float
            Cell centre coordinates
        dx : float
            Cell widths
        """
        
        x_L, x_R = interval[0]
        y_L, y_R = interval[1]
        self.x_L = [x_L, y_L]
        self.x_R = [x_R, y_R]
        self.nx = nx
        self.ngz = ngz
        self.dx = [(x_R - x_L) / nx[0], (y_R - y_L) / nx[1]]
        self.dt = cfl_factor * min(self.dx)
        self.x = [numpy.arange(x_L - (ngz - 0.5) * self.dx[0], 
                               x_R + ngz * self.dx[0], self.dx[0]),
                  numpy.arange(y_L - (ngz - 0.5) * self.dx[1], 
                               y_R + ngz * self.dx[1], self.dx[1])]


