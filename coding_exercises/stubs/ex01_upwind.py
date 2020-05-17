"""
Exercise 1: advection with the upwind scheme

These code stubs give a lot of the boilerplate code. Look for lines

  #! To be completed

for the bits to fill in. Note that all code in the "lib" subdirectory is
complete, but a lot in the "systems" subdirectory is not.

"""

import numpy
from matplotlib import pyplot

from lib import errors, grids, boundary_conditions


def initial_data_exp_sine(x):
    """
    Set up simple initial data

    Parameters
    ----------
    x : array of float
        Cell centre coordinates

    Returns
    -------
    q : array of float
        The initial data
    """
    
    q = numpy.exp(numpy.sin(2 * numpy.pi * x))
    
    return q


def upwind_step(q, g):
    """
    Update the solution one timestep using the upwind method

    Parameters
    ----------
    q : array of float
        The solution at t^n
    g : Grid
        Information about the grid

    Returns
    -------
    q_rhs : array of float
        The update to take the solution to q at t^{n+1}
    """
    
    q_rhs = numpy.zeros_like(q)
    
    for i in range(g.ngz, g.nx + g.ngz):
        q_rhs[i] = #! To be completed

    return q_rhs


def upwind(nx, t_end, cfl_factor=0.9, initial_data=initial_data_exp_sine):
    """
    Solve the advection equation on [0, 1] using N gridpoints to t=t_end.

    Parameters
    ----------
    nx : int
        Number of interior gridpoints.
    t_end : float
        Final time.
    cfl_factor : float
        The ratio of dt to dx
    initial_data : function
        The initial data as a function of x

    Returns
    -------
    x : array of float
        Coordinates
    q : array of float
        Solution at the final time.
    """
    
    ngz = 1  # This is all we need for upwind
    
    g = grids.Grid([0, 1], nx, ngz, cfl_factor)
    q = initial_data(g.x)
    
    t = 0
    while t < t_end:
        # Update current time
        if t + g.dt > t_end:
            g.dt = t_end - t
        t += g.dt
        # Take a single step
        q += g.dt * upwind_step(q, g)
        q = boundary_conditions.boundaries(q, g)
    # Done
    return g.x, q


if __name__ == "__main__":
    
    # Solve once at low resolution to see how it does
    x, q = upwind(20, 1.0)
    x_exact = numpy.linspace(0, 1, 1000)
    q_exact = initial_data_exp_sine(x_exact)
    pyplot.plot(x_exact, q_exact, 'k-', label="Exact")
    pyplot.plot(x, q, 'bo', label=r"$N=20$")
    pyplot.xlabel(r"$x$")
    pyplot.ylabel(r"$q$")
    pyplot.legend()
    pyplot.title("Upwind, t=1")
    pyplot.show()
    
    # Solve for ten periods to see phase and dissipation errors
    x, q = upwind(20, 10.0)
    x_exact = numpy.linspace(0, 1, 1000)
    q_exact = initial_data_exp_sine(x_exact)
    pyplot.plot(x_exact, q_exact, 'k-', label="Exact")
    pyplot.plot(x, q, 'bo', label=f"$N=20$")
    pyplot.xlabel(r"$x$")
    pyplot.ylabel(r"$q$")
    pyplot.legend()
    pyplot.title("Upwind, t=10")
    pyplot.show()
    
    # Solve with various grid spacings.
    # Check convergence
    
    Ns = 2**numpy.arange(4, 10)
    errors_1norm = numpy.zeros(Ns.shape)
    errors_2norm = numpy.zeros(Ns.shape)
    errors_infnorm = numpy.zeros(Ns.shape)
    dxs = 1.0 / Ns
    
    for i, N in enumerate(Ns):
        x, q = upwind(N, 1.0)
        q_exact = initial_data_exp_sine(x)
        errors_1norm[i] = errors.error_norm(q, q_exact, dxs[i], p=1)
        errors_2norm[i] = errors.error_norm(q, q_exact, dxs[i], p=2)
        errors_infnorm[i] = errors.error_norm(q, q_exact, dxs[i], p='inf')
    
    pyplot.loglog(dxs, errors_1norm, 'x',
                  label=r"|$q_{exact} - q_{numerical}|_1$")
    pyplot.loglog(dxs, errors_2norm, 'o',
                  label=r"|$q_{exact} - q_{numerical}|_2$")
    pyplot.loglog(dxs, errors_infnorm, '^',
                  label=r"|$q_{exact} - q_{numerical}|_{\infty}$")
    pyplot.loglog(dxs, errors_1norm[-1] * (dxs / dxs[-1]),
                  label=r"$\propto (\Delta x)$")
    pyplot.xlabel(r"$\Delta x$")
    pyplot.ylabel("Errors")
    pyplot.legend()
    pyplot.show()
    