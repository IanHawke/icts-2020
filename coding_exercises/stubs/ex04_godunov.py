"""
Exercise 4: Burgers equation with the Godunov scheme, Rusanov.
"""

import numpy
from matplotlib import pyplot

from lib import grids, boundary_conditions
from systems.burgers import Burgers
from ex01_upwind import initial_data_exp_sine
from ex02_discontinuities import initial_data_hat


def godunov_rusanov_step(q, g, flux):
    """
    Update the solution one timestep using the Godunov method

    Parameters
    ----------
    q : array of float
        The solution at t^n
    g : Grid
        Information about the grid
    flux : function
        The flux function to use

    Returns
    -------
    q_rhs : array of float
        The update to take the solution to q at t^{n+1}
        
    Notes
    -----
    Compare the update formula carefully with the upwind scheme.
    The f_rusanov is now an intercell flux, so f_rusanov[i] is the flux through x[i-1/2].
    This means the indics are off-by-one compared to the upwind scheme.
    """
    
    q_rhs = numpy.zeros_like(q)
    f = flux(q)
    f_rusanov = numpy.zeros_like(q)
    
    for i in range(g.ngz, g.nx + g.ngz + 1):
        f_rusanov[i] = #! To be completed
    
    for i in range(g.ngz, g.nx + g.ngz):
        q_rhs[i] = #! To be completed

    return q_rhs


def godunov(nx, t_end, system,
            cfl_factor=0.9, initial_data=initial_data_exp_sine,
            bc_form="outflow"):
    """
    Solve a conservation law on [0, 1] using N gridpoints to t=t_end.

    Parameters
    ----------
    nx : int
        Number of interior gridpoints.
    t_end : float
        Final time.
    system : class
        The system defining the conservation law
    cfl_factor : float
        The ratio of dt to dx
    initial_data : function
        The initial data as a function of x
    bc_form : string
        The form of the boundary condition

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
        # Compute the maximum safe timestep
        g.dt = cfl_factor * g.dx / system.max_lambda(q)
        # Update current time
        if t + g.dt > t_end:
            g.dt = t_end - t
        t += g.dt
        # Take a single step
        q += g.dt * godunov_rusanov_step(q, g, flux=system.flux)
        q = boundary_conditions.boundaries(q, g, form=bc_form)
    # Done
    return g.x, q


if __name__ == "__main__":
    
    
    system = Burgers()
    
    # Solve for the smooth initial data at different times
    nx = 200
    Ts = [0.05, 0.1, 0.15, 0.2]
    mpl_formats=['b-', 'r-', 'g-', 'c-']
    for t_end, fmt in zip(Ts, mpl_formats):
        x, q = godunov(nx, t_end, system, initial_data=initial_data_exp_sine)
        pyplot.plot(x, q, fmt, label=fr"$t={t_end}$")
    pyplot.xlabel(r"$x$")
    pyplot.ylabel(r"$q$")
    pyplot.legend()
    pyplot.title("Godunov-Rusanov, Burgers, nx=200")
    pyplot.show()
    
    mpl_formats=['bo--', 'rx--', 'g^--', 'c*--']
    Ns = [10, 20, 40, 80]
    # Solve once at various resolutions to see how it does
    for nx, fmt in zip(Ns, mpl_formats):
        x, q = godunov(nx, 0.2, system, initial_data=initial_data_exp_sine)
        pyplot.plot(x, q, fmt, label=fr"$N={nx}$")
    pyplot.xlabel(r"$x$")
    pyplot.ylabel(r"$q$")
    pyplot.legend()
    pyplot.title("Godunov-Rusanov, Burgers, t=0.2")
    pyplot.show()
    
    # Repeat for the hat initial data to see the wave structure
    for nx, fmt in zip(Ns, mpl_formats):
        x, q = godunov(nx, 0.2, system, initial_data=initial_data_hat)
        pyplot.plot(x, q, fmt, label=fr"$N={nx}$")
    pyplot.xlabel(r"$x$")
    pyplot.ylabel(r"$q$")
    pyplot.legend()
    pyplot.title("Godunov-Rusanov, Burgers, t=0.2")
    pyplot.show()
    