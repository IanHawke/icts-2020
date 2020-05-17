"""
Exercise 5: Euler equations with the Godunov scheme, Rusanov.
"""

import numpy
from matplotlib import pyplot

from lib import grids, boundary_conditions
from systems.euler import Euler


def initial_data_sod(system, x):
    """
    Sod shock tube data. Assumes gamma=1.4

    Parameters
    ----------
    system : Class
        The instance of the Euler class
    x : array of float
        Cell centre coordinates

    Returns
    -------
    q : array of float
        The initial data
    """
    gamma = system.gamma
    assert numpy.allclose(gamma, 1.4)
    
    rho = numpy.where(x < 0.5,
                      1.0 * numpy.ones_like(x),
                      0.125 * numpy.ones_like(x))
    v = numpy.zeros_like(x)
    p = numpy.where(x < 0.5,
                    1.0 * numpy.ones_like(x),
                    0.1 * numpy.ones_like(x))
    e = p / rho / (gamma - 1.0)
    return system.p2c(rho, v, e)
    

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
        f_rusanov[i, :] = #! To be completed
    
    for i in range(g.ngz, g.nx + g.ngz):
        q_rhs[i, :] = #! To be completed

    return q_rhs


def godunov_hlle_step(q, g, flux, c2p=None):
    """
    Update the solution one timestep using the Godunov method, HLLE flux

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
    f_HLLE = numpy.zeros_like(q)
    
    _, v, _, _, cs = c2p(q)
    
    for i in range(g.ngz, g.nx + g.ngz + 1):
        #! To be completed
        f_HLLE[i, :] = #! To be completed
    
    for i in range(g.ngz, g.nx + g.ngz):
        q_rhs[i, :] = 1.0 / g.dx * (f_HLLE[i, :] - f_HLLE[i+1, :])

    return q_rhs


def godunov(nx, t_end, system,
            cfl_factor=0.9, initial_data=initial_data_sod,
            bc_form="outflow", rp_case="rusanov"):
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
    rp_case : string
        The Riemann problem solver ('rusanov' or 'hlle')

    Returns
    -------
    x : array of float
        Coordinates
    q : array of float
        Solution at the final time.
    """
    
    ngz = 1  # This is all we need for upwind
    
    g = grids.Grid([0, 1], nx, ngz, cfl_factor)
    q = initial_data(system, g.x)
    
    t = 0
    while t < t_end:
        # Compute the maximum safe timestep
        g.dt = cfl_factor * g.dx / system.max_lambda(q)
        # Update current time
        if t + g.dt > t_end:
            g.dt = t_end - t
        t += g.dt
        # Take a single step
        if rp_case == 'rusanov':
            q += g.dt * godunov_rusanov_step(q, g, flux=system.flux)
        else:
            q += g.dt * godunov_hlle_step(q, g, flux=system.flux, c2p=system.c2p)
        q = boundary_conditions.boundaries(q, g, form=bc_form)
    # Done
    return g.x, q


if __name__ == "__main__":
    
    system = Euler(gamma=1.4)
    
    # Solve for the Sod problem: low resolution
    nx = 40
    x, q = godunov(nx, 0.2, system, initial_data=initial_data_sod)
    prims = system.c2p(q)
    labels = [r"$\rho$", r"$v$", r"$e$", r"$p$", r"$c_s$"]
    fig, axes = pyplot.subplots(2, 2)
    for p, ax, label in zip(prims, axes.flatten(), labels):
        ax.plot(x, p)
        ax.set_xlabel(r"$x$")
        ax.set_ylabel(label)
    fig.suptitle("Godunov, Euler Sod problem, nx=40")
    fig.tight_layout()
    pyplot.show()
    
    # Now get the results for high resolution
    fig, axes = pyplot.subplots(2, 2)
    for nx in [100, 400, 1600]:
        x, q = godunov(nx, 0.2, system, initial_data=initial_data_sod)
        prims = system.c2p(q)
        for p, ax, label in zip(prims, axes.flatten(), labels):
            ax.plot(x, p, label=f"nx={nx}")
            ax.set_xlabel(r"$x$")
            ax.set_ylabel(label)
    fig.suptitle("Godunov, Euler Sod problem")
    fig.tight_layout()
    pyplot.legend()
    pyplot.show()
    
    # Repeat, for HLLE
    
    # Solve for the Sod problem: low resolution
    nx = 40
    x, q = godunov(nx, 0.2, system, initial_data=initial_data_sod, rp_case="hlle")
    prims = system.c2p(q)
    labels = [r"$\rho$", r"$v$", r"$e$", r"$p$", r"$c_s$"]
    fig, axes = pyplot.subplots(2, 2)
    for p, ax, label in zip(prims, axes.flatten(), labels):
        ax.plot(x, p)
        ax.set_xlabel(r"$x$")
        ax.set_ylabel(label)
    fig.suptitle("Godunov, HLLE, Euler Sod problem, nx=40")
    fig.tight_layout()
    pyplot.show()
    
    # Now get the results for high resolution
    fig, axes = pyplot.subplots(2, 2)
    for nx in [100, 400, 1600]:
        x, q = godunov(nx, 0.2, system, initial_data=initial_data_sod, rp_case="hlle")
        prims = system.c2p(q)
        for p, ax, label in zip(prims, axes.flatten(), labels):
            ax.plot(x, p, label=f"nx={nx}")
            ax.set_xlabel(r"$x$")
            ax.set_ylabel(label)
    fig.suptitle("Godunov, HLLE, Euler Sod problem")
    fig.tight_layout()
    pyplot.legend()
    pyplot.show()