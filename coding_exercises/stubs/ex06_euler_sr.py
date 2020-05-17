"""
Exercise 6: SR Euler equations with the Godunov scheme, Rusanov.
"""

import numpy
from matplotlib import pyplot

from lib import grids, boundary_conditions
from systems.euler_sr import EulerSR
from ex05_euler import initial_data_sod, godunov_rusanov_step


def initial_data_blast(system, x):
    """
    Sod shock tube data. Assumes gamma=5/3

    Parameters
    ----------
    system : Class
        The instance of the EulerSR class
    x : array of float
        Cell centre coordinates

    Returns
    -------
    q : array of float
        The initial data
    """
    gamma = system.gamma
    assert numpy.allclose(gamma, 5/3)
    
    rho = numpy.where(x < 0.5,
                      1.0 * numpy.ones_like(x),
                      0.125 * numpy.ones_like(x))
    v = numpy.zeros_like(x)
    p = numpy.where(x < 0.5,
                    1000.0 * numpy.ones_like(x),
                    0.01 * numpy.ones_like(x))
    epsilon = p / rho / (gamma - 1.0)
    return system.p2c(rho, v, epsilon)
    

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
            raise NotImplementedError
        q = boundary_conditions.boundaries(q, g, form=bc_form)
    # Done
    return g.x, q


if __name__ == "__main__":
    
    system = EulerSR(gamma=1.4)
    
    # Set up some "random" data
    
    nx = 1000
    rho = numpy.ones(nx) + 0.1 * numpy.random.randn(nx)
    v = 0.7 * 2 * (numpy.random.rand(nx) - 0.5)
    epsilon = numpy.ones(nx) + 0.3 * numpy.random.randn(nx)
    q = system.p2c(rho, v, epsilon)
    rho_c2p, v_c2p, epsilon_c2p, _, _, _, _ = system.c2p(q)
    print("Check P->C->P for rho, v, epsilon:",
          numpy.allclose(rho, rho_c2p),
          numpy.allclose(v, v_c2p),
          numpy.allclose(epsilon, epsilon_c2p))
    
    # Solve for the Sod problem: low resolution
    nx = 40
    x, q = godunov(nx, 0.4, system, initial_data=initial_data_sod)
    prims = system.c2p(q)
    labels = [r"$\rho$", r"$v$", r"$\epsilon$", r"$p$", r"$W$", r"$h$", r"$c_s$"]
    fig, axes = pyplot.subplots(2, 3)
    for p, ax, label in zip(prims, axes.flatten(), labels):
        ax.plot(x, p)
        ax.set_xlabel(r"$x$")
        ax.set_ylabel(label)
    fig.suptitle("Godunov, SR Euler Sod problem, nx=40")
    fig.tight_layout()
    pyplot.show()
    
    # Now get the results for high resolution
    fig, axes = pyplot.subplots(2, 3)
    for nx in [100, 200, 400]:
        x, q = godunov(nx, 0.4, system, initial_data=initial_data_sod)
        prims = system.c2p(q)
        for p, ax, label in zip(prims, axes.flatten(), labels):
            ax.plot(x, p, label=f"nx={nx}")
            ax.set_xlabel(r"$x$")
            ax.set_ylabel(label)
    fig.suptitle("Godunov, SR Euler Sod problem")
    fig.tight_layout()
    pyplot.legend()
    pyplot.show()
    
    # Repeat, for blast wave
    system = EulerSR(gamma=5/3)
    nx = 40
    x, q = godunov(nx, 0.4, system, initial_data=initial_data_blast)
    prims = system.c2p(q)
    labels = [r"$\rho$", r"$v$", r"$\epsilon$", r"$p$", r"$W$", r"$h$", r"$c_s$"]
    fig, axes = pyplot.subplots(2, 3)
    for p, ax, label in zip(prims, axes.flatten(), labels):
        ax.plot(x, p)
        ax.set_xlabel(r"$x$")
        ax.set_ylabel(label)
    fig.suptitle("Godunov, SR Euler blast wave problem, nx=40")
    fig.tight_layout()
    pyplot.show()
    
    # Now get the results for high resolution
    fig, axes = pyplot.subplots(2, 3)
    for nx in [100, 200, 400]:
        x, q = godunov(nx, 0.4, system, initial_data=initial_data_blast)
        prims = system.c2p(q)
        for p, ax, label in zip(prims, axes.flatten(), labels):
            ax.plot(x, p, label=f"nx={nx}")
            ax.set_xlabel(r"$x$")
            ax.set_ylabel(label)
    fig.suptitle("Godunov, SR Euler blast wave problem")
    fig.tight_layout()
    pyplot.legend()
    pyplot.show()
    