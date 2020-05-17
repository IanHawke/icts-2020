"""
Exercise 11: Flux vector splitting applied to Euler equations
"""

import numpy
from matplotlib import pyplot

from lib import grids, boundary_conditions, weno
from systems.euler import Euler
from systems.euler_sr import EulerSR
from ex05_euler import initial_data_sod
from ex06_euler_sr import initial_data_blast


def weno_fvs_step(q, g, flux, order=2):
    """
    Update the solution one timestep using the slope limited minmod method

    Parameters
    ----------
    q : array of float
        The solution at t^n
    g : Grid
        Information about the grid
    flux : function
        The flux function to use
    order : int
        WENO order to use

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
    alpha = #! To be completed
    f_p = #! To be completed
    f_m = #! To be completed
    f_p_L = numpy.zeros_like(q)
    f_m_R = numpy.zeros_like(q)
    f_fvs = numpy.zeros_like(q)
    
    for i in range(g.ngz - 1, g.nx + g.ngz + 1):
        for k in range(q.shape[1]):
            # Reconstruct f plus to the right to get the state to the Left of the interface
            f_p_L[i+1, k] = #! To be completed
            # Reconstruct f minus to the left to get the state to the Right of the interface
            f_m_R[i, k] = #! To be completed
        
    for i in range(g.ngz, g.nx + g.ngz + 1):
        f_fvs[i, :] = #! To be completed
    
    for i in range(g.ngz, g.nx + g.ngz):
        q_rhs[i, :] = #! To be completed

    return q_rhs


def weno_fvs(nx, t_end, system,
               cfl_factor=0.7, initial_data=initial_data_sod,
               bc_form="outflow", weno_order=2):
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
    weno_order : int
        Order of accuracy of the WENO reconstruction

    Returns
    -------
    x : array of float
        Coordinates
    q : array of float
        Solution at the final time.
    """
    
    ngz = weno_order  # WENO schemes of order 2k-1 need k ghostzones
    
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
        # Take first step
        q1 = q + g.dt * weno_fvs_step(q, g, flux=system.flux, order=weno_order)
        q1 = boundary_conditions.boundaries(q1, g, form=bc_form)
        # Second RK2 step
        q2 = (3 * q + q1 + g.dt * weno_fvs_step(q1, g, flux=system.flux, order=weno_order)) / 4
        q2 = boundary_conditions.boundaries(q2, g, form=bc_form)
        # Second RK2 step
        q = (q + 2 * q2 + 2 * g.dt * weno_fvs_step(q2, g, flux=system.flux, order=weno_order)) / 3
        q = boundary_conditions.boundaries(q, g, form=bc_form)
    # Done
    return g.x, q


if __name__ == "__main__":
    
    system = Euler(gamma=1.4)
    
    # Solve for the Sod problem: low resolution
    nx = 40
    x, q = weno_fvs(nx, 0.2, system, initial_data=initial_data_sod, weno_order=2)
    prims = system.c2p(q)
    labels = [r"$\rho$", r"$v$", r"$e$", r"$p$", r"$c_s$"]
    fig, axes = pyplot.subplots(2, 2)
    for p, ax, label in zip(prims, axes.flatten(), labels):
        ax.plot(x, p)
        ax.set_xlabel(r"$x$")
        ax.set_ylabel(label)
    fig.suptitle("WENO2, Euler Sod problem, nx=40")
    fig.tight_layout()
    pyplot.show()
    
    # Now get the results for high resolution
    fig, axes = pyplot.subplots(2, 2)
    for nx in [100, 200, 400]:
        x, q = weno_fvs(nx, 0.2, system, initial_data=initial_data_sod, weno_order=2)
        prims = system.c2p(q)
        for p, ax, label in zip(prims, axes.flatten(), labels):
            ax.plot(x, p, label=f"nx={nx}")
            ax.set_xlabel(r"$x$")
            ax.set_ylabel(label)
    fig.suptitle("WENO2, Euler Sod problem")
    fig.tight_layout()
    pyplot.legend()
    pyplot.show()
    
    # Now do the SR case
    
    system = EulerSR(gamma=1.4)
    
    # Solve for the Sod problem: low resolution
    nx = 40
    x, q = weno_fvs(nx, 0.4, system, initial_data=initial_data_sod, weno_order=2)
    prims = system.c2p(q)
    labels = [r"$\rho$", r"$v$", r"$\epsilon$", r"$p$", r"$W$", r"$h$", r"$c_s$"]
    fig, axes = pyplot.subplots(2, 3)
    for p, ax, label in zip(prims, axes.flatten(), labels):
        ax.plot(x, p)
        ax.set_xlabel(r"$x$")
        ax.set_ylabel(label)
    fig.suptitle("WENO2, SR Euler Sod problem, nx=40")
    fig.tight_layout()
    pyplot.show()
    
    # Now get the results for high resolution
    fig, axes = pyplot.subplots(2, 3)
    for nx in [100, 200, 400]:
        x, q = weno_fvs(nx, 0.4, system, initial_data=initial_data_sod, weno_order=2)
        prims = system.c2p(q)
        for p, ax, label in zip(prims, axes.flatten(), labels):
            ax.plot(x, p, label=f"nx={nx}")
            ax.set_xlabel(r"$x$")
            ax.set_ylabel(label)
    fig.suptitle("WENO2, SR Euler Sod problem")
    fig.tight_layout()
    pyplot.legend()
    pyplot.show()
    
    # Repeat, for blast wave
    system = EulerSR(gamma=5/3)
    nx = 40
    x, q = weno_fvs(nx, 0.4, system, initial_data=initial_data_blast, weno_order=2)
    prims = system.c2p(q)
    labels = [r"$\rho$", r"$v$", r"$\epsilon$", r"$p$", r"$W$", r"$h$", r"$c_s$"]
    fig, axes = pyplot.subplots(2, 3)
    for p, ax, label in zip(prims, axes.flatten(), labels):
        ax.plot(x, p)
        ax.set_xlabel(r"$x$")
        ax.set_ylabel(label)
    fig.suptitle("WENO2, SR Euler blast wave problem, nx=40")
    fig.tight_layout()
    pyplot.show()
    
    # Now get the results for high resolution
    fig, axes = pyplot.subplots(2, 3)
    for nx in [100, 200, 400]:
        x, q = weno_fvs(nx, 0.4, system, initial_data=initial_data_blast, weno_order=2)
        prims = system.c2p(q)
        for p, ax, label in zip(prims, axes.flatten(), labels):
            ax.plot(x, p, label=f"nx={nx}")
            ax.set_xlabel(r"$x$")
            ax.set_ylabel(label)
    fig.suptitle("WENO2, SR Euler blast wave problem")
    fig.tight_layout()
    pyplot.legend()
    pyplot.show()
    
    # Now do WENO3
    
    system = Euler(gamma=1.4)
    
    # Solve for the Sod problem: low resolution
    nx = 40
    x, q = weno_fvs(nx, 0.2, system, initial_data=initial_data_sod, weno_order=3)
    prims = system.c2p(q)
    labels = [r"$\rho$", r"$v$", r"$e$", r"$p$", r"$c_s$"]
    fig, axes = pyplot.subplots(2, 2)
    for p, ax, label in zip(prims, axes.flatten(), labels):
        ax.plot(x, p)
        ax.set_xlabel(r"$x$")
        ax.set_ylabel(label)
    fig.suptitle("WENO3, Euler Sod problem, nx=40")
    fig.tight_layout()
    pyplot.show()
    
    # Now get the results for high resolution
    fig, axes = pyplot.subplots(2, 2)
    for nx in [100, 200, 400]:
        x, q = weno_fvs(nx, 0.2, system, initial_data=initial_data_sod, weno_order=3)
        prims = system.c2p(q)
        for p, ax, label in zip(prims, axes.flatten(), labels):
            ax.plot(x, p, label=f"nx={nx}")
            ax.set_xlabel(r"$x$")
            ax.set_ylabel(label)
    fig.suptitle("WENO3, Euler Sod problem")
    fig.tight_layout()
    pyplot.legend()
    pyplot.show()
    
    # Now do the SR case
    
    system = EulerSR(gamma=1.4)
    
    # Solve for the Sod problem: low resolution
    nx = 40
    x, q = weno_fvs(nx, 0.4, system, initial_data=initial_data_sod, weno_order=3)
    prims = system.c2p(q)
    labels = [r"$\rho$", r"$v$", r"$\epsilon$", r"$p$", r"$W$", r"$h$", r"$c_s$"]
    fig, axes = pyplot.subplots(2, 3)
    for p, ax, label in zip(prims, axes.flatten(), labels):
        ax.plot(x, p)
        ax.set_xlabel(r"$x$")
        ax.set_ylabel(label)
    fig.suptitle("WENO3, SR Euler Sod problem, nx=40")
    fig.tight_layout()
    pyplot.show()
    
    # Now get the results for high resolution
    fig, axes = pyplot.subplots(2, 3)
    for nx in [100, 200, 400]:
        x, q = weno_fvs(nx, 0.4, system, initial_data=initial_data_sod, weno_order=3)
        prims = system.c2p(q)
        for p, ax, label in zip(prims, axes.flatten(), labels):
            ax.plot(x, p, label=f"nx={nx}")
            ax.set_xlabel(r"$x$")
            ax.set_ylabel(label)
    fig.suptitle("WENO3, SR Euler Sod problem")
    fig.tight_layout()
    pyplot.legend()
    pyplot.show()
    
    # Repeat, for blast wave
    system = EulerSR(gamma=5/3)
    nx = 40
    x, q = weno_fvs(nx, 0.4, system, initial_data=initial_data_blast, weno_order=3)
    prims = system.c2p(q)
    labels = [r"$\rho$", r"$v$", r"$\epsilon$", r"$p$", r"$W$", r"$h$", r"$c_s$"]
    fig, axes = pyplot.subplots(2, 3)
    for p, ax, label in zip(prims, axes.flatten(), labels):
        ax.plot(x, p)
        ax.set_xlabel(r"$x$")
        ax.set_ylabel(label)
    fig.suptitle("WENO3, SR Euler blast wave problem, nx=40")
    fig.tight_layout()
    pyplot.show()
    
    # Now get the results for high resolution
    fig, axes = pyplot.subplots(2, 3)
    for nx in [100, 200, 400]:
        x, q = weno_fvs(nx, 0.4, system, initial_data=initial_data_blast, weno_order=3)
        prims = system.c2p(q)
        for p, ax, label in zip(prims, axes.flatten(), labels):
            ax.plot(x, p, label=f"nx={nx}")
            ax.set_xlabel(r"$x$")
            ax.set_ylabel(label)
    fig.suptitle("WENO3, SR Euler blast wave problem")
    fig.tight_layout()
    pyplot.legend()
    pyplot.show()
    