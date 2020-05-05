"""
Exercise 9: RK2 applied to systems.
"""

import numpy
from matplotlib import pyplot

from lib import grids, boundary_conditions
from systems.euler import Euler
from systems.euler_sr import EulerSR
from ex05_euler import initial_data_sod
from ex06_euler_sr import initial_data_blast
from ex07_minmod import minmod


def minmod_rusanov_step(q, g, flux):
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
    f_rusanov = numpy.zeros_like(q)
    q_L = numpy.zeros_like(q)
    q_R = numpy.zeros_like(q)
    
    for i in range(g.ngz - 1, g.nx + g.ngz + 1):
        for k in range(q.shape[1]):
            sigma_up = q[i+1, k] - q[i, k]
            sigma_do = q[i, k] - q[i-1, k]
            sigma_bar = minmod(sigma_up, sigma_do)
            q_R[i, k] = q[i, k] - 0.5 * sigma_bar
            q_L[i+1, k] = q[i, k] + 0.5 * sigma_bar
        
    f_L = flux(q_L)
    f_R = flux(q_R)
    for i in range(g.ngz, g.nx + g.ngz + 1):
        f_rusanov[i, :] = (f_L[i, :] + f_R[i, :] + g.dx / g.dt * (q_L[i, :] - q_R[i, :])) / 2
    
    for i in range(g.ngz, g.nx + g.ngz):
        q_rhs[i, :] = 1.0 / g.dx * (f_rusanov[i, :] - f_rusanov[i+1, :])

    return q_rhs


def rk2_minmod(nx, t_end, system,
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

    Returns
    -------
    x : array of float
        Coordinates
    q : array of float
        Solution at the final time.
    """
    
    ngz = 2  # A second ghostzone is needed for slope limiting
    
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
        if rp_case == 'rusanov':
            q1 = q + g.dt * minmod_rusanov_step(q, g, flux=system.flux)
        else:
            raise NotImplementedError
        q1 = boundary_conditions.boundaries(q1, g, form=bc_form)
        # Take second step
        if rp_case == 'rusanov':
            q = 0.5 * (q + q1 + g.dt * minmod_rusanov_step(q1, g, flux=system.flux))
        else:
            raise NotImplementedError
        q = boundary_conditions.boundaries(q, g, form=bc_form)
    # Done
    return g.x, q


if __name__ == "__main__":
    
    system = Euler(gamma=1.4)
    
    # Solve for the Sod problem: low resolution
    nx = 40
    x, q = rk2_minmod(nx, 0.2, system, initial_data=initial_data_sod)
    prims = system.c2p(q)
    labels = [r"$\rho$", r"$v$", r"$e$", r"$p$", r"$c_s$"]
    fig, axes = pyplot.subplots(2, 2)
    for p, ax, label in zip(prims, axes.flatten(), labels):
        ax.plot(x, p)
        ax.set_xlabel(r"$x$")
        ax.set_ylabel(label)
    fig.suptitle("RK2, Minmod-Rusanov, Euler Sod problem, nx=40")
    fig.tight_layout()
    pyplot.show()
    
    # Now get the results for high resolution
    fig, axes = pyplot.subplots(2, 2)
    for nx in [100, 400, 1600]:
        x, q = rk2_minmod(nx, 0.2, system, initial_data=initial_data_sod)
        prims = system.c2p(q)
        for p, ax, label in zip(prims, axes.flatten(), labels):
            ax.plot(x, p, label=f"nx={nx}")
            ax.set_xlabel(r"$x$")
            ax.set_ylabel(label)
    fig.suptitle("RK2, Minmod-Rusanov, Euler Sod problem")
    fig.tight_layout()
    pyplot.legend()
    pyplot.show()
    
    # Now do the SR case
    
    system = EulerSR(gamma=1.4)
    
    # Solve for the Sod problem: low resolution
    nx = 40
    x, q = rk2_minmod(nx, 0.4, system, initial_data=initial_data_sod)
    prims = system.c2p(q)
    labels = [r"$\rho$", r"$v$", r"$\epsilon$", r"$p$", r"$W$", r"$h$", r"$c_s$"]
    fig, axes = pyplot.subplots(2, 3)
    for p, ax, label in zip(prims, axes.flatten(), labels):
        ax.plot(x, p)
        ax.set_xlabel(r"$x$")
        ax.set_ylabel(label)
    fig.suptitle("RK2, Minmod-Rusanov, SR Euler Sod problem, nx=40")
    fig.tight_layout()
    pyplot.show()
    
    # Now get the results for high resolution
    fig, axes = pyplot.subplots(2, 3)
    for nx in [100, 200, 400]:
        x, q = rk2_minmod(nx, 0.4, system, initial_data=initial_data_sod)
        prims = system.c2p(q)
        for p, ax, label in zip(prims, axes.flatten(), labels):
            ax.plot(x, p, label=f"nx={nx}")
            ax.set_xlabel(r"$x$")
            ax.set_ylabel(label)
    fig.suptitle("RK2, Minmod-Rusanov, SR Euler Sod problem")
    fig.tight_layout()
    pyplot.legend()
    pyplot.show()
    
    # Repeat, for blast wave
    system = EulerSR(gamma=5/3)
    nx = 40
    x, q = rk2_minmod(nx, 0.4, system, initial_data=initial_data_blast)
    prims = system.c2p(q)
    labels = [r"$\rho$", r"$v$", r"$\epsilon$", r"$p$", r"$W$", r"$h$", r"$c_s$"]
    fig, axes = pyplot.subplots(2, 3)
    for p, ax, label in zip(prims, axes.flatten(), labels):
        ax.plot(x, p)
        ax.set_xlabel(r"$x$")
        ax.set_ylabel(label)
    fig.suptitle("RK2, Minmod-Rusanov, SR Euler blast wave problem, nx=40")
    fig.tight_layout()
    pyplot.show()
    
    # Now get the results for high resolution
    fig, axes = pyplot.subplots(2, 3)
    for nx in [100, 200, 400]:
        x, q = rk2_minmod(nx, 0.4, system, initial_data=initial_data_blast)
        prims = system.c2p(q)
        for p, ax, label in zip(prims, axes.flatten(), labels):
            ax.plot(x, p, label=f"nx={nx}")
            ax.set_xlabel(r"$x$")
            ax.set_ylabel(label)
    fig.suptitle("RK2, Minmod-Rusanov, SR Euler blast wave problem")
    fig.tight_layout()
    pyplot.legend()
    pyplot.show()
    