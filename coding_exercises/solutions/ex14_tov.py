"""
Exercise 14: TOV star
"""

import numpy
from matplotlib import pyplot

from lib import grids, boundary_conditions, weno, tov
from systems.euler_tov import EulerTOV


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

    # Then evolve, computing the source first
    q_rhs = system.source(q)
    f = flux(q)
    alpha = abs(system.max_lambda(q))
    f_p = (f + alpha * q) / 2
    f_m = (f - alpha * q) / 2
    f_p_L = numpy.zeros_like(q)
    f_m_R = numpy.zeros_like(q)
    f_fvs = numpy.zeros_like(q)

    for i in range(g.ngz - 1, g.nr + g.ngz + 1):
        for k in range(q.shape[1]):
            # Reconstruct f plus to the right to get the state to the Left of the interface
            f_p_L[i+1, k] = weno.weno(f_p[i-(order-1):i+order, k], order)
            # Reconstruct f minus to the left to get the state to the Right of the interface
            f_m_R[i, k] = weno.weno(f_m[i-(order-1):i+(order), k][::-1], order)

    for i in range(g.ngz, g.nr + g.ngz + 1):
        f_fvs[i, :] = f_p_L[i, :] + f_m_R[i, :]

    for i in range(g.ngz, g.nr + g.ngz):
        q_rhs[i, :] += 1.0 / g.dr * (f_fvs[i, :] - f_fvs[i+1, :])

    return q_rhs


def weno_fvs(nr, t_end, system,
             g,
             bc_form="spherical_euler", weno_order=2, cfl_factor=0.9):
    """
    Solve a conservation law on [0, 10] using N gridpoints to t=t_end.

    Parameters
    ----------
    nr : int
        Number of interior gridpoints.
    t_end : float
        Final time.
    system : class
        The system defining the conservation law
    cfl_factor : float
        The ratio of dt to dx
    initial_data : function
        The initial data as a function of r
    bc_form : string
        The form of the boundary condition
    weno_order : int
        Order of accuracy of the WENO reconstruction

    Returns
    -------
    r : array of float
        Coordinates
    q : array of float
        Solution at the final time.
    """

    q = system.q0.copy()

    t = 0
    while t < t_end:
        print("Step", g.dt)
        # Compute the maximum safe timestep
        g.dt = cfl_factor * g.dr / system.max_lambda(q)
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
    return g.r, q


if __name__ == "__main__":

    nr = 100
    g = grids.GridSpherical([0, 10], nr, ngz=2, cfl_factor=0.5)
    system = EulerTOV(2, g.r, tov.initial_data_tov)

    # Solve for the TOV problem
    for t_end in [0.5]:
        r, q = weno_fvs(nr, t_end, system, g, weno_order=2, cfl_factor=0.5)
        prims = system.c2p(q)
        labels = [r"$\rho$", r"$v$", r"$\epsilon$", r"$p$", r"$c_s$"]
        fig, axes = pyplot.subplots(2, 2)
        for p, ax, label in zip(prims, axes.flatten(), labels):
            ax.plot(r, p)
            ax.set_xlabel(r"$r$")
            ax.set_ylabel(label)
#            ax.set_xlim(0, 1)
        fig.suptitle(f"WENO2, Spherical Sod problem, nr={nr}, t={t_end}")
        fig.tight_layout()
        pyplot.show()
