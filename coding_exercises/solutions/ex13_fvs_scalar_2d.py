"""
Exercise 13: Flux vector splitting, 2d
"""

import numpy
from matplotlib import pyplot

from lib import grids, boundary_conditions, errors, weno
from systems.advection import Advection2d


def initial_data_gauss(x):
    """
    A gaussian pulse centred at 0.5
    """
    q = numpy.zeros((len(x[0]), len(x[1])))
    for i in range(len(x[0])):
        for j in range(len(x[1])):
            r = numpy.sqrt((x[0][i] - 0.5)**2 + (x[1][j] - 0.5)**2)
            q[i, j] = numpy.exp(-100*r**2)
    return q


def weno_fvs_1d_step(q, g, flux, dirn, order=2):
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
    f = flux(q, dirn)
    alpha = abs(system.max_lambda(q))
    f_p = (f + alpha * q) / 2
    f_m = (f - alpha * q) / 2
    f_p_L = numpy.zeros_like(q)
    f_m_R = numpy.zeros_like(q)
    f_fvs = numpy.zeros_like(q)
    
    for i in range(g.ngz - 1, g.nx[dirn] + g.ngz + 1):
        # Reconstruct f plus to the right to get the state to the Left of the interface
        f_p_L[i+1] = weno.weno(f_p[i-(order-1):i+order], order)
        # Reconstruct f minus to the left to get the state to the Right of the interface
        f_m_R[i] = weno.weno(f_m[i-(order-1):i+(order)][::-1], order)
        
    for i in range(g.ngz, g.nx[dirn] + g.ngz + 1):
        f_fvs[i] = f_p_L[i] + f_m_R[i]
    
    for i in range(g.ngz, g.nx[dirn] + g.ngz):
        q_rhs[i] = 1.0 / g.dx[dirn] * (f_fvs[i] - f_fvs[i+1])

    return q_rhs

def weno_fvs_step(q, g, flux, order=2):
    """
    A single step in 2d using FVS.
    """
    q_rhs = numpy.zeros_like(q)
    # Compute in x direction
    for j in range(g.ngz, g.nx[1]+g.ngz):
        q_rhs[:, j] += weno_fvs_1d_step(q[:, j], g, flux, 0, order)
    # Compute in y direction
    for i in range(g.ngz, g.nx[0]+g.ngz):
        q_rhs[i, :] += weno_fvs_1d_step(q[i, :], g, flux, 1, order)
    
    return q_rhs

def weno_fvs(nx, t_end, system,
               cfl_factor=0.4, initial_data=initial_data_gauss,
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
    
    g = grids.Grid2d([[0, 1], [0,1]], [nx, nx], ngz, cfl_factor)
    q = initial_data(g.x)
    
    t = 0
    while t < t_end:
        # Compute the maximum safe timestep
        g.dt = cfl_factor * min(g.dx) / system.max_lambda(q)
        # Update current time
        if t + g.dt > t_end:
            g.dt = t_end - t
        t += g.dt
        # Take first step
        q1 = q + g.dt * weno_fvs_step(q, g, flux=system.flux, order=weno_order)
        q1 = boundary_conditions.boundaries2d(q1, g, form=bc_form)
        # Second RK2 step
        q2 = (3 * q + q1 + g.dt * weno_fvs_step(q1, g, flux=system.flux, order=weno_order)) / 4
        q2 = boundary_conditions.boundaries2d(q2, g, form=bc_form)
        # Second RK2 step
        q = (q + 2 * q2 + 2 * g.dt * weno_fvs_step(q2, g, flux=system.flux, order=weno_order)) / 3
        q = boundary_conditions.boundaries2d(q, g, form=bc_form)
    # Done
    return g.x, q


if __name__ == "__main__":
    
    system = Advection2d()
    
    # Solve once at low resolution to see how it does
    x, q = weno_fvs(20, 1.0, system, bc_form="periodic")
    x_exact = numpy.linspace(0, 1, 1000)
    q_exact = initial_data_gauss([x_exact, x_exact])
    pyplot.plot(x_exact, q_exact[:, 500], 'k-', label="Exact")
    pyplot.plot(x[0], q[:, 12], 'bo', label=r"$N=20$")
    pyplot.xlabel(r"$x$")
    pyplot.ylabel(r"$q$")
    pyplot.legend()
    pyplot.title("FVS, t=1")
    pyplot.show()
    
    fig, axes = pyplot.subplots(1, 2)
    im = axes[0].imshow(q_exact, vmin=0, vmax=1)
    im = axes[1].imshow(q, vmin=0, vmax=1)
    fig.colorbar(im, ax=axes.ravel().tolist())
    pyplot.show()
    
    
    # Solve once at high resolution to see how it does
    x, q = weno_fvs(100, 1.0, system, bc_form="periodic")
    x_exact = numpy.linspace(0, 1, 1000)
    q_exact = initial_data_gauss([x_exact, x_exact])
    pyplot.plot(x_exact, q_exact[:, 500], 'k-', label="Exact")
    pyplot.plot(x[0], q[:, 52], 'bo', label=r"$N=100$")
    pyplot.xlabel(r"$x$")
    pyplot.ylabel(r"$q$")
    pyplot.legend()
    pyplot.title("FVS, t=1")
    pyplot.show()
    
    fig, axes = pyplot.subplots(1, 2)
    im = axes[0].imshow(q_exact, vmin=0, vmax=1)
    im = axes[1].imshow(q, vmin=0, vmax=1)
    fig.colorbar(im, ax=axes.ravel().tolist())
    pyplot.show()
    
    
    # Repeat at higher order
    # Solve once at low resolution to see how it does
    x, q = weno_fvs(20, 1.0, system, bc_form="periodic", weno_order=3)
    x_exact = numpy.linspace(0, 1, 1000)
    q_exact = initial_data_gauss([x_exact, x_exact])
    pyplot.plot(x_exact, q_exact[:, 500], 'k-', label="Exact")
    pyplot.plot(x[0], q[:, 12], 'bo', label=r"$N=20$")
    pyplot.xlabel(r"$x$")
    pyplot.ylabel(r"$q$")
    pyplot.legend()
    pyplot.title("FVS, t=1")
    pyplot.show()
    
    fig, axes = pyplot.subplots(1, 2)
    im = axes[0].imshow(q_exact, vmin=0, vmax=1)
    im = axes[1].imshow(q, vmin=0, vmax=1)
    fig.colorbar(im, ax=axes.ravel().tolist())
    pyplot.show()
    
    
    # Solve once at high resolution to see how it does
    x, q = weno_fvs(100, 1.0, system, bc_form="periodic", weno_order=3)
    x_exact = numpy.linspace(0, 1, 1000)
    q_exact = initial_data_gauss([x_exact, x_exact])
    pyplot.plot(x_exact, q_exact[:, 500], 'k-', label="Exact")
    pyplot.plot(x[0], q[:, 52], 'bo', label=r"$N=100$")
    pyplot.xlabel(r"$x$")
    pyplot.ylabel(r"$q$")
    pyplot.legend()
    pyplot.title("FVS, t=1")
    pyplot.show()
    
    fig, axes = pyplot.subplots(1, 2)
    im = axes[0].imshow(q_exact, vmin=0, vmax=1)
    im = axes[1].imshow(q, vmin=0, vmax=1)
    fig.colorbar(im, ax=axes.ravel().tolist())
    pyplot.show()
    
    