"""
Exercise 10: Flux vector splitting
"""

import numpy
from matplotlib import pyplot

from lib import grids, boundary_conditions, errors, weno
from systems.advection import Advection
from ex01_upwind import initial_data_exp_sine
from ex07_minmod import minmod


def minmod_fvs_step(q, g, flux):
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
    f = flux(q)
    alpha = abs(system.max_lambda(q))
    f_p = #! To be completed
    f_m = #! To be completed
    f_p_L = numpy.zeros_like(q)
    f_m_R = numpy.zeros_like(q)
    f_fvs = numpy.zeros_like(q)
    
    for i in range(g.ngz - 1, g.nx + g.ngz + 1):
        # Reconstruct f plus to the right to get the state to the Left of the interface
        sigma_up = f_p[i+1] - f_p[i]
        sigma_do = f_p[i] - f_p[i-1]
        sigma_bar = minmod(sigma_up, sigma_do)
        f_p_L[i+1] = f_p[i] + 0.5 * sigma_bar
        # Reconstruct f minus to the left to get the state to the Right of the interface
        sigma_up = f_m[i+1] - f_m[i]
        sigma_do = f_m[i] - f_m[i-1]
        sigma_bar = minmod(sigma_up, sigma_do)
        f_m_R[i] = f_m[i] - 0.5 * sigma_bar
        
    for i in range(g.ngz, g.nx + g.ngz + 1):
        f_fvs[i] = #! To be completed
    
    for i in range(g.ngz, g.nx + g.ngz):
        q_rhs[i] = #! To be completed

    return q_rhs


def minmod_fvs(nx, t_end, system,
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
    
    ngz = 2  # A second ghostzone is needed for slope limiting
    
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
        # Take first step
        q1 = q + g.dt * minmod_fvs_step(q, g, flux=system.flux)
        q1 = boundary_conditions.boundaries(q1, g, form=bc_form)
        # Second RK2 step
        q = 0.5 * (q + q1 + g.dt * minmod_fvs_step(q1, g, flux=system.flux))
        q = boundary_conditions.boundaries(q, g, form=bc_form)
    # Done
    return g.x, q


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
    alpha = abs(system.max_lambda(q))
    f_p = #! To be completed
    f_m = #! To be completed
    f_p_L = numpy.zeros_like(q)
    f_m_R = numpy.zeros_like(q)
    f_fvs = numpy.zeros_like(q)
    
    for i in range(g.ngz - 1, g.nx + g.ngz + 1):
        # Reconstruct f plus to the right to get the state to the Left of the interface
        f_p_L[i+1] = weno.weno(f_p[i-(order-1):i+order], order)
        # Reconstruct f minus to the left to get the state to the Right of the interface
        f_m_R[i] = weno.weno(f_m[i-(order-1):i+(order)][::-1], order)
        
    for i in range(g.ngz, g.nx + g.ngz + 1):
        f_fvs[i] = #! To be completed
    
    for i in range(g.ngz, g.nx + g.ngz):
        q_rhs[i] = #! To be completed

    return q_rhs


def weno_fvs(nx, t_end, system,
               cfl_factor=0.7, initial_data=initial_data_exp_sine,
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
    q = initial_data(g.x)
    
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
    
    system = Advection()
    
    # Solve once at low resolution to see how it does
    x, q = minmod_fvs(20, 1.0, system, bc_form="periodic")
    x_exact = numpy.linspace(0, 1, 1000)
    q_exact = initial_data_exp_sine(x_exact)
    pyplot.plot(x_exact, q_exact, 'k-', label="Exact")
    pyplot.plot(x, q, 'bo', label=r"$N=20$")
    pyplot.xlabel(r"$x$")
    pyplot.ylabel(r"$q$")
    pyplot.legend()
    pyplot.title("FVS, t=1")
    pyplot.show()
    
    # Solve for ten periods to see phase and dissipation errors
    x, q = minmod_fvs(20, 10.0, system, bc_form="periodic")
    x_exact = numpy.linspace(0, 1, 1000)
    q_exact = initial_data_exp_sine(x_exact)
    pyplot.plot(x_exact, q_exact, 'k-', label="Exact")
    pyplot.plot(x, q, 'bo', label=f"$N=20$")
    pyplot.xlabel(r"$x$")
    pyplot.ylabel(r"$q$")
    pyplot.legend()
    pyplot.title("FVS, t=10")
    pyplot.show()
    
    # Solve with various grid spacings.
    # Check convergence
    
    Ns = 2**numpy.arange(4, 10)
    errors_1norm = numpy.zeros(Ns.shape)
    errors_2norm = numpy.zeros(Ns.shape)
    errors_infnorm = numpy.zeros(Ns.shape)
    dxs = 1.0 / Ns
    
    for i, N in enumerate(Ns):
        x, q = minmod_fvs(N, 1.0, system, bc_form="periodic")
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
    pyplot.loglog(dxs, errors_1norm[-1] * (dxs / dxs[-1])**1.5,
                  label=r"$\propto (\Delta x)^{3/2}$")
    pyplot.xlabel(r"$\Delta x$")
    pyplot.ylabel("Errors")
    pyplot.legend()
    pyplot.show()
    
    # Now do it with WENO2, RK3
    
    # Solve once at low resolution to see how it does
    x, q = weno_fvs(20, 1.0, system, bc_form="periodic", weno_order=2)
    x_exact = numpy.linspace(0, 1, 1000)
    q_exact = initial_data_exp_sine(x_exact)
    pyplot.plot(x_exact, q_exact, 'k-', label="Exact")
    pyplot.plot(x, q, 'bo', label=r"$N=20$")
    pyplot.xlabel(r"$x$")
    pyplot.ylabel(r"$q$")
    pyplot.legend()
    pyplot.title("WENO2 FVS, t=1")
    pyplot.show()
    
    # Solve for ten periods to see phase and dissipation errors
    x, q = weno_fvs(20, 10.0, system, bc_form="periodic", weno_order=2)
    x_exact = numpy.linspace(0, 1, 1000)
    q_exact = initial_data_exp_sine(x_exact)
    pyplot.plot(x_exact, q_exact, 'k-', label="Exact")
    pyplot.plot(x, q, 'bo', label=f"$N=20$")
    pyplot.xlabel(r"$x$")
    pyplot.ylabel(r"$q$")
    pyplot.legend()
    pyplot.title("WENO2 FVS, t=10")
    pyplot.show()
    
    # Solve with various grid spacings.
    # Check convergence
    
    Ns = 2**numpy.arange(3, 9)
    errors_1norm = numpy.zeros(Ns.shape)
    errors_2norm = numpy.zeros(Ns.shape)
    errors_infnorm = numpy.zeros(Ns.shape)
    dxs = 1.0 / Ns
    
    for i, N in enumerate(Ns):
        x, q = weno_fvs(N, 1.0, system, bc_form="periodic", weno_order=2)
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
    pyplot.loglog(dxs, errors_1norm[-1] * (dxs / dxs[-1])**2,
                  label=r"$\propto (\Delta x)^{2}$")
    pyplot.xlabel(r"$\Delta x$")
    pyplot.ylabel("Errors")
    pyplot.legend()
    pyplot.show()
    
    # Now do it with WENO3, RK3
    
    # Solve once at low resolution to see how it does
    x, q = weno_fvs(20, 1.0, system, bc_form="periodic", weno_order=3)
    x_exact = numpy.linspace(0, 1, 1000)
    q_exact = initial_data_exp_sine(x_exact)
    pyplot.plot(x_exact, q_exact, 'k-', label="Exact")
    pyplot.plot(x, q, 'bo', label=r"$N=20$")
    pyplot.xlabel(r"$x$")
    pyplot.ylabel(r"$q$")
    pyplot.legend()
    pyplot.title("WENO3 FVS, t=1")
    pyplot.show()
    
    # Solve for ten periods to see phase and dissipation errors
    x, q = weno_fvs(20, 10.0, system, bc_form="periodic", weno_order=3)
    x_exact = numpy.linspace(0, 1, 1000)
    q_exact = initial_data_exp_sine(x_exact)
    pyplot.plot(x_exact, q_exact, 'k-', label="Exact")
    pyplot.plot(x, q, 'bo', label=f"$N=20$")
    pyplot.xlabel(r"$x$")
    pyplot.ylabel(r"$q$")
    pyplot.legend()
    pyplot.title("WENO3 FVS, t=10")
    pyplot.show()
    
    # Solve with various grid spacings.
    # Check convergence
    
    Ns = 2**numpy.arange(3, 9)
    errors_1norm = numpy.zeros(Ns.shape)
    errors_2norm = numpy.zeros(Ns.shape)
    errors_infnorm = numpy.zeros(Ns.shape)
    dxs = 1.0 / Ns
    
    for i, N in enumerate(Ns):
        x, q = weno_fvs(N, 1.0, system, bc_form="periodic", weno_order=3)
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
    pyplot.loglog(dxs, errors_1norm[-1] * (dxs / dxs[-1])**3,
                  label=r"$\propto (\Delta x)^{3}$")
    pyplot.xlabel(r"$\Delta x$")
    pyplot.ylabel("Errors")
    pyplot.legend()
    pyplot.show()
    