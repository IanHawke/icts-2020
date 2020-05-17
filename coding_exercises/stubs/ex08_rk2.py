"""
Exercise 8: RK2.
"""

import numpy
from matplotlib import pyplot

from lib import grids, boundary_conditions, errors
from systems.advection import Advection
from systems.burgers import Burgers
from ex01_upwind import initial_data_exp_sine
from ex02_discontinuities import initial_data_hat
from ex07_minmod import minmod, minmod_rusanov_step


def minmod_upwind_step(q, g, flux):
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
    f_upwind = numpy.zeros_like(q)
    q_L = numpy.zeros_like(q)
    q_R = numpy.zeros_like(q)
    
    for i in range(g.ngz - 1, g.nx + g.ngz + 1):
        sigma_up = q[i+1] - q[i]
        sigma_do = q[i] - q[i-1]
        sigma_bar = minmod(sigma_up, sigma_do)
        q_R[i] = q[i] - 0.5 * sigma_bar
        q_L[i+1] = q[i] + 0.5 * sigma_bar
        
    for i in range(g.ngz, g.nx + g.ngz + 1):
        f_upwind[i] = #! To be completed
        
    for i in range(g.ngz, g.nx + g.ngz):
        q_rhs[i] = #! To be completed

    return q_rhs



def rk2_minmod(nx, t_end, system,
               cfl_factor=0.9, initial_data=initial_data_exp_sine,
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
        if rp_case == 'rusanov':
            q1 = q + g.dt * minmod_rusanov_step(q, g, flux=system.flux)
        elif rp_case == 'upwind':
            q1 = q + g.dt * minmod_upwind_step(q, g, flux=system.flux)
        else:
            raise NotImplementedError
        q1 = boundary_conditions.boundaries(q1, g, form=bc_form)
        # Take second step
        if rp_case == 'rusanov':
            q = 0.5 * (q + q1 + g.dt * minmod_rusanov_step(q1, g, flux=system.flux))
        elif rp_case == 'upwind':
            q = 0.5 * (q + q1 + g.dt * minmod_upwind_step(q1, g, flux=system.flux))
        else:
            raise NotImplementedError
        q = boundary_conditions.boundaries(q, g, form=bc_form)
    # Done
    return g.x, q


if __name__ == "__main__":
    
    system = Advection()
    
    # Solve once at low resolution to see how it does
    x, q = rk2_minmod(20, 1.0, system, bc_form="periodic")
    x_exact = numpy.linspace(0, 1, 1000)
    q_exact = initial_data_exp_sine(x_exact)
    pyplot.plot(x_exact, q_exact, 'k-', label="Exact")
    pyplot.plot(x, q, 'bo', label=r"$N=20$")
    pyplot.xlabel(r"$x$")
    pyplot.ylabel(r"$q$")
    pyplot.legend()
    pyplot.title("RK2, Minmod, t=1")
    pyplot.show()
    
    # Solve for ten periods to see phase and dissipation errors
    x, q = rk2_minmod(20, 10.0, system, bc_form="periodic")
    x_exact = numpy.linspace(0, 1, 1000)
    q_exact = initial_data_exp_sine(x_exact)
    pyplot.plot(x_exact, q_exact, 'k-', label="Exact")
    pyplot.plot(x, q, 'bo', label=f"$N=20$")
    pyplot.xlabel(r"$x$")
    pyplot.ylabel(r"$q$")
    pyplot.legend()
    pyplot.title("RK2, Minmod, t=10")
    pyplot.show()
    
    # Solve with various grid spacings.
    # Check convergence
    
    Ns = 2**numpy.arange(4, 10)
    errors_1norm = numpy.zeros(Ns.shape)
    errors_2norm = numpy.zeros(Ns.shape)
    errors_infnorm = numpy.zeros(Ns.shape)
    dxs = 1.0 / Ns
    
    for i, N in enumerate(Ns):
        x, q = rk2_minmod(N, 1.0, system, bc_form="periodic")
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
    
    # Now do Burgers
    
    system = Burgers()
    
    # Solve for the smooth initial data at different times
    nx = 200
    Ts = [0.05, 0.1, 0.15, 0.2]
    mpl_formats=['b-', 'r-', 'g-', 'c-']
    for t_end, fmt in zip(Ts, mpl_formats):
        x, q = rk2_minmod(nx, t_end, system, initial_data=initial_data_exp_sine)
        pyplot.plot(x, q, fmt, label=fr"$t={t_end}$")
    pyplot.xlabel(r"$x$")
    pyplot.ylabel(r"$q$")
    pyplot.legend()
    pyplot.title("RK2 Minmod-Rusanov, Burgers, nx=200")
    pyplot.show()
    
    mpl_formats=['bo--', 'rx--', 'g^--', 'c*--']
    Ns = [10, 20, 40, 80]
    # Solve once at various resolutions to see how it does
    for nx, fmt in zip(Ns, mpl_formats):
        x, q = rk2_minmod(nx, 0.2, system, initial_data=initial_data_exp_sine)
        pyplot.plot(x, q, fmt, label=fr"$N={nx}$")
    pyplot.xlabel(r"$x$")
    pyplot.ylabel(r"$q$")
    pyplot.legend()
    pyplot.title("RK2 Minmod-Rusanov, Burgers, t=0.2")
    pyplot.show()
    
    # Repeat for the hat initial data to see the wave structure
    for nx, fmt in zip(Ns, mpl_formats):
        x, q = rk2_minmod(nx, 0.2, system, initial_data=initial_data_hat)
        pyplot.plot(x, q, fmt, label=fr"$N={nx}$")
    pyplot.xlabel(r"$x$")
    pyplot.ylabel(r"$q$")
    pyplot.legend()
    pyplot.title("RK2 Minmod-Rusanov, Burgers, t=0.2")
    pyplot.show()
    
    # Now back to advection using the upwind solver
    
    system = Advection()
    
    # Solve once at low resolution to see how it does
    x, q = rk2_minmod(20, 1.0, system, bc_form="periodic", rp_case="upwind")
    x_exact = numpy.linspace(0, 1, 1000)
    q_exact = initial_data_exp_sine(x_exact)
    pyplot.plot(x_exact, q_exact, 'k-', label="Exact")
    pyplot.plot(x, q, 'bo', label=r"$N=20$")
    pyplot.xlabel(r"$x$")
    pyplot.ylabel(r"$q$")
    pyplot.legend()
    pyplot.title("RK2, Minmod, upwind, t=1")
    pyplot.show()
    
    # Solve for ten periods to see phase and dissipation errors
    x, q = rk2_minmod(20, 10.0, system, bc_form="periodic", rp_case="upwind")
    x_exact = numpy.linspace(0, 1, 1000)
    q_exact = initial_data_exp_sine(x_exact)
    pyplot.plot(x_exact, q_exact, 'k-', label="Exact")
    pyplot.plot(x, q, 'bo', label=f"$N=20$")
    pyplot.xlabel(r"$x$")
    pyplot.ylabel(r"$q$")
    pyplot.legend()
    pyplot.title("RK2, Minmod, upwind, t=10")
    pyplot.show()
    
    # Solve with various grid spacings.
    # Check convergence
    
    Ns = 2**numpy.arange(4, 10)
    errors_1norm = numpy.zeros(Ns.shape)
    errors_2norm = numpy.zeros(Ns.shape)
    errors_infnorm = numpy.zeros(Ns.shape)
    dxs = 1.0 / Ns
    
    for i, N in enumerate(Ns):
        x, q = rk2_minmod(N, 1.0, system, bc_form="periodic", rp_case="upwind")
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
    