"""
Exercise 2: discontinuities with the upwind scheme
"""

import numpy
from matplotlib import pyplot

from lib import errors
from ex01_upwind import upwind


def initial_data_hat(x):
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
    
    q = numpy.where(numpy.abs(x-0.5) < 0.25,
                    numpy.ones_like(x), numpy.zeros_like(x))
    
    return q


if __name__ == "__main__":
    
    # Solve once at low resolution to see how it does
    x, q = upwind(20, 1.0, initial_data=initial_data_hat)
    x_exact = numpy.linspace(0, 1, 1000)
    q_exact = initial_data_hat(x_exact)
    pyplot.plot(x_exact, q_exact, 'k-', label="Exact")
    pyplot.plot(x, q, 'bo', label=f"$N=20$")
    pyplot.xlabel(r"$x$")
    pyplot.ylabel(r"$q$")
    pyplot.legend()
    pyplot.title("Upwind, t=1")
    pyplot.show()
    
    # Solve for ten periods to see phase and dissipation errors
    x, q = upwind(20, 10.0, initial_data=initial_data_hat)
    x_exact = numpy.linspace(0, 1, 1000)
    q_exact = initial_data_hat(x_exact)
    pyplot.plot(x_exact, q_exact, 'k-', label="Exact")
    pyplot.plot(x, q, 'bo', label=r"$N=20$")
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
        x, q = upwind(N, 1.0, initial_data=initial_data_hat)
        q_exact = initial_data_hat(x)
        errors_1norm[i] = errors.error_norm(q, q_exact, dxs[i], p=1)
        errors_2norm[i] = errors.error_norm(q, q_exact, dxs[i], p=2)
        errors_infnorm[i] = errors.error_norm(q, q_exact, dxs[i], p='inf')
    
    pyplot.loglog(dxs, errors_1norm, 'x',
                  label=r"|$q_{exact} - q_{numerical}|_1$")
    pyplot.loglog(dxs, errors_2norm, 'o',
                  label=r"|$q_{exact} - q_{numerical}|_2$")
    pyplot.loglog(dxs, errors_infnorm, '^',
                  label=r"|$q_{exact} - q_{numerical}|_{\infty}$")
    pyplot.loglog(dxs, errors_1norm[-1] * (dxs / dxs[-1])**(1/2),
                  label=r"$\propto (\Delta x)^{1/2}$")
    pyplot.xlabel(r"$\Delta x$")
    pyplot.ylabel("Errors")
    pyplot.legend()
    pyplot.show()
    
    # Compare a single run after one period with different CFLs
    cfls = [0.1, 0.5, 0.9, 0.99]
    mpl_formats=['bo', 'rx', 'g^', 'c*']
    for cfl, fmt in zip(cfls, mpl_formats):
        x, q = upwind(20, 1.0, cfl_factor=cfl,
                      initial_data=initial_data_hat)
        pyplot.plot(x, q, fmt,
                    label=f"CFL={cfl}")
    x_exact = numpy.linspace(0, 1, 1000)
    q_exact = initial_data_hat(x_exact)
    pyplot.plot(x_exact, q_exact, 'k-', label="Exact")
    pyplot.xlabel(r"$x$")
    pyplot.ylabel(r"$q$")
    pyplot.legend()
    pyplot.title("Upwind, t=1")
    pyplot.show()