"""
TOV solution in polar areal coordinates
"""

import numpy
from scipy.integrate import solve_ivp


def initial_data_tov(r):
    """
    Generate a TOV, polytrope, p = 100 rho^2

    Parameters
    ----------
    r : array of float
        Coordinates

    Returns
    -------
    q : array of float
        Conserved variables at r
    alpha, a : array of float
        Metric variables at r
    """
    K = 100
    gamma = 2
    rho_c = 1.28e-3
    p_c = K * rho_c**gamma
    
    q = numpy.zeros((len(r), 3))
    alpha = numpy.zeros(len(r))
    a = numpy.zeros(len(r))
    
    def tov_rhs(r, tov_q):
        m, Phi, p = tov_q
        p = max(p, 0)
        rho0 = (p / K)**(1.0 / gamma)
        epsilon = K * rho0**(gamma - 1.0) / (gamma - 1.0)
        
        dqdr = numpy.zeros_like(tov_q)
        dqdr[0] = 4.0 * numpy.pi * r**2 * rho0 * (1.0 + epsilon)
        if r < 1e-10:
            dqdr[1] = 0
        else:
            dqdr[1] = (m + 4.0 * numpy.pi * r**3 * p) / (r * (r - 2.0 * m))
        dqdr[2] = -(rho0 * (1.0 + epsilon) + p) * dqdr[1]
        
        return dqdr
    
    sol = solve_ivp(tov_rhs, [0, 1000], [0, 0, p_c], dense_output=True,
                    max_step=0.1)
    solution = sol.sol(numpy.abs(r))
    p = solution[2, :]
    p[p < 0.0] = 0.0
    rho0 = (p / K)**(1.0 / gamma)
    epsilon = K * rho0**(gamma - 1.0) / (gamma - 1.0)
    q[:, 0] = rho0
    q[:, 1] = 0.0
    q[:, 2] = rho0 * epsilon
    
    m = solution[0, :]
    a = 1.0 / numpy.sqrt(1 - 2 * m / r)
    
    Phi_infinity = sol.y[1, -1]
    Phi = solution[1, :] - Phi_infinity
    alpha = numpy.exp(Phi)
    
    return q * a[:, numpy.newaxis] * r[:, numpy.newaxis]**2, alpha, a
