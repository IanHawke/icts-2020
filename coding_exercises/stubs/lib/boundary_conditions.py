"""
Simple boundary conditions
"""

def boundaries(q, g, form="periodic"):
    """
    Apply boundary conditions

    Parameters
    ----------
    q : array of float
        The solution, correct in the interior
    g : Grid
        Information about the grid
    form : string, optional
        The type of boundary condition to impose. The default is "periodic".

    Returns
    -------
    q : array of float
        The solution, correct everywhere
    """
    
    if q.ndim == 1:
        if form == "periodic":
            q[:g.ngz] = q[-2*g.ngz:-g.ngz]
            q[-g.ngz:] = q[g.ngz:2*g.ngz]
        elif form == "outflow":
            q[:g.ngz] = q[g.ngz]
            q[-g.ngz:] = q[-(g.ngz+1)]
        else:
            raise NotImplementedError
    elif q.ndim == 2:
        if form == "periodic":
            q[:g.ngz, :] = q[-2*g.ngz:-g.ngz, :]
            q[-g.ngz:, :] = q[g.ngz:2*g.ngz, :]
        elif form == "outflow":
            q[:g.ngz, :] = q[g.ngz, :]
            q[-g.ngz:, :] = q[-(g.ngz+1), :]
        elif form == "spherical_euler":
            q[:g.ngz, :] =  q[g.ngz:2*g.ngz, :][::-1]
            q[:g.ngz, 1] = -q[g.ngz:2*g.ngz, 1][::-1]  # Velocity flips sign
            q[-g.ngz:, :] = q[-(g.ngz+1), :]
        else:
            raise NotImplementedError
    else:
        raise NotImplementedError
    
    return q


def boundaries2d(q, g, form="periodic"):
    """
    Apply boundary conditions

    Parameters
    ----------
    q : array of float
        The solution, correct in the interior
    g : Grid
        Information about the grid
    form : string, optional
        The type of boundary condition to impose. The default is "periodic".

    Returns
    -------
    q : array of float
        The solution, correct everywhere
    """
    
    if q.ndim == 2:
        if form == "periodic":
            q[:, :g.ngz] = q[:, -2*g.ngz:-g.ngz]
            q[:, -g.ngz:] = q[:, g.ngz:2*g.ngz]
            q[:g.ngz, :] = q[-2*g.ngz:-g.ngz, :]
            q[-g.ngz:, :] = q[g.ngz:2*g.ngz, :]
        else:
            raise NotImplementedError
    else:
        raise NotImplementedError
    
    return q
