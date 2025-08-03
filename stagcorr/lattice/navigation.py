"""Lattice coordinate navigation utilities.

This module provides functions for converting between lattice coordinates
and linear matrix indices, essential for lattice QCD calculations on
hypercubic grids.
"""

import numpy as np

def lattice_coord_iso(position, vol):
    """Convert linear matrix index to lattice coordinates.
    
    Transforms a linear index (used for matrix addressing) back to 
    4D lattice coordinates [t, x, y, z].
    
    Parameters
    ----------
    position : int
        Linear matrix index (0 to prod(vol)-1)
    vol : tuple
        Lattice dimensions (T, X, Y, Z)
        
    Returns
    -------
    list
        Lattice coordinates [t, x, y, z] corresponding to the matrix index
        
    Examples
    --------
    >>> vol = (4, 4, 4, 4)
    >>> coords = lattice_coord_iso(85, vol)  
    >>> print(coords)  # [1, 1, 1, 1] for example
    
    Notes
    -----
    Uses row-major ordering with time as the slowest varying index.
    Inverse operation of lattice_coord_iso_reverse().
    """
    a = [0 for i in vol]
    i = 0
    while position:
        position,r = divmod(position, vol[::-1][i])
        i += 1
        a[-i] = r
    return a

def lattice_coord_iso_reverse(position, vol):
    """Convert lattice coordinates to linear matrix index.
    
    Transforms 4D lattice coordinates [t, x, y, z] to a linear index
    suitable for matrix addressing.
    
    Parameters
    ----------
    position : list or array_like
        Lattice coordinates [t, x, y, z]
    vol : tuple
        Lattice dimensions (T, X, Y, Z)
        
    Returns
    -------
    int
        Linear matrix index (0 to prod(vol)-1)
        
    Examples
    --------
    >>> vol = (4, 4, 4, 4)
    >>> index = lattice_coord_iso_reverse([1, 1, 1, 1], vol)
    >>> print(index)  # 85 for example
    
    Notes
    -----
    Uses row-major ordering with time as the slowest varying index.
    Inverse operation of lattice_coord_iso().
    """
    res = position[::-1][0]
    for i, coord in enumerate(position[::-1][1:]):
        res += coord * np.prod(vol[::-1][:i+1])
    return res


def add_lattice_periodic(coord1, coord2, vol):
    """Add lattice coordinates with periodic boundary conditions.
    
    Parameters
    ----------
    coord1, coord2 : array_like
        Lattice coordinates [t, x, y, z]
    vol : tuple
        Lattice dimensions (T, X, Y, Z)
        
    Returns
    -------
    list
        Sum of coordinates with periodic wrapping
    """
    lattice_coord = [(a + b) % v for a, b, v in zip(coord1, coord2, vol)]
    
    return lattice_coord

def subtract_lattice_periodic(coord1, coord2, vol):
    """Subtract lattice coordinates with periodic boundary conditions.
    
    Parameters
    ----------
    coord1, coord2 : array_like
        Lattice coordinates [t, x, y, z]
    vol : tuple
        Lattice dimensions (T, X, Y, Z)
        
    Returns
    -------
    list
        Difference of coordinates with periodic wrapping
    """
    lattice_coord = [(a - b) % v for a, b, v in zip(coord1, coord2, vol)]

    return lattice_coord

def add_lattice_antiperiodic_txyz(coord1, coord2, vol):
    """Add lattice coordinates with antiperiodic boundary conditions in time.
    
    Parameters
    ----------
    coord1, coord2 : array_like
        Lattice coordinates [t, x, y, z]
    vol : tuple
        Lattice dimensions (T, X, Y, Z)
        
    Returns
    -------
    tuple
        (lattice_coord, sign) where sign accounts for antiperiodic BC
    """
    lattice_coord = [(a + b) % v for a, b, v in zip(coord1, coord2, vol)]
    
    if (coord1[0] + coord2[0]) // vol[0] % 2 ==0:
        return lattice_coord, 1
    else:
        return lattice_coord, -1

def subtract_lattice_antiperiodic_txyz(coord1, coord2, vol):

    lattice_coord = [(a - b) % v for a, b, v in zip(coord1, coord2, vol)]
    
    if (coord1[0] - coord2[0]) // vol[0] % 2 ==0:
        return lattice_coord, 1
    else:
        return lattice_coord, -1
