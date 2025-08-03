"""Free field staggered fermion propagator computation.

This module handles the construction and inversion of staggered Dirac operators
to generate fermion propagators for free field lattice QCD calculations.
"""

import time
from stagcorr.propagatorUtils.propIOFree import saveProp, loadProp
from stagcorr.propagatorUtils.StagDiracOPFree import staggered_operator
from numpy.linalg import inv


def propagator(vol, mass, naikeps, antiperiodic=True, save=True):
    """Generate staggered fermion propagator for free field theory.
    
    Constructs the staggered Dirac operator and inverts it to obtain the 
    fermion propagator. Supports caching to disk for reuse with identical parameters.
    
    Parameters
    ----------
    vol : tuple
        Lattice dimensions (T, X, Y, Z)
    mass : float
        Bare quark mass in lattice units
    naikeps : float
        Naik improvement parameter epsilon
    antiperiodic : bool, default=True
        Use antiperiodic boundary conditions in time direction
    save : bool, default=True
        Save computed propagator to disk for future reuse
        
    Returns
    -------
    numpy.ndarray
        Fermion propagator matrix of shape (V, V) where V = prod(vol)
        Element [i,j] gives propagation amplitude from lattice site j to site i
        
    Notes
    -----
    - First attempts to load existing propagator from disk cache
    - If not found, constructs Dirac operator and inverts using numpy.linalg.inv
    - Prints timing information for operator construction and inversion
    - Uses 2×stag_op convention: propagator = (2×stag_op)^(-1)
    
    Examples
    --------
    >>> from stagcorr.propagatorUtils.propcoreFree import propagator
    >>> vol = (4, 4, 4, 4)
    >>> prop = propagator(vol=vol, mass=0.1, naikeps=0, antiperiodic=True)
    >>> print(prop.shape)  # (256, 256) for 4^4 lattice
    """
    gField = "Free"
    try:
        propagator = loadProp(vol=vol, mass=mass, naikeps=naikeps, antiperiodic=antiperiodic)
        #print("Propagator loaded")
    
    except:
        op_start = time.time()
        stag_op = staggered_operator(vol=vol, mass=mass, naikeps=naikeps, gField=gField, antiperiodic=antiperiodic)
        op_stop = time.time()
        print("Operator Computed, time: " +str(op_stop-op_start))

        prop_start = time.time()
        propagator = inv(stag_op)
        prop_stop = time.time()

        print("Operator Inverted, time: " +str(prop_stop-prop_start))
        
        if save == True:
            saveProp(prop=propagator, vol=vol, mass=mass, naikeps=naikeps, antiperiodic=antiperiodic)
    
    return propagator
