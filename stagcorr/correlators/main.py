"""
Main interface for staggered fermion correlation function calculations.

This module provides the primary entry point for computing n-point correlation 
functions in staggered fermion lattice QCD.
"""

import stagcorr.correlators.npoint as npoint
import stagcorr.correlators.npointFree as npointFree

def generate_npt(spinTasteMassNaikMomDag1, prop1=None, volume=(4,4,4,4), 
                 spinTasteMassNaikMomDag2=None, prop2=None, 
                 spinTasteMassNaikMomDag3=None, prop3=None, 
                 spinTasteMassNaikMomDag4=None, prop4=None, 
                 antiperiodic=True, gField=None):
    """
    Generate n-point staggered fermion correlation functions.
    
    Primary interface for computing 2-point, 3-point, and 4-point correlation 
    functions using staggered fermion operators. Automatically routes to free 
    field theory or interacting calculations based on gauge field presence.
    
    Parameters
    ----------
    spinTasteMassNaikMomDag1 : list
        First operator specification: [spin, taste, mass, naik_eps, momentum, dagger]
        - spin (str): Dirac structure ('G5', 'GX', 'GY', 'GZ', 'GT', etc.)
        - taste (str): Taste structure ('G5', 'G1', 'GX', 'GY', 'GZ', etc.)
        - mass (float): Bare quark mass in lattice units
        - naik_eps (float): Naik improvement parameter (typically 0)
        - momentum (list): [px, py, pz] in lattice units (2π/L)
        - dagger (bool): Whether the operator to be contracted has a dagger applied
    
    prop1 : numpy.ndarray, optional
        Pre-computed propagator for operator 1. If None, computed automatically.
    
    volume : tuple, default=(4,4,4,4)
        Lattice dimensions as (T, X, Y, Z)
    
    spinTasteMassNaikMomDag2 : list, optional
        Second operator specification (same format as first). 
        If None, creates 2-point function with operator 1.
    
    prop2 : numpy.ndarray, optional
        Pre-computed propagator for operator 2.
    
    spinTasteMassNaikMomDag3 : list, optional
        Third operator specification for 3-point or 4-point functions.
    
    prop3 : numpy.ndarray, optional
        Pre-computed propagator for operator 3.
    
    spinTasteMassNaikMomDag4 : list, optional  
        Fourth operator specification for 4-point functions.
    
    prop4 : numpy.ndarray, optional
        Pre-computed propagator for operator 4.
    
    antiperiodic : bool, default=True
        Use antiperiodic boundary conditions in time direction.
    
    gField : numpy.ndarray, optional
        Gauge field configuration. If None, uses free field theory.
    
    Returns
    -------
    numpy.ndarray
        Correlation function tensor:
        - 2-point: shape (nt, nt), indexed as C[t2, t1]
        - 3-point: shape (nt, nt, nt), indexed as C[t3, t2, t1]  
        - 4-point: shape (nt, nt, nt, nt), indexed as C[t4, t3, t2, t1]
        
        Complex values with real/imaginary parts depending on operator symmetries.
    
    Examples
    --------
    >>> import stagcorr.correlators.main as corr
    >>> 
    >>> # Pion 2-point function
    >>> vol = (6, 6, 6, 6)
    >>> Ct = corr.generate_npt(
    ...     spinTasteMassNaikMomDag1=["G5", "G5", 0.1, 0, [0,0,0], True],
    ...     spinTasteMassNaikMomDag2=["G5", "G5", 0.1, 0, [0,0,0], False],
    ...     volume=vol
    ... )
    >>> print(Ct.shape)  # (6, 6)
    >>>
    >>> # Two-pion to vector 3-point function  
    >>> Ct = corr.generate_npt(
    ...     spinTasteMassNaikMomDag1=["G5", "G5", 0.1, 0, [0,0,1], True],   # π⁺
    ...     spinTasteMassNaikMomDag2=["G5", "G5", 0.1, 0, [0,0,-1], True],  # π⁻
    ...     spinTasteMassNaikMomDag3=["GZ", "G1", 0.1, 0, [0,0,0], False],   # ρ
    ...     volume=vol
    ... )
    >>> print(Ct.shape)  # (6, 6, 6)
    
    Notes
    -----
    - Uses MILC convention: 2×stag_op = 1/2 × propagator
    - Results validated against MILC lattice QCD code  
    - Free field calculations use analytical methods for speed
    - Gauge field calculations use HISQ staggered action
    """

    if gField == None:

        corrArr = npointFree.npt(spinTasteMassNaikMomDag1=spinTasteMassNaikMomDag1, prop1 =prop1, volume=volume, spinTasteMassNaikMomDag2=spinTasteMassNaikMomDag2, prop2 = prop2, spinTasteMassNaikMomDag3=spinTasteMassNaikMomDag3, prop3 = prop3, spinTasteMassNaikMomDag4 = spinTasteMassNaikMomDag4, prop4 = prop4, antiperiodic=antiperiodic, gField = gField)
        
    else:
        corrArr =  npoint.npt(spinTasteMassNaikMomDag1=spinTasteMassNaikMomDag1, prop1 =prop1, volume=volume, spinTasteMassNaikMomDag2=spinTasteMassNaikMomDag2, prop2 = prop2, spinTasteMassNaikMomDag3=spinTasteMassNaikMomDag3, prop3 = prop3, spinTasteMassNaikMomDag4 = spinTasteMassNaikMomDag4, prop4 = prop4, antiperiodic=antiperiodic, gField = gField)
        
    return corrArr