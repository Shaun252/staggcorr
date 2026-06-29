"""Free field n-point correlation function calculations.

This module handles the computation of correlation functions for free staggered 
fermions, routing to appropriate 2-point, 3-point, or 4-point calculations
based on the number of operators specified.
"""

import stagcorr.stagFuncs as opT
from stagcorr.propagatorUtils.propcoreFree import propagator
from stagcorr.correlators.twoPointsFreeFixed import tieup2pt_fullProp
from stagcorr.correlators.threePointsFreeFixed import tieup3pt_fullProp
from stagcorr.correlators.fourPointsFreeFixed import tieup4pt_fullProp

def npt(spinTasteMassNaikMomDag1, prop1 =None, volume=(4,4,4,4), spinTasteMassNaikMomDag2=None, prop2 =None, spinTasteMassNaikMomDag3=None, prop3 =None, spinTasteMassNaikMomDag4 = None, prop4 =None, antiperiodic=True, gField = None):
    """Compute n-point correlation functions for free staggered fermions.
    
    Internal routing function that determines whether to compute 2-point, 3-point,
    or 4-point correlation functions based on provided operator specifications.
    Handles propagator generation and operator phase calculations.
    
    Parameters
    ----------
    spinTasteMassNaikMomDag1 : list
        First operator: [spin, taste, mass, naik_eps, momentum, dagger]
    prop1 : numpy.ndarray, optional
        Pre-computed propagator for operator 1
    volume : tuple, default=(4,4,4,4)
        Lattice dimensions (T, X, Y, Z)
    spinTasteMassNaikMomDag2 : list, optional
        Second operator specification
    prop2 : numpy.ndarray, optional
        Pre-computed propagator for operator 2
    spinTasteMassNaikMomDag3 : list, optional
        Third operator specification  
    prop3 : numpy.ndarray, optional
        Pre-computed propagator for operator 3
    spinTasteMassNaikMomDag4 : list, optional
        Fourth operator specification
    prop4 : numpy.ndarray, optional
        Pre-computed propagator for operator 4
    antiperiodic : bool, default=True
        Use antiperiodic boundary conditions in time
    gField : unused
        Gauge field parameter (unused in free field theory)
        
    Returns
    -------
    numpy.ndarray
        Correlation function tensor with appropriate dimensions for n-point function
        
    Notes
    -----
    This is an internal function typically called via generate_npt() interface.
    """
    phase1, shift1 = opT.phase_shift_operator(spin=spinTasteMassNaikMomDag1[0], taste=spinTasteMassNaikMomDag1[1])
    if type(prop1) == type(None):
        prop1 = propagator(vol=volume, mass=spinTasteMassNaikMomDag1[2], naikeps=spinTasteMassNaikMomDag1[3], antiperiodic=antiperiodic)
    mom1 = spinTasteMassNaikMomDag1[4]
    dag1 = spinTasteMassNaikMomDag1[5]
    
    if spinTasteMassNaikMomDag2 == None:
        phase2, shift2, prop2, mom2 = phase1, shift1, prop1, mom1
        dag2 = not dag1
        
        correlator_arr = tieup2pt_fullProp(prop1=prop1, prop2=prop2, 
                                           mom1=mom1, mom2=mom2, 
                                           phase1=phase1, phase2=phase2, 
                                           shift1=shift1, shift2=shift2, 
                                           dag1=dag1, dag2=dag2, 
                                           vol=volume)
        
    if spinTasteMassNaikMomDag2 != None and spinTasteMassNaikMomDag3==None:
    
        phase2, shift2 = opT.phase_shift_operator(spin=spinTasteMassNaikMomDag2[0], taste=spinTasteMassNaikMomDag2[1])
        if type(prop2) == type(None):
            prop2 = propagator(vol=volume, mass=spinTasteMassNaikMomDag2[2], naikeps=spinTasteMassNaikMomDag2[3], antiperiodic=antiperiodic)
        mom2 = spinTasteMassNaikMomDag2[4]
        dag2 = spinTasteMassNaikMomDag2[5]
        
        correlator_arr = tieup2pt_fullProp(prop1=prop1, prop2=prop2, 
                                           mom1=mom1, mom2=mom2, 
                                           phase1=phase1, phase2=phase2, 
                                           shift1=shift1, shift2=shift2,
                                           dag1=dag1, dag2=dag2,
                                           vol=volume)
        
        
    if  spinTasteMassNaikMomDag3 != None and spinTasteMassNaikMomDag4 == None:
        
        phase2, shift2 = opT.phase_shift_operator(spin=spinTasteMassNaikMomDag2[0], taste=spinTasteMassNaikMomDag2[1])
        if type(prop2) == type(None):
            prop2 = propagator(vol=volume, mass=spinTasteMassNaikMomDag2[2], naikeps=spinTasteMassNaikMomDag2[3], antiperiodic=antiperiodic)
            
        mom2 = spinTasteMassNaikMomDag2[4]
        dag2 = spinTasteMassNaikMomDag2[5]
        
        phase3, shift3 = opT.phase_shift_operator(spin=spinTasteMassNaikMomDag3[0], taste=spinTasteMassNaikMomDag3[1])
        if type(prop3) == type(None):
            prop3 = propagator(vol=volume, mass=spinTasteMassNaikMomDag3[2], naikeps=spinTasteMassNaikMomDag3[3], antiperiodic=antiperiodic)
        mom3 = spinTasteMassNaikMomDag3[4]
        dag3 = spinTasteMassNaikMomDag3[5]
        
        
        correlator_arr = tieup3pt_fullProp(prop1=prop1, prop2=prop2, prop3=prop3, 
                                            mom1=mom1, mom2=mom2, mom3=mom3, 
                                            phase1=phase1, phase2=phase2, phase3=phase3,
                                            shift1=shift1, shift2=shift2, shift3=shift3,
                                            dag1=dag1, dag2=dag2, dag3=dag3,
                                            vol=volume)
        
    
    if  spinTasteMassNaikMomDag4 != None:
        
        phase2, shift2 = opT.phase_shift_operator(spin=spinTasteMassNaikMomDag2[0], taste=spinTasteMassNaikMomDag2[1])
        if type(prop2) == type(None):
            prop2 = propagator(vol=volume, mass=spinTasteMassNaikMomDag2[2], naikeps=spinTasteMassNaikMomDag2[3], antiperiodic=antiperiodic)
            
        mom2 = spinTasteMassNaikMomDag2[4]
        dag2 = spinTasteMassNaikMomDag2[5]

        
        phase3, shift3 = opT.phase_shift_operator(spin=spinTasteMassNaikMomDag3[0], taste=spinTasteMassNaikMomDag3[1])
        if type(prop3) == type(None):
            prop3 = propagator(vol=volume, mass=spinTasteMassNaikMomDag3[2], naikeps=spinTasteMassNaikMomDag3[3], antiperiodic=antiperiodic)
            
        mom3 = spinTasteMassNaikMomDag3[4]
        dag3 = spinTasteMassNaikMomDag3[5]

        
        phase4, shift4 = opT.phase_shift_operator(spin=spinTasteMassNaikMomDag4[0], taste=spinTasteMassNaikMomDag4[1])
        if type(prop4) == type(None):
            prop4 = propagator(vol=volume, mass=spinTasteMassNaikMomDag4[2], naikeps=spinTasteMassNaikMomDag4[3], antiperiodic=antiperiodic)
        mom4 = spinTasteMassNaikMomDag4[4]
        dag4 = spinTasteMassNaikMomDag4[5]

        
        correlator_arr = tieup4pt_fullProp(prop1=prop1, prop2=prop2, prop3=prop3, prop4=prop4,
                                            mom1=mom1, mom2=mom2, mom3=mom3, mom4=mom4, 
                                            phase1=phase1, phase2=phase2, phase3=phase3, phase4=phase4,
                                            shift1=shift1, shift2=shift2, shift3=shift3, shift4=shift4,
                                            dag1=dag1, dag2=dag2, dag3=dag3, dag4=dag4,
                                            vol=volume)
        
        
    return correlator_arr
