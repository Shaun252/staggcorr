import stagcorr.stagFuncs as opT
from stagcorr.propagatorUtils.propcore import propagator
from stagcorr.correlators.twoPoints import tieup2pt_fullProp
from stagcorr.correlators.threePoints import tieup3pt_fullProp
from stagcorr.correlators.fourPoints import tieup4pt_fullProp

def npt(spinTasteMassNaikMomSymShift1, prop1 =None, volume=(4,4,4,4), spinTasteMassNaikMomSymShift2=None, prop2 =None, spinTasteMassNaikMomSymShift3=None, prop3 =None, spinTasteMassNaikMomSymShift4 = None, prop4 =None, antiperiodic=True, gField = None):
    
    phase1, shift1 = opT.phase_shift_operator(spin=spinTasteMassNaikMomSymShift1[0], taste=spinTasteMassNaikMomSymShift1[1])
    if type(prop1) == type(None):
        prop1 = propagator(vol=volume, mass=spinTasteMassNaikMomSymShift1[2], naikeps=spinTasteMassNaikMomSymShift1[3], antiperiodic=antiperiodic)
    mom1 = spinTasteMassNaikMomSymShift1[4]
    sym1 = spinTasteMassNaikMomSymShift1[5]
    
    if spinTasteMassNaikMomSymShift2 == None:
        phase2, shift2, prop2, mom2, one_sided2 = phase1, shift1, prop1, mom1, one_sided1
        
        correlator_arr = tieup2pt_fullProp(prop1=prop1, prop2=prop2, mom1=mom1, 
                                            mom2=mom2, phase1=phase1, phase2=phase2, 
                                            shift1=shift1, shift2=shift2,
                                            sym1=sym1, sym2=sym2,
                                            vol=volume)
        
    if spinTasteMassNaikMomSymShift2 != None and spinTasteMassNaikMomSymShift3==None:
    
        phase2, shift2 = opT.phase_shift_operator(spin=spinTasteMassNaikMomSymShift2[0], taste=spinTasteMassNaikMomSymShift2[1])
        if type(prop2) == type(None):
            prop2 = propagator(vol=volume, mass=spinTasteMassNaikMomSymShift2[2], naikeps=spinTasteMassNaikMomSymShift2[3], antiperiodic=antiperiodic)
        mom2 = spinTasteMassNaikMomSymShift2[4]
        sym2 = spinTasteMassNaikMomSymShift2[5]
        
        correlator_arr = tieup2pt_fullProp(prop1=prop1, prop2=prop2, mom1=mom1, 
                                            mom2=mom2, phase1=phase1, phase2=phase2, 
                                            shift1=shift1, shift2=shift2,
                                            sym1=sym1, sym2=sym2,
                                            vol=volume)
        
        
    if  spinTasteMassNaikMomSymShift3 != None and spinTasteMassNaikMomSymShift4 == None:
        
        phase2, shift2 = opT.phase_shift_operator(spin=spinTasteMassNaikMomSymShift2[0], taste=spinTasteMassNaikMomSymShift2[1])
        if type(prop2) == type(None):
            prop2 = propagator(vol=volume, mass=spinTasteMassNaikMomSymShift2[2], naikeps=spinTasteMassNaikMomSymShift2[3], antiperiodic=antiperiodic)
            
        mom2 = spinTasteMassNaikMomSymShift2[4]
        sym2 = spinTasteMassNaikMomSymShift2[5]
        
        phase3, shift3 = opT.phase_shift_operator(spin=spinTasteMassNaikMomSymShift3[0], taste=spinTasteMassNaikMomSymShift3[1])
        if type(prop3) == type(None):
            prop3 = propagator(vol=volume, mass=spinTasteMassNaikMomSymShift3[2], naikeps=spinTasteMassNaikMomSymShift3[3], antiperiodic=antiperiodic)
        mom3 = spinTasteMassNaikMomSymShift3[4]
        sym3 = spinTasteMassNaikMomSymShift3[5]
        
        
        correlator_arr = tieup3pt_fullProp(prop1=prop1, prop2=prop2, prop3=prop3, 
                                            mom1=mom1, mom2=mom2, mom3=mom3, 
                                            phase1=phase1, phase2=phase2, phase3=phase3,
                                            shift1=shift1, shift2=shift2, shift3=shift3,
                                            sym1=sym1, sym2=sym2, sym3=sym3,
                                            vol=volume)
        
    
    if  spinTasteMassNaikMomSymShift4 != None:
        
        phase2, shift2 = opT.phase_shift_operator(spin=spinTasteMassNaikMomSymShift2[0], taste=spinTasteMassNaikMomSymShift2[1])
        if type(prop2) == type(None):
            prop2 = propagator(vol=volume, mass=spinTasteMassNaikMomSymShift2[2], naikeps=spinTasteMassNaikMomSymShift2[3], antiperiodic=antiperiodic)
            
        mom2 = spinTasteMassNaikMomSymShift2[4]
        sym2 = spinTasteMassNaikMomSymShift2[5]

        
        phase3, shift3 = opT.phase_shift_operator(spin=spinTasteMassNaikMomSymShift3[0], taste=spinTasteMassNaikMomSymShift3[1])
        if type(prop3) == type(None):
            prop3 = propagator(vol=volume, mass=spinTasteMassNaikMomSymShift3[2], naikeps=spinTasteMassNaikMomSymShift3[3], antiperiodic=antiperiodic)
            
        mom3 = spinTasteMassNaikMomSymShift3[4]
        sym3 = spinTasteMassNaikMomSymShift3[5]

        
        phase4, shift4 = opT.phase_shift_operator(spin=spinTasteMassNaikMomSymShift4[0], taste=spinTasteMassNaikMomSymShift4[1])
        if type(prop4) == type(None):
            prop4 = propagator(vol=volume, mass=spinTasteMassNaikMomSymShift4[2], naikeps=spinTasteMassNaikMomSymShift4[3], antiperiodic=antiperiodic)
        mom4 = spinTasteMassNaikMomSymShift4[4]
        sym4 = spinTasteMassNaikMomSymShift4[5]

        
        correlator_arr = tieup4pt_fullProp(prop1=prop1, prop2=prop2, prop3=prop3, prop4=prop4,
                                            mom1=mom1, mom2=mom2, mom3=mom3, mom4=mom4, 
                                            phase1=phase1, phase2=phase2, phase3=phase3, phase4=phase4,
                                            shift1=shift1, shift2=shift2, shift3=shift3, shift4=shift4,
                                            sym1=sym1, sym2=sym2, sym3=sym3, sym4=sym4,
                                            vol=volume)
        
        
    return correlator_arr
