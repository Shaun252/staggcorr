import stagcorr.correlators.npoint as npoint
import stagcorr.correlators.npointFree as npointFree

def generate_npt(spinTasteMassNaikMomSymShift1, prop1 =None, volume=(4,4,4,4), spinTasteMassNaikMomSymShift2=None, prop2 =None, spinTasteMassNaikMomSymShift3=None, prop3 =None, spinTasteMassNaikMomSymShift4 = None, prop4 =None, antiperiodic=True, gField = None):
    
    if gField == None:
        corrArr = npointFree.npt(spinTasteMassNaikMomSymShift1=spinTasteMassNaikMomSymShift1, prop1 =prop1, volume=volume, spinTasteMassNaikMomSymShift2=spinTasteMassNaikMomSymShift2, prop2 = prop2, spinTasteMassNaikMomSymShift3=spinTasteMassNaikMomSymShift3, prop3 = prop3, spinTasteMassNaikMomSymShift4 = spinTasteMassNaikMomSymShift4, prop4 = prop4, antiperiodic=antiperiodic, gField = gField)
        
    else:
        corrArr =  npoint.npt(spinTasteMassNaikMomSymShift1=spinTasteMassNaikMomSymShift1, prop1 =prop1, volume=volume, spinTasteMassNaikMomSymShift2=spinTasteMassNaikMomSymShift2, prop2 = prop2, spinTasteMassNaikMomSymShift3=spinTasteMassNaikMomSymShift3, prop3 = prop3, spinTasteMassNaikMomSymShift4 = spinTasteMassNaikMomSymShift4, prop4 = prop4, antiperiodic=antiperiodic, gField = gField)
        
    return corrArr