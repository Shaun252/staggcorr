import time
from stagcorr.propagatorUtils.propIO import saveProp, loadProp
from stagcorr.propagatorUtils.StagDiracOP import staggered_operator
from numpy.linalg import inv


def propagator(vol, mass, naikeps, gField=None, antiperiodic=True, save = True):
    
    try:
        propagator = loadProp(vol=vol, mass=mass, naikeps=naikeps, gField=gField, antiperiodic=antiperiodic)
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
            saveProp(prop=propagator, vol=vol, mass=mass, naikeps=naikeps, gField=gField, antiperiodic=antiperiodic)

    return propagator
