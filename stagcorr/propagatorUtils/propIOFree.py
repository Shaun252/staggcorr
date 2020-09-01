import numpy as np
import os

def saveProp(prop, vol, mass, naikeps, antiperiodic):
    volstring = "".join(map(str,vol))
    mstring = "m"+str(float(mass))
    naikstring = "naik" + str(float(naikeps))

    gField = "Free" 
    if antiperiodic == True:
        pstring = "antiperiodic"
    else:
        pstring = "periodic"
        
    propdir = "./props/" + str(volstring) + "/" + gField + "/"
    
    propname= "_".join([mstring, naikstring, pstring])
    
    saveloc = propdir + propname
    
    if not os.path.exists(propdir):
        os.makedirs(propdir)
    
    np.save(saveloc, prop)
    
    return "Prop saved"

def loadProp(vol, mass, naikeps, antiperiodic):
    volstring = "".join(map(str,vol))
    mstring = "m"+str(float(mass))
    naikstring = "naik" + str(float(naikeps))
    gField = "Free"  
    if antiperiodic == True:
        pstring = "antiperiodic"
    else:
        pstring = "periodic"
        
    propdir = "./props/" + str(volstring) + "/" + gField + "/"
    
    propname= "_".join([mstring, naikstring, pstring])
    
    loadloc = propdir + propname +".npy"
    
    prop = np.load(loadloc)
    
    return prop