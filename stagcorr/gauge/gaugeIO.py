from operator import methodcaller
import numpy as np
from stagcorr.lattice.navigation import *

def load_gauge_field_vec(filename, vol, SUN=3):
    filename = "gauge_config/" + "_".join(map(str,vol)) + '/'  + filename
    with open(filename, "r") as f:
        lines = f.readlines()
        formatted_lines = list(map(methodcaller("split", "\t"), map(methodcaller("strip"), lines[3:])))

        physVol = np.prod(vol)
        Nlinks = len(vol)
        gVecVol = vol+(Nlinks,)+(SUN,SUN)
        gfieldArrList = []
        
        for i in range(physVol):
            gfieldArrList.append([])
            physCoord = lattice_coord_iso(position=i, vol=vol)
            for mu in range(Nlinks):
                gfieldArr = np.zeros((SUN, SUN), dtype=np.complex128)
                for j in range(SUN):
                    for k in range(SUN):
                        
                        
                        gfieldVecCoord = lattice_coord_iso_reverse(position=physCoord+[mu,j,k], vol=gVecVol)
                        gfieldArr[j,k] = float(formatted_lines[gfieldVecCoord][0]) + float(formatted_lines[gfieldVecCoord][1])*1j
                                
                gfieldArrList[i].append(gfieldArr)
        
    return gfieldArrList, vol

def save_gauge_field(filename, gfieldArrList, vol, SUN=3):
    gaugedir = "gauge_config/" + "_".join(map(str,vol)) + '/' 
    if not os.path.exists(gaugedir):
        os.makedirs(gaugedir)
    filename = gaugedir+ filename
    with open(filename, "w") as f:
        f.write("XXXX")
        f.write("\n")
        f.write(str(datetime.datetime.today()))
        f.write("\n")
        f.write("\t".join(map(str,vol)))
        f.write("\n")
        
        physVol = np.prod(vol)
        Nlinks = len(vol)
        gVecVol = vol+(Nlinks,)+(SUN,SUN)
        
        for i in range(physVol):
            physCoord = lattice_coord_iso(position=i, vol=vol)
            for mu in range(Nlinks):
                for j in range(SUN):
                    for k in range(SUN):
                        gfieldArr = gfieldArrList[i][mu]
                        x, y = gfieldArr[j,k].real , gfieldArr[j,k].imag
                        f.write("\t".join([str(x), str(y)]))
                        f.write("\n")
                        
    return "File saved at " + filename
        

def generate_unit_gauge_field_vec(vol, SUN=3):
    physVol = np.prod(vol)
    Nlinks = len(vol)
    gfieldArrList = []

    for i in range(physVol):
        gfieldArrList.append([])
        for mu in range(Nlinks):
            gfieldArrList[i].append(np.identity(3))
            
    return gfieldArrList