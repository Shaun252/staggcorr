from math import sqrt
import numpy as np

def reuniterise(gMatrix):
    
    g00r = gMatrix[0,0].real
    g00i = gMatrix[0,0].imag
    g01r = gMatrix[0,1].real
    g01i = gMatrix[0,1].imag
    g02r = gMatrix[0,2].real
    g02i = gMatrix[0,2].imag
    
    row1Norm = sqrt(g00r**2 + g00i**2 + g01r**2 + g01i**2 + g02r**2 + g02i**2)
    
    gMatrix[0,:] = gMatrix[0,:] / row1Norm
    
    row1ProjectRow2 = np.dot(gMatrix[1,:], gMatrix[0,:].conj())
    
    gMatrix[1,:] = gMatrix[1,:] - gMatrix[0,:] * row1ProjectRow2    
    
    g10r = gMatrix[1,0].real
    g10i = gMatrix[1,0].imag
    g11r = gMatrix[1,1].real
    g11i = gMatrix[1,1].imag
    g12r = gMatrix[1,2].real
    g12i = gMatrix[1,2].imag
    
    row2Norm = sqrt(g10r**2 + g10i**2 + g11r**2 + g11i**2 + g12r**2 + g12i**2)
    
    gMatrix[1,:] = gMatrix[1,:] / row2Norm
    
    gMatrix[2,0] = gMatrix[0,1]*gMatrix[1,2] - gMatrix[0,2]*gMatrix[1,1]
    gMatrix[2,1] = gMatrix[0,2]*gMatrix[1,0] - gMatrix[0,0]*gMatrix[1,2]
    gMatrix[2,2] = gMatrix[0,0]*gMatrix[1,1] - gMatrix[0,1]*gMatrix[1,0]
    gMatrix[2,:] = gMatrix[2,:].conj()
    
    return gMatrix


def matrixDistance(gMatrixOld, gMatrixNew):
    diff = gMatrixOld- gMatrixNew
    sqrtMat = diff * diff.conj()
    dist = sqrtMat.sum()
    return dist

def checkUnitarity(gFieldArr):
    for i, row1 in enumerate(gFieldArr):
        for j, row2 in enumerate(gFieldArr):
            print(i,j)
            print(np.dot(row1, row2.conj()))