import numpy as np
import stagcorr.lattice.navigation as lattice
import stagcorr.stagFuncs as opT
from stagcorr.propagatorUtils.nnLists import generateAntiPeriodicNearestNeighbourLists, generatePeriodicNearestNeighbourLists
from stagcorr.gauge.HISQSmear import constructHISQGFields

########################################################################################################################################
#### Generating the Staggered Operator

def staggered_operator(vol, mass, gField=None, naikeps=0, antiperiodic=True):
    ### Operator Params
    SUN = 3
    physVol = np.prod(vol)
    physRange = range(physVol)
    
    matrix_dim = np.prod(vol+(SUN,))
    matrix_range = range(matrix_dim)
    
    dim = len(vol)
    dim_range = range(dim)
    spatial_dim_range = range(1,dim)
    
    ### Initialise Operator Matrix
    stag_op = np.zeros((matrix_dim, matrix_dim), dtype=np.complex128)
    identityArr = np.identity(SUN)
    
    ### Generate staggered eta phase and coordinate tables to link between physical lattice coords and stag op matrix indices
    coord_list = []
    eta_phase_list = []
    STIndexToMatrixIndexIsoList = []
    for i in physRange:
        coordi = lattice.lattice_coord_iso(i, vol=vol)
        coord_list.append(coordi)
        
        eta_phase_list.append(opT.eta_phase(coordi[:dim]))
        
        lowerMatrixIndex = i*SUN
        upperMatrixIndex = (i+1)*SUN
        
        STIndexToMatrixIndexIsoList.append((lowerMatrixIndex, upperMatrixIndex))
        
    ### Unit vector for delta functions in operator    
    unit_vec_list = []
    unit_vec_list2 = []
    unit_vec_list3 = []
    for mu in dim_range:
        unit_vec = [0 for i in dim_range]
        unit_vec[mu] = 1
        unit_vec_list.append(unit_vec[:])
        unit_vec[mu] = 2
        unit_vec_list2.append(unit_vec[:])
        unit_vec[mu] = 3
        unit_vec_list3.append(unit_vec[:])
    
    ### Functions to generate nearest neighbour coordinates in terms of matrix indices (not lattice coords)
    if antiperiodic == True:
        NNListFunc = generateAntiPeriodicNearestNeighbourLists
    else:
        NNListFunc = generatePeriodicNearestNeighbourLists
    
    ### Load Gauge
    XGaugeArrList, WGaugeArrList = constructHISQGFields(gField=gField, vol=vol)
            
    ### Generate Operator SUNXSUN block at a time       
    for n in physRange:
        STCoordn = coord_list[n]
        eta_phase = eta_phase_list[n]
        MatrixLowerN, MatrixUpperN = STIndexToMatrixIndexIsoList[n]

        ### NNeighbour lists
        nPlus1IndexList = []
        nMinus1IndexList = []
        nPlus2IndexList = []
        nMinus2IndexList = []
        nPlus3IndexList = []
        nMinus3IndexList = []
        
        bdryPhaseP1, bdryPhaseM1, bdryPhaseP3, bdryPhaseM3 = NNListFunc(STCoordn=STCoordn, unit_vec_list=unit_vec_list, unit_vec_list2=unit_vec_list2, unit_vec_list3=unit_vec_list3, nPlus1IndexList=nPlus1IndexList, nMinus1IndexList=nMinus1IndexList, nPlus2IndexList=nPlus2IndexList, nMinus2IndexList=nMinus2IndexList, nPlus3IndexList=nPlus3IndexList, nMinus3IndexList=nMinus3IndexList, vol=vol)
        
        ### Mass Delta
        stag_op[MatrixLowerN:MatrixUpperN, MatrixLowerN:MatrixUpperN] += mass * identityArr
        
        ### Time part of delta_m,n+1, separated to includes antiperiodic bd phases, includes naik terms
        m = nPlus1IndexList[0]
        MatrixLowerM,MatrixUpperM = STIndexToMatrixIndexIsoList[m]
        
        stag_op[MatrixLowerN:MatrixUpperN, MatrixLowerM:MatrixUpperM] += bdryPhaseP1 * 0.5 * eta_phase[0] * ( WGaugeArrList[n][0] + (1+naikeps)/8 * XGaugeArrList[n][0])
        
        ### Time part of delta_m,n-1, separated to includes antiperiodic bd phases, includes naik terms
        m = nMinus1IndexList[0]
        MatrixLowerM,MatrixUpperM = STIndexToMatrixIndexIsoList[m]
        
        stag_op[MatrixLowerN:MatrixUpperN, MatrixLowerM:MatrixUpperM] -= bdryPhaseM1 * 0.5 * eta_phase[0] * ( WGaugeArrList[m][0].conj().T + (1+naikeps)/8 *  XGaugeArrList[m][0].conj().T)
        
        ### Time part of delta_m,n+3, separated to includes antiperiodic bd phases, naik term
        m = nPlus3IndexList[0]
        MatrixLowerM,MatrixUpperM = STIndexToMatrixIndexIsoList[m]    
        
        stag_op[MatrixLowerN:MatrixUpperN, MatrixLowerM:MatrixUpperM] -= bdryPhaseP3 * (1+naikeps)*eta_phase[0]/48 * (XGaugeArrList[nPlus2IndexList[0]][0] @ XGaugeArrList[nPlus1IndexList[0]][0] @ XGaugeArrList[n][0])
        
        ### Time part of delta_m,n-3, separated to includes antiperiodic bd phases, naik term
        m = nMinus3IndexList[0]
        MatrixLowerM,MatrixUpperM = STIndexToMatrixIndexIsoList[m]
        
        stag_op[MatrixLowerN:MatrixUpperN, MatrixLowerM:MatrixUpperM] += bdryPhaseM3 * (1+naikeps)*eta_phase[0]/48 * (XGaugeArrList[nMinus1IndexList[0]][0] @ XGaugeArrList[nMinus2IndexList[0]][0] @ XGaugeArrList[m][0]).conj().T
        

        for mu in spatial_dim_range:
            ### Spatial part of delta_m,n+1, includes naik term
            m = nPlus1IndexList[mu]
            MatrixLowerM,MatrixUpperM = STIndexToMatrixIndexIsoList[m]
            
            stag_op[MatrixLowerN:MatrixUpperN, MatrixLowerM:MatrixUpperM] += 0.5 * eta_phase[mu] * (WGaugeArrList[n][mu] + (1+naikeps)/8 * XGaugeArrList[n][mu])
            
            ### Spatial part of delta_m,n-1, includes naik term
            m = nMinus1IndexList[mu]
            MatrixLowerM,MatrixUpperM = STIndexToMatrixIndexIsoList[m]
            
            stag_op[MatrixLowerN:MatrixUpperN, MatrixLowerM:MatrixUpperM] -= 0.5 * eta_phase[mu] * (WGaugeArrList[m][mu].conj().T + (1+naikeps)/8 *  XGaugeArrList[m][mu].conj().T)
            
            ### Spatial part of delta_m,n+3, naik term
            m = nPlus3IndexList[mu]
            MatrixLowerM,MatrixUpperM = STIndexToMatrixIndexIsoList[m]
            
            stag_op[MatrixLowerN:MatrixUpperN, MatrixLowerM:MatrixUpperM] -= (1+naikeps)*eta_phase[mu]/48 * (XGaugeArrList[nPlus2IndexList[mu]][mu] @ XGaugeArrList[nPlus1IndexList[mu]][mu] @ XGaugeArrList[n][mu])
            
            ### Spatial part of delta_m,n-3, naik term
            m = nMinus3IndexList[mu]
            MatrixLowerM,MatrixUpperM = STIndexToMatrixIndexIsoList[m]
            
            stag_op[MatrixLowerN:MatrixUpperN, MatrixLowerM:MatrixUpperM] += (1+naikeps)*eta_phase[mu]/48 * (XGaugeArrList[nMinus1IndexList[mu]][mu] @ XGaugeArrList[nMinus2IndexList[mu]][mu] @ XGaugeArrList[m][mu]).conj().T
            
            
    return  2*stag_op