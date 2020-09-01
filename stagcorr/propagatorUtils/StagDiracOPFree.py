import numpy as np
import stagcorr.lattice.navigation as lattice
import stagcorr.stagFuncs as opT
from stagcorr.propagatorUtils.nnLists import generateAntiPeriodicNearestNeighbourLists, generatePeriodicNearestNeighbourLists
########################################################################################################################################
#### Generating the Staggered Operator

def staggered_operator(vol, mass, gField=None, naikeps=0, antiperiodic=True):
    ### Operator Params
    
    matrix_dim = np.prod(vol)
    matrix_range = range(matrix_dim)
    
    dim = len(vol)
    dim_range = range(dim)
    spatial_dim_range = range(1,dim)
    
    ### Initialise Operator Matrix
    stag_op = np.zeros((matrix_dim, matrix_dim), dtype=np.complex128)
    identityArr = 1
    
    ### Generate staggered eta phase and coordinate tables to link between physical lattice coords and stag op matrix indices
    coord_list = []
    eta_phase_list = []
    STIndexToMatrixIndexIsoList = []
    for i in matrix_range:
        coordi = lattice.lattice_coord_iso(i, vol=vol)
        coord_list.append(coordi)
        
        eta_phase_list.append(opT.eta_phase(coordi))
        
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
            
    ### Generate Operator SUNXSUN block at a time       
    for n in matrix_range:
        STCoordn = coord_list[n]
        eta_phase = eta_phase_list[n]

        ### NNeighbour lists
        nPlus1IndexList = []
        nMinus1IndexList = []
        nPlus2IndexList = []
        nMinus2IndexList = []
        nPlus3IndexList = []
        nMinus3IndexList = []
        
        bdryPhaseP1, bdryPhaseM1, bdryPhaseP3, bdryPhaseM3 = NNListFunc(STCoordn=STCoordn, unit_vec_list=unit_vec_list, unit_vec_list2=unit_vec_list2, unit_vec_list3=unit_vec_list3, nPlus1IndexList=nPlus1IndexList, nMinus1IndexList=nMinus1IndexList, nPlus2IndexList=nPlus2IndexList, nMinus2IndexList=nMinus2IndexList, nPlus3IndexList=nPlus3IndexList, nMinus3IndexList=nMinus3IndexList, vol=vol)
        
        ### Mass Delta
        stag_op[n, n] += mass * identityArr
        
        ### Time part of delta_m,n+1, separated to includes antiperiodic bd phases, includes naik terms
        m = nPlus1IndexList[0]
        
        stag_op[n, m] += bdryPhaseP1 * 0.5 * eta_phase[0] * ( 1 + (1+naikeps)/8)
        
        ### Time part of delta_m,n-1, separated to includes antiperiodic bd phases, includes naik terms
        m = nMinus1IndexList[0]
        
        stag_op[n, m] -= bdryPhaseM1 * 0.5 * eta_phase[0] * ( 1 + (1+naikeps)/8)
        
        ### Time part of delta_m,n+3, separated to includes antiperiodic bd phases, naik term
        m = nPlus3IndexList[0]

        stag_op[n, m] -= bdryPhaseP3 * (1+naikeps)*eta_phase[0]/48
        
        ### Time part of delta_m,n-3, separated to includes antiperiodic bd phases, naik term
        m = nMinus3IndexList[0]
        
        stag_op[n, m] += bdryPhaseM3 * (1+naikeps)*eta_phase[0]/48
        

        for mu in spatial_dim_range:
            ### Spatial part of delta_m,n+1, includes naik term
            m = nPlus1IndexList[mu]
            
            stag_op[n, m] += 0.5 * eta_phase[mu] * (1 + (1+naikeps)/8)
            
            ### Spatial part of delta_m,n-1, includes naik term
            m = nMinus1IndexList[mu]
            
            stag_op[n, m] -= 0.5 * eta_phase[mu] * (1 + (1+naikeps)/8)
            
            ### Spatial part of delta_m,n+3, naik term
            m = nPlus3IndexList[mu]
            
            stag_op[n, m] -= (1+naikeps)*eta_phase[mu]/48
            
            ### Spatial part of delta_m,n-3, naik term
            m = nMinus3IndexList[mu]
            
            stag_op[n, m] += (1+naikeps)*eta_phase[mu]/48
            
            
    return  2*stag_op