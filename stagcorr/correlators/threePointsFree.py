from stagcorr.stagFuncs import *
from stagcorr.lattice.navigation import *
from stagcorr.correlators.symmetricLinks import LinkSymList

from collections import defaultdict

########################################################################################################################################
#### 3pt function

def tieup3pt_fullProp(prop1, prop2, prop3, mom1, mom2, mom3, phase1, phase2, phase3, shift1, shift2, shift3, sym1, sym2, sym3, vol, SUN=3):
    ##Assuming all spatial dimensions are equal
    
    N=vol[-1]
    matrix_dim = prop1.shape[0]
    mat_range = range(matrix_dim)
    
    temporal_dim =vol[0]
    temporal_range = range(temporal_dim)
    
    spatial_dim = matrix_dim // temporal_dim
    
    phase_func1 = phase_func_txyz(phase1)
    phase_func2 = phase_func_txyz(phase2)
    phase_func3 = phase_func_txyz(phase3)
    
    phase_mom_arr1 = np.empty(matrix_dim, dtype=np.complex128)
    phase_mom_arr2 = np.empty(matrix_dim, dtype=np.complex128)
    phase_mom_arr3 = np.empty(matrix_dim, dtype=np.complex128)
    
    third_spatial_sum = np.zeros((temporal_dim, temporal_dim, temporal_dim), dtype=np.complex128)

    coord_list = []
    
    shift1_indices_list = [i for i, x in enumerate(shift1) if x == 1]
    shift2_indices_list = [i for i, x in enumerate(shift2) if x == 1]
    shift3_indices_list = [i for i, x in enumerate(shift3) if x == 1]
    
    Nlinks3 = len(shift3_indices_list)
    Nlinks2 = len(shift2_indices_list)
    Nlinks1 = len(shift1_indices_list)
    
    coord_iso_dict1 = defaultdict(list)
    coord_iso_dict2 = defaultdict(list)
    coord_iso_dict3 = defaultdict(list)
                            
    shift_bd_dict1 = defaultdict(list)
    shift_bd_dict2 = defaultdict(list)
    shift_bd_dict3 = defaultdict(list)
                         
    if sym1 == 0:
        sym1 = LinkSymList[Nlinks1-1]
        no_terms1 = 2 ** Nlinks1
    else:
        no_terms1 = len(sym1)
        
    if sym2 == 0:
        sym2 = LinkSymList[Nlinks2-1]
        no_terms2 = 2 ** Nlinks2
    else:
        no_terms2 = len(sym2)
        
    if sym3 == 0:
        sym3 = LinkSymList[Nlinks3-1]
        no_terms3 = 2 ** Nlinks3
    else:
        no_terms3 = len(sym3)
    
    shift_list1 = []
    for link_dirs in sym1:
        shift = np.zeros((4), dtype=int)
        for i, k in enumerate(shift1_indices_list):
            shift[k] = int(link_dirs[i])
        
        shift_list1.append(shift)
    
    shift_list2 = []
    for link_dirs in sym2:
        shift = np.zeros((4), dtype=int)
        for i, k in enumerate(shift2_indices_list):
            shift[k] = int(link_dirs[i])
        
        shift_list2.append(shift)

    shift_list3 = []
    for link_dirs in sym3:
        shift = np.zeros((4), dtype=int)
        for i, k in enumerate(shift3_indices_list):
            shift[k] = int(link_dirs[i])
        
        shift_list3.append(shift)

    for i in mat_range:
        
        coordi = lattice_coord_iso(i, vol=vol)
        coord_list.append(coordi)
        
        stagPhase1 = phase_func1(coordi)
        stagPhase2 = phase_func2(coordi)
        stagPhase3 = phase_func3(coordi)
       
        momPhase1 = mom_character(mom1, coordi[1:], N)
        momPhase2 = mom_character(mom2, coordi[1:], N)
        momPhase3 = mom_character(mom3, coordi[1:], N)
        
        phase_mom_arr1[i] = stagPhase1 * momPhase1
        phase_mom_arr2[i] = stagPhase2 * momPhase2
        phase_mom_arr3[i] = stagPhase3 * momPhase3
        
        for i, link_dirs in enumerate(sym1):
            
            shift = shift_list1[i]
            
            coordi_shift, bdphase = add_lattice_antiperiodic_txyz(coordi, shift, vol=vol)            
            coord_iso_dict1[link_dirs].append(lattice_coord_iso_reverse(coordi_shift, vol=vol))
            shift_bd_dict1[link_dirs].append(bdphase)

            
        for i, link_dirs in enumerate(sym2):
            
            shift = shift_list2[i]
            
            coordi_shift, bdphase = add_lattice_antiperiodic_txyz(coordi, shift, vol=vol)
            coord_iso_dict2[link_dirs].append(lattice_coord_iso_reverse(coordi_shift, vol=vol))
            shift_bd_dict2[link_dirs].append(bdphase)
            
        for i, link_dirs in enumerate(sym3):
            
            shift = shift_list3[i]
            
            coordi_shift, bdphase = add_lattice_antiperiodic_txyz(coordi, shift, vol=vol)
            coord_iso_dict3[link_dirs].append(lattice_coord_iso_reverse(coordi_shift, vol=vol))
            shift_bd_dict3[link_dirs].append(bdphase)

###Idea here is that the 3 point is D(n1+d1|n2)e^{ip2n2} D(n2+d2|n3)e^{ip3n3} D(n3+d3|n1)e^{ip1n1} where prop1 = D(n1+d1|n2). This means we are free to multiply the whole propagator elementwise by the momentum (and staggered phase) corresponding to the second index but not first as this has a delta shift di


    if Nlinks3 == 0:
        Prop3PhaseMom1 = prop3 * phase_mom_arr1
    else: 
        Prop3PhaseMom1 = np.zeros((matrix_dim, matrix_dim), dtype=np.complex128)
        for link_dirs in sym3:
            Prop3PhaseMom1 += (shift_bd_dict3[link_dirs] * prop3[coord_iso_dict3[link_dirs],:].T).T
        Prop3PhaseMom1 /= no_terms3
        Prop3PhaseMom1 *= phase_mom_arr1
        
    del prop3
    del phase_mom_arr1
    
    if Nlinks2 == 0:
        Prop2PhaseMom3 = prop2 * phase_mom_arr3
    else: 
        Prop2PhaseMom3 = np.zeros((matrix_dim, matrix_dim), dtype=np.complex128)
        for link_dirs in sym2:
            Prop2PhaseMom3 += (shift_bd_dict2[link_dirs] * prop2[coord_iso_dict2[link_dirs],:].T).T
        Prop2PhaseMom3 /= no_terms2
        Prop2PhaseMom3 *= phase_mom_arr3
        
    del prop2
    del phase_mom_arr3
    
    
    if Nlinks1 ==0:
        Prop1PhaseMom2 = prop1 * phase_mom_arr2     
    else:
        Prop1PhaseMom2 = np.zeros((matrix_dim, matrix_dim), dtype=np.complex128)
        for link_dirs in sym1: 
            Prop1PhaseMom2 += (shift_bd_dict1[link_dirs] * prop1[coord_iso_dict1[link_dirs],:].T).T
        Prop1PhaseMom2 /= no_terms1  
        Prop1PhaseMom2 *= phase_mom_arr2
        
    del prop1
    del phase_mom_arr2

    for timej in temporal_range:
        lowerbnd = timej * spatial_dim
        upperbnd = (timej+1) * spatial_dim

        firstSpatialSumPhaseMom1 = Prop2PhaseMom3[:,lowerbnd:upperbnd] @ Prop3PhaseMom1[lowerbnd:upperbnd,:]
        
        for timek in temporal_range:
            lowerbnd = timek * spatial_dim
            upperbnd = (timek+1) * spatial_dim

            secondSpatialSumPhaseMom1 = Prop1PhaseMom2[:,lowerbnd:upperbnd] @ firstSpatialSumPhaseMom1[lowerbnd:upperbnd,:]

            for timel in temporal_range:
                lowerbnd = timel * spatial_dim
                upperbnd = (timel+1) * spatial_dim
                    
                third_spatial_sum[timej,timek,timel]= np.trace(secondSpatialSumPhaseMom1[lowerbnd:upperbnd,lowerbnd:upperbnd])    
       
    ### result is ntxntxnt array with [t3,t2,t1] where t3 is the tie up time, t2 is extended source time 
    ### and t1 is source time in the milc code. Also when using point sources any momentum is possible and the tie up 
    
    return third_spatial_sum *3
