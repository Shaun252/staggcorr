from stagcorr.stagFuncs import *
from stagcorr.lattice.navigation import *
from stagcorr.correlators.symmetricLinks import LinkSymList
from stagcorr.propagatorUtils.propcore import apply_shift_bd_and_phase

from collections import defaultdict

########################################################################################################################################
#### 4pt function

def tieup4pt_fullProp(prop1, prop2, prop3, prop4, mom1, mom2, mom3, mom4, phase1, phase2, phase3, phase4, shift1, shift2, shift3, shift4, dag1, dag2, dag3, dag4, vol, SUN=3):
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
    phase_func4 = phase_func_txyz(phase4)
    
    phase_mom_arr1 = np.empty(matrix_dim, dtype=np.complex128)
    phase_mom_arr2 = np.empty(matrix_dim, dtype=np.complex128)
    phase_mom_arr3 = np.empty(matrix_dim, dtype=np.complex128)
    phase_mom_arr4 = np.empty(matrix_dim, dtype=np.complex128)
    
    fourth_spatial_sum = np.zeros((temporal_dim, temporal_dim, temporal_dim, temporal_dim), dtype=np.complex128)

    coord_list = []
    
    shift1_indices_list = [i for i, x in enumerate(shift1) if x == 1]
    shift2_indices_list = [i for i, x in enumerate(shift2) if x == 1]
    shift3_indices_list = [i for i, x in enumerate(shift3) if x == 1]
    shift4_indices_list = [i for i, x in enumerate(shift4) if x == 1]
    
    Nlinks4 = len(shift4_indices_list)
    Nlinks3 = len(shift3_indices_list)
    Nlinks2 = len(shift2_indices_list)
    Nlinks1 = len(shift1_indices_list)
    
    coord_iso_dict1 = defaultdict(list)
    coord_iso_dict2 = defaultdict(list)
    coord_iso_dict3 = defaultdict(list)
    coord_iso_dict4 = defaultdict(list)
                            
    shift_bd_dict1 = defaultdict(list)
    shift_bd_dict2 = defaultdict(list)
    shift_bd_dict3 = defaultdict(list)
    shift_bd_dict4 = defaultdict(list)
                         
    sym1 = LinkSymList[Nlinks1]
    no_terms1 = 2 ** Nlinks1
 
    sym2 = LinkSymList[Nlinks2]
    no_terms2 = 2 ** Nlinks2
        
    sym3 = LinkSymList[Nlinks3]
    no_terms3 = 2 ** Nlinks3

    sym4 = LinkSymList[Nlinks4]
    no_terms4 = 2 ** Nlinks4
    
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
    
    shift_list4 = []
    for link_dirs in sym4:
        shift = np.zeros((4), dtype=int)
        for i, k in enumerate(shift4_indices_list):
            shift[k] = int(link_dirs[i])
        
        shift_list4.append(shift)

    for i in mat_range:
        
        coordi = lattice_coord_iso(i, vol=vol)
        coord_list.append(coordi)
        
        stagPhase1 = phase_func1(coordi)
        stagPhase2 = phase_func2(coordi)
        stagPhase3 = phase_func3(coordi)
        stagPhase4 = phase_func4(coordi)
       
        momPhase1 = mom_character(mom1, coordi[1:], N, dag1)
        momPhase2 = mom_character(mom2, coordi[1:], N, dag2)
        momPhase3 = mom_character(mom3, coordi[1:], N, dag3)
        momPhase4 = mom_character(mom4, coordi[1:], N, dag4)
        
        phase_mom_arr1[i] = stagPhase1 * momPhase1
        phase_mom_arr2[i] = stagPhase2 * momPhase2
        phase_mom_arr3[i] = stagPhase3 * momPhase3
        phase_mom_arr4[i] = stagPhase4 * momPhase4
        
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
            
        for i, link_dirs in enumerate(sym4):
            
            shift = shift_list4[i]
            
            coordi_shift, bdphase = add_lattice_antiperiodic_txyz(coordi, shift, vol=vol)
            coord_iso_dict4[link_dirs].append(lattice_coord_iso_reverse(coordi_shift, vol=vol))
            shift_bd_dict4[link_dirs].append(bdphase)
        
    
    prop1 = apply_shift_bd_and_phase(prop=prop1, coord_iso_dict_l=coord_iso_dict2, coord_iso_dict_r=coord_iso_dict1,
                                     bd_dict_l = shift_bd_dict2, bd_dict_r = shift_bd_dict1,
                                     no_terms_l = no_terms2, no_terms_r = no_terms1,
                                     dag_l = dag2, dag_r = dag1,
                                     phase_l = phase_mom_arr2, phase_r = phase_mom_arr1,
                                     matrix_dim=matrix_dim)
    
    
    prop2 = apply_shift_bd_and_phase(prop=prop2, coord_iso_dict_l=coord_iso_dict3, coord_iso_dict_r=coord_iso_dict2,
                                     bd_dict_l = shift_bd_dict3, bd_dict_r = shift_bd_dict2,
                                     no_terms_l = no_terms3, no_terms_r = no_terms2,
                                     dag_l = dag3, dag_r = dag2,
                                     phase_l = phase_mom_arr3, phase_r = phase_mom_arr2,
                                     matrix_dim=matrix_dim)
    
    prop3 = apply_shift_bd_and_phase(prop=prop3, coord_iso_dict_l=coord_iso_dict4, coord_iso_dict_r=coord_iso_dict3,
                                     bd_dict_l = shift_bd_dict4, bd_dict_r = shift_bd_dict3,
                                     no_terms_l = no_terms4, no_terms_r = no_terms3,
                                     dag_l = dag4, dag_r = dag3,
                                     phase_l = phase_mom_arr4, phase_r = phase_mom_arr3,
                                     matrix_dim=matrix_dim)
    
    prop4 = apply_shift_bd_and_phase(prop=prop4, coord_iso_dict_l=coord_iso_dict1, coord_iso_dict_r=coord_iso_dict4,
                                     bd_dict_l = shift_bd_dict1, bd_dict_r = shift_bd_dict4,
                                     no_terms_l = no_terms1, no_terms_r = no_terms4,
                                     dag_l = dag1, dag_r = dag4,
                                     phase_l = phase_mom_arr1, phase_r = phase_mom_arr4,
                                     matrix_dim=matrix_dim)
    

    for timej in temporal_range:
        lowerbnd = timej * spatial_dim
        upperbnd = (timej+1) * spatial_dim

        firstSpatialSum = prop4[:,lowerbnd:upperbnd] @ prop3[lowerbnd:upperbnd,:]
        
        for timek in temporal_range:
            lowerbnd = timek * spatial_dim
            upperbnd = (timek+1) * spatial_dim
            
            secondSpatialSum = firstSpatialSum[:,lowerbnd:upperbnd] @ prop2[lowerbnd:upperbnd,:]

            for timel in temporal_range:
                lowerbnd = timel * spatial_dim
                upperbnd = (timel+1) * spatial_dim
                    
                thirdSpatialSum = secondSpatialSum[:,lowerbnd:upperbnd] @ prop1[lowerbnd:upperbnd,:] 
                
                for timem in temporal_range:
                    lowerbnd = timem * spatial_dim
                    upperbnd = (timem+1) * spatial_dim
                    
                    fourth_spatial_sum[timej,timek,timel,timem]= np.trace(thirdSpatialSum[lowerbnd:upperbnd,lowerbnd:upperbnd])    
       
    ### result is ntxntxnt array with [t3,t2,t1] where t3 is the tie up time, t2 is extended source time 
    ### and t1 is source time in the milc code. Also when using point sources any momentum is possible and the tie up 
    
    return fourth_spatial_sum*3
