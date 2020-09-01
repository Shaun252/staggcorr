from stagcorr.stagFuncs import *
from stagcorr.lattice.navigation import *
from stagcorr.correlators.symmetricLinks import LinkSymList

from collections import defaultdict

########################################################################################################################################
#### 2pt function


def tieup2pt_fullProp(prop1, prop2, mom1, mom2, phase1, phase2, shift1, shift2, sym1, sym2, vol, SUN = 3):
    ##Assuming all spatial dimensions are equal
    
    N=vol[-1]
    matrix_dim = prop1.shape[0]
    mat_range = range(matrix_dim)
    mat_vol = vol + (SUN,)
    
    phys_dim  =  np.prod(vol)
    phys_range  =  range(phys_dim)
    
    temporal_dim =vol[0]
    temporal_range = range(temporal_dim)
    
    spatial_dim = phys_dim // temporal_dim
    spatial_dimxSUN = spatial_dim * SUN
    
    phase_func1 = phase_func_txyz(phase1)
    phase_func2 = phase_func_txyz(phase2)
    
    phase_mom_arr1 = np.empty(matrix_dim, dtype=np.complex128)
    phase_mom_arr2 = np.empty(matrix_dim, dtype=np.complex128)
    
    second_spatial_sum = np.zeros((temporal_dim, temporal_dim), dtype=np.complex128)

    coord_list = []
    
    shift1_indices_list = [i for i, x in enumerate(shift1) if x == 1]
    shift2_indices_list = [i for i, x in enumerate(shift2) if x == 1]
    
    Nlinks2 = len(shift2_indices_list)
    Nlinks1 = len(shift1_indices_list)
    
    Nlinks1 = len(shift1_indices_list)
    Nlinks2 = len(shift2_indices_list)
    
    coord_iso_dict1 = defaultdict(list)
    coord_iso_dict2 = defaultdict(list)
                            
    shift_bd_dict1 = defaultdict(list)
    shift_bd_dict2 = defaultdict(list)
                         
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

    for i in mat_range:
        
        coordi = lattice_coord_iso(i, vol=mat_vol)
        phys_coord = coordi[:-1]
        coord_list.append(phys_coord)
        
        stagPhase1 = phase_func1(phys_coord)
        stagPhase2 = phase_func2(phys_coord)
       
        momPhase1 = mom_character(mom1, phys_coord[1:], N)
        momPhase2 = mom_character(mom2, phys_coord[1:], N)
        
        phase_mom_arr1[i] = stagPhase1 * momPhase1
        phase_mom_arr2[i] = stagPhase2 * momPhase2
       
        
        for i, link_dirs in enumerate(sym1):
            
            shift = shift_list1[i]
            
            coordi_shift, bdphase = add_lattice_antiperiodic_txyz(phys_coord, shift, vol=vol)
            coordi_shift.append(coordi[-1])
            coord_iso_dict1[link_dirs].append(lattice_coord_iso_reverse(coordi_shift, vol=mat_vol))
            shift_bd_dict1[link_dirs].append(bdphase)
            

            
        for i, link_dirs in enumerate(sym2):
            
            shift = shift_list2[i]
            
            coordi_shift, bdphase = add_lattice_antiperiodic_txyz(phys_coord, shift, vol=vol)
            coordi_shift.append(coordi[-1])
            coord_iso_dict2[link_dirs].append(lattice_coord_iso_reverse(coordi_shift, vol=mat_vol))
            shift_bd_dict2[link_dirs].append(bdphase)
            
            
    if Nlinks2 == 0:
        Prop2PhaseMom1 = prop2 * phase_mom_arr1
    else: 
        Prop2PhaseMom1 = np.zeros((matrix_dim, matrix_dim), dtype=np.complex128)
        for link_dirs in sym2:
            Prop2PhaseMom1 += (shift_bd_dict2[link_dirs] * prop2[coord_iso_dict2[link_dirs],:].T).T
        Prop2PhaseMom1 /= no_terms2
        Prop2PhaseMom1 *= phase_mom_arr1
        
    del prop2
    del phase_mom_arr1
    
    
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
        lowerbnd = timej * spatial_dimxSUN
        upperbnd = (timej+1) * spatial_dimxSUN

        firstSpatialSumPhaseMom1 = Prop1PhaseMom2[:,lowerbnd:upperbnd] @ Prop2PhaseMom1[lowerbnd:upperbnd,:]
        
        for timek in temporal_range:
            lowerbnd = timek * spatial_dimxSUN
            upperbnd = (timek+1) * spatial_dimxSUN

            second_spatial_sum[timej,timek]= np.trace(firstSpatialSumPhaseMom1[lowerbnd:upperbnd,lowerbnd:upperbnd])    
       
    ### result is ntxntxnt array with [t3,t2,t1] where t3 is the tie up time, t2 is extended source time 
    ### and t1 is source time in the milc code. Also when using point sources any momentum is possible and the tie up 
    
    return second_spatial_sum