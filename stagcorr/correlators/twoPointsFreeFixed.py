"""Two-point correlation function calculations for free staggered fermions.

This module implements the core Wick contraction algorithms for computing
2-point correlation functions from staggered fermion propagators.
"""

from stagcorr.stagFuncs import *
from stagcorr.lattice.navigation import *
from stagcorr.correlators.symmetricLinks import LinkSymList
from stagcorr.propagatorUtils.propcore import apply_shift_bd_and_phase

from collections import defaultdict

########################################################################################################################################
#### 2pt function


def tieup2pt_fullProp(prop1, prop2, mom1, mom2, phase1, phase2, shift1, shift2, dag1, dag2, vol, SUN = 3):
    """Compute 2-point correlation function from propagators.
    
    Performs Wick contraction to compute 2-point correlation function
    C(t2,t1) = <O2(t2) O1(t1)> using pre-computed propagators.
    
    Parameters
    ----------
    prop1, prop2 : numpy.ndarray
        Fermion propagators for operators 1 and 2
    mom1, mom2 : list
        3-momentum vectors [px, py, pz] for operators 1 and 2
    phase1, phase2 : sympy expression
        Staggered phase expressions for operators 1 and 2
    shift1, shift2 : numpy.ndarray
        Coordinate shifts for operators 1 and 2
    dag1, dag2 : int
        Symmetric shift parameters for operators 1 and 2
    vol : tuple
        Lattice dimensions (T, X, Y, Z)
    SUN : int, default=3
        SU(N) color group dimension (typically 3 for QCD)
        
    Returns
    -------
    numpy.ndarray
        2-point correlation function matrix of shape (nt, nt)
        where C[t2, t1] = <O2(t2) O1(t1)>
        
    Notes
    -----
    - Implements full Wick contraction with momentum projection
    - Applies staggered fermion phases and coordinate shifts
    - Uses MILC normalization conventions
    """
    ##Assuming all spatial dimensions are equal

    N=vol[-1]
    matrix_dim = prop1.shape[0]
    mat_range = range(matrix_dim)

    temporal_dim =vol[0]
    temporal_range = range(temporal_dim)
    
    spatial_dim = matrix_dim // temporal_dim
    
    phase_func1 = phase_func_txyz(phase1)
    phase_func2 = phase_func_txyz(phase2)
    
    phase_mom_arr1 = np.empty(matrix_dim, dtype=np.complex128)
    phase_mom_arr2 = np.empty(matrix_dim, dtype=np.complex128)
    
    second_spatial_sum = np.zeros((temporal_dim, temporal_dim), dtype=np.complex128)

    coord_list = []
    
    shift1_indices_list = [i for i, x in enumerate(shift1) if x == 1]
    shift2_indices_list = [i for i, x in enumerate(shift2) if x == 1]
           
    Nlinks1 = len(shift1_indices_list)
    Nlinks2 = len(shift2_indices_list)
    
    coord_iso_dict1 = defaultdict(list)
    coord_iso_dict2 = defaultdict(list)
                            
    shift_bd_dict1 = defaultdict(list)
    shift_bd_dict2 = defaultdict(list)
                         
    sym1 = LinkSymList[Nlinks1]
    no_terms1 = 2 ** Nlinks1
    
    sym2 = LinkSymList[Nlinks2]
    no_terms2 = 2 ** Nlinks2
    
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
        
        coordi = lattice_coord_iso(i, vol=vol)
        coord_list.append(coordi)
        
        stagPhase1 = phase_func1(coordi)
        stagPhase2 = phase_func2(coordi)
       
        momPhase1 = mom_character(mom1, coordi[1:], N, dag1)
        momPhase2 = mom_character(mom2, coordi[1:], N, dag2)
        
        phase_mom_arr1[i] = stagPhase1 * momPhase1
        phase_mom_arr2[i] = stagPhase2 * momPhase2

            
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

                   
    ###Idea here is that the 2 point is 
    
    ### Idea here is that the contraction is  e^{ip1n1}stagphase(n1)e^{ip2n2}stagphase(n2) D(n1{+d1}|n2{+d2}) D(n2{+d2}|n1{+d1})
    ### The shift d1/d2 is on RHS of the prop(1 or 2) when the operator had a dagger, LHS when it did not
    ### prop1(n2|n1), prop2(n1|n2). The standard is to have the dagger at the source x, and no dagger at sink. 
    ### Hence prop1(n2+d2|n1+d1), prop2(n1|n2)
    
    
    prop1 = apply_shift_bd_and_phase(prop=prop1, coord_iso_dict_l=coord_iso_dict2, coord_iso_dict_r=coord_iso_dict1,
                                     bd_dict_l = shift_bd_dict2, bd_dict_r = shift_bd_dict1,
                                     no_terms_l = no_terms2, no_terms_r = no_terms1,
                                     dag_l = dag2, dag_r = dag1,
                                     phase_l = phase_mom_arr2, phase_r = phase_mom_arr1,
                                     matrix_dim=matrix_dim)
    
    
    prop2 = apply_shift_bd_and_phase(prop=prop2, coord_iso_dict_l=coord_iso_dict1, coord_iso_dict_r=coord_iso_dict2,
                                     bd_dict_l = shift_bd_dict1, bd_dict_r = shift_bd_dict2,
                                     no_terms_l = no_terms1, no_terms_r = no_terms2,
                                     dag_l = dag1, dag_r = dag2,
                                     phase_l = phase_mom_arr1, phase_r = phase_mom_arr2,
                                     matrix_dim=matrix_dim)

    
    
    for timej in temporal_range:
        lowerbnd = timej * spatial_dim
        upperbnd = (timej+1) * spatial_dim

        firstSpatialSumPhaseMom1 = prop1[:,lowerbnd:upperbnd] @ prop2[lowerbnd:upperbnd,:]
        
        for timek in temporal_range:
            lowerbnd = timek * spatial_dim
            upperbnd = (timek+1) * spatial_dim

            second_spatial_sum[timej,timek]= np.trace(firstSpatialSumPhaseMom1[lowerbnd:upperbnd,lowerbnd:upperbnd])    
       
    ### result is ntxntxnt array with [t2,t1] where t2 is the tie up time and t1 is source time in the milc code.
    
    return second_spatial_sum*3
