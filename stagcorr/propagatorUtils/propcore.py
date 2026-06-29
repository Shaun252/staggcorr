import time
from stagcorr.propagatorUtils.propIO import saveProp, loadProp
from stagcorr.propagatorUtils.StagDiracOP import staggered_operator
import numpy as np
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


def apply_shift_bd_and_phase(
        prop,                   # D(x,y) matrix with (matrix_dim x matrix_dim), rows = first index, cols = second index
        coord_iso_dict_l,       # dict: link_dirs -> list_of_shifted_indices (for the left index)
        coord_iso_dict_r,       # dict: link_dirs -> list_of_shifted_indices (for the right index)
        bd_dict_l,              # dict: link_dirs -> list_of_bdphases aligned with those left shifted indices
        bd_dict_r,              # dict: link_dirs -> list_of_bdphases aligned with those left shifted indices
        no_terms_l,             # number of symmetric terms (2^Nlinks) or 1 for left shift
        no_terms_r,             # number of symmetric terms (2^Nlinks) or 1 for right shift
        dag_l,                  # True  => operator corresponding to left index was daggered (no shift)
                                # False => operator corresponding to left index was not daggered (shift)
        dag_r,                  # True  => operator corresponding to right index was daggered (shift)
                                # False  => operator corresponding to right index was daggered (shift)                      
        phase_l,                # 1D array length N: staggered*momentum phase for row (first-index) position (unshifted row phase)
        phase_r,                # 1D array length N: staggered*momentum phase for col (second-index) position (unshifted col phase)
        matrix_dim              #
    ):
    """
    Applys shifts to propagator, then applies bd-phase to the SHIFTED index and apply phase_rows/phase_cols to the UNSHIFTED index.
    Returns (complex128) array same shape as prop.

    """
    
    prop_tmp = np.zeros((matrix_dim, matrix_dim), dtype=np.complex128)
    
    ### left index phase and shifting
    if dag_l and no_terms_l != 1:
        # no shift left + left multiply by momentum
        prop_tmp = (phase_l * prop.T).T
        
    elif no_terms_l != 1:
        # shift left index
        for link_dirs, sym_shifted_index_list_l in coord_iso_dict_l.items():
            prop_tmp += (bd_dict_l[link_dirs] * prop[sym_shifted_index_list_l,:].T).T
        prop_tmp /= no_terms_l
    else:
        # do nothing
        prop_tmp = prop

    ### right index phase and shifting, and local treatmeant
    if no_terms_r == 1 or not dag_r:
        # no shift right + right multiply by momentum
        prop_tmp2 = phase_r * prop_tmp
        
    else:
        # shift right
        prop_tmp2 = np.zeros((matrix_dim, matrix_dim), dtype=np.complex128)
        for link_dirs, sym_shifted_index_list_r in coord_iso_dict_r.items():
            prop_tmp2 += bd_dict_r[link_dirs] * prop_tmp[:,sym_shifted_index_list_r]
        prop_tmp2 /= no_terms_r

    return prop_tmp2