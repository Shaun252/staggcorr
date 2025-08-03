"""Staggered fermion phase functions and operator utilities.

This module implements the staggered fermion formalism including phase functions,
operator mappings, and momentum characters for lattice QCD calculations.
Follows the TXYZ coordinate convention.
"""

import sympy as sp
import numpy as np
########################################################################################################################################
#### Staggered Phase functions with convention txyz

def bar(gamma_vec, mod = True):
    """Compute the bar operation for staggered phases.
    
    For each component i, sums all other components j≠i.
    Used in staggered fermion phase calculations.
    
    Parameters
    ----------
    gamma_vec : array_like
        Input gamma vector (typically coordinate or momentum)
    mod : bool, default=True
        Apply modulo 2 operation to result
        
    Returns
    -------
    numpy.ndarray
        Bar operation result: bar_i = sum_{j≠i} gamma_j
    """
    new_gamma_vec = np.zeros((len(gamma_vec)), dtype = type(gamma_vec[0]))
    for i, val1 in enumerate(gamma_vec):
        for j, val2 in enumerate(gamma_vec):
            if j != i:
                new_gamma_vec[i]+=val2
                
    if mod == True:
        return np.mod(new_gamma_vec, 2)
    else:
        return new_gamma_vec


def less_than(gamma_vec, mod = True):
    """Compute the less-than operation for staggered phases.
    
    For each component i, sums all components j<i.
    Essential for eta phase calculations in staggered fermions.
    
    Parameters
    ----------
    gamma_vec : array_like
        Input gamma vector (typically coordinate or momentum)
    mod : bool, default=True
        Apply modulo 2 operation to result
        
    Returns
    -------
    numpy.ndarray
        Less-than operation result: lt_i = sum_{j<i} gamma_j
    """
    new_gamma_vec = np.zeros((len(gamma_vec)), dtype = type(gamma_vec[0]))
    for i, val1 in enumerate(gamma_vec):
        for j, val2 in enumerate(gamma_vec):
            if j < i:
                new_gamma_vec[i]+=val2
                
    if mod == True:
        return np.mod(new_gamma_vec, 2)
    else:
        return new_gamma_vec

def greater_than(gamma_vec, mod = True):
    new_gamma_vec = np.zeros((len(gamma_vec)), dtype = type(gamma_vec[0]))
    for i, val1 in enumerate(gamma_vec):
        for j, val2 in enumerate(gamma_vec):
            if j > i:
                new_gamma_vec[i]+=val2
                
    if mod == True:
        return np.mod(new_gamma_vec, 2)
    else:
        return new_gamma_vec

def epsilon_phase(n_vec):
    exponent = np.sum(n_vec)
    
    return np.power(-1, exponent)

def eta_phase(n_vec):
    """Compute the eta phase factor for staggered fermions.
    
    Calculates (-1)^(less_than operation), fundamental to
    staggered fermion Dirac operator construction.
    
    Parameters
    ----------
    n_vec : array_like
        Coordinate vector [t, x, y, z]
        
    Returns
    -------
    numpy.ndarray
        Phase factors: (-1)^(less_than(n_vec))
    """
    exponent =  less_than(n_vec)
    
    return np.power(-1, exponent)

def zeta_phase(n_vec):
    exponent = greater_than(n_vec)
    
    return np.power(-1, exponent)

def phase_func_txyz(phase_exp):
    
    n_vec = sp.var("n_0 n_1 n_2 n_3")
    
    phase_exp= sp.sympify(phase_exp)
    
    phase_exponent = sum(phase_exp.free_symbols)
    
    phase_exp_func = sp.lambdify([n_vec], np.power(-1, phase_exponent))
    
    return phase_exp_func


##########################################################################################################################################
#### Useful functions

def sympy_mod(expression):
    #### For mod 2 of sympy expressions eg (-1)^(3n_1 + 2n_2 + 4_n3) = (-1)^(n_1)
    new_expression = 0
    if expression ==sp.sympify(1) or expression == sp.sympify(-1):
        new_expression += sp.sympify(1)
        
    elif expression != 0:
        if len(str(expression).split(" "))==1:
            if len(expression.args) == 2: 
                if -1 in  expression.args:
                    new_expression = -1 * expression
                else:
                    for symbol in expression.free_symbols:
                        
                        coeff = expression.coeff(symbol)
                        if coeff % 2 == 1:
                            new_expression += symbol
            else:
                new_expression = expression
        else:
            for arg in expression.args:
                if arg.is_constant():
                    new_expression += arg % 2

            for symbol in expression.free_symbols:
                coeff = expression.coeff(symbol)
                if coeff % 2 == 1:
                    new_expression += symbol
    else:
        new_expression = 0
    return new_expression

def add(gamma_vec1, gamma_vec2):
    """Add two gamma vectors with modulo 2 arithmetic.
    
    Parameters
    ----------
    gamma_vec1, gamma_vec2 : array_like
        Input gamma vectors to add
        
    Returns
    -------
    numpy.ndarray
        Sum of vectors with mod 2 applied to each component
        
    Notes
    -----
    Used for staggered fermion phase calculations where
    operations are performed in Z_2 arithmetic.
    """
    gamma_new = np.add(gamma_vec1, gamma_vec2)
    gamma_new = np.mod(gamma_new, 2)
    return gamma_new


##########################################################################################################################################
#### Converting from gamma spin taste operators to phase and shift / link operators


def phase_shift_operator(spin, taste):
    """Convert gamma operators to staggered phase and shift.
    
    Transforms continuum Dirac gamma matrices (spin ⊗ taste) into
    staggered fermion phase factors and coordinate shifts.
    
    Parameters
    ----------
    spin : str
        Spin Dirac structure ('G1', 'G5', 'GX', 'GY', 'GZ', 'GT', etc.)
    taste : str
        Taste structure ('G1', 'G5', 'GX', 'GY', 'GZ', 'GT', etc.)
        
    Returns
    -------
    tuple
        (phase_expression, shift_vector)
        - phase_expression: SymPy expression for staggered phase
        - shift_vector: numpy array of coordinate shifts [t,x,y,z]
        
    Notes
    -----
    Uses TXYZ convention [x_0, x_1, x_2, x_3] = [t, x, y, z].
    Based on Follana HISQ paper appendix formulation.
    """
    ### convention = [x_0, x_1, x_2, x_3]
    ### See Follana HISQ paper appendix
    
    n_vec = sp.var("n_0 n_1 n_2 n_3")
    
    ### Any gamma with G"AB", i.e two Letters/numbers I don't under the definitions and the extra phases from
    ### anti commuting in the dict are not understood fully
    gamma_dict = {
    "G1":  [[0,0,0,0],1],
    "GX":  [[0,1,0,0],1],
    "GY":  [[0,0,1,0],1],
    "GZ":  [[0,0,0,1],1],
    "GT":  [[1,0,0,0],1],
    "G5":  [[1,1,1,1],1],
    "GYZ": [[0,0,1,1],1], 
    "GZX": [[0,1,0,1],-1],
    "GXY": [[0,1,1,0],1],
    "GXT": [[1,1,0,0],-1],
    "GYT": [[1,0,1,0],-1],
    "GZT": [[1,0,0,1],-1],
    "G5X": [[1,0,1,1],-1],
    "G5Y": [[1,1,0,1],1],
    "G5Z": [[1,1,1,0],-1],
    "G5T": [[0,1,1,1],1]}
    
    
    m, sign_m = gamma_dict[spin]
    s, sign_s = gamma_dict[taste]
    
    phase_exp1 = np.dot(n_vec,add(less_than(s), greater_than(m)))
    
    phase_exp2 = np.dot(m,less_than(add(s,m)))

    var_indp_phase = (-1) ** phase_exp2 * sign_m * sign_s
    
    if var_indp_phase == 1:
        total_phase_exp = phase_exp1
    else:
        total_phase_exp = phase_exp1 + sp.sympify(1)
    
    shift = add(m, s)
    
    
    return total_phase_exp, shift

##########################################################################################################################################
#### mom character


def mom_character(mom_three_vec, group_element, N):
    
    pdotx = np.dot(mom_three_vec, group_element)
    
    char_val = np.exp(1j * 2 *np.pi * pdotx /  N)
    
    return char_val 
