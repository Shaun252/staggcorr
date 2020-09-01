from stagcorr.gauge.gaugeIO import *
import time
from math import sqrt, factorial
from itertools import permutations
import sympy as sp
from stagcorr.gauge.reuniterise import reuniterise

class U:
   
    def __init__(self, position, direction, dagger=False):
        
        self.position = sp.sympify(position)
        self.direction = sp.sympify(direction)
        self.dagger = dagger
    def display_form(self):
        if self.dagger == False:
            return r'$U_\{}({})$'.format(self.direction, self.position)
        else:
            return r'$U^\dagger_\{}({})$'.format(self.direction, self.position)
        
    

def second_deriv(Umu, direction):

    n = Umu.position
    mu = Umu.direction
    
    nPrho  = n + direction
    nPmu  = n + mu
    nMrho = n-direction
    nMrhoPmu = n-direction+mu
    
    Urho_n = U(n,direction)
    Umu_nPrho = U(n+direction,mu)
    Urho_nPmu_dag = U(n+mu,direction, dagger=True)
    
    Urho_nMrho_dag = U(n-direction,direction, dagger=True)
    Umu_nMrho = U(n-direction,mu)
    Urho_nMrhoPmu = U(n-direction+mu, direction)
    
    
    
    
    return [[1/4, Urho_n, Umu_nPrho, Urho_nPmu_dag], [1/4,Urho_nMrho_dag, Umu_nMrho, Urho_nMrhoPmu], [-2/4, Umu]]



def ApplySecondDeriv(MOP, direction):
    
    full_list = []
    for term in MOP:
        coeff = term[0]
        
        for i, Link in enumerate(term[1:]):
            second_derivLink = second_deriv(Link, direction=direction)
            
            for secondderivTerm in second_derivLink:
                coeffSD = secondderivTerm[0]
                secondderivTermCompleteList = [coeff*coeffSD]
            
                
                full_list.append(secondderivTermCompleteList + term[1:i+1] + secondderivTerm[1:] + term[i+2:])
                
    return full_list


def applyFmuSinglePerm():
    
    Umu = U(sp.sympify('n'), sp.sympify('mu'))
    
    full_list = [[1, Umu],]
    
    rho = sp.sympify('rho')
    sigma = sp.sympify('sigma')
    nu = sp.sympify('nu')
    
    second_deriv_rho = second_deriv(Umu=Umu, direction=rho)
    second_deriv_sigma = second_deriv(Umu=Umu, direction=sigma)
    second_deriv_nu = second_deriv(Umu=Umu, direction=nu)
    
    second_deriv_nu_rho = ApplySecondDeriv(MOP=second_deriv_rho, direction=nu)
    second_deriv_sigma_rho = ApplySecondDeriv(MOP=second_deriv_rho, direction=sigma)
    second_deriv_nu_sigma = ApplySecondDeriv(MOP=second_deriv_sigma, direction=nu)
    
    second_deriv_nu_sigma_rho = ApplySecondDeriv(MOP=second_deriv_nu_rho, direction=sigma)
    
    full_list = full_list + second_deriv_rho + second_deriv_sigma + second_deriv_nu + second_deriv_nu_rho + second_deriv_sigma_rho + second_deriv_nu_sigma + second_deriv_nu_sigma_rho
    
    return full_list
    

def DisplayMultiLinkOP(MOP):
    
    for term in MOP:
        coeff = term[0]
        
        displayStr = str(coeff)
        
        for Link in term[1:]:
            displayStr+=Link.display_form()
            
        display(Latex(displayStr))
        

def convert_sympy_op(gFieldOperation):
    full_list = []
    n_symbol = sp.sympify('n')
    
    for term in gFieldOperation:
        coeff = [term[0]]
        term_list = []
        for U in term[1:]:
            if U.position ==  n_symbol:
                shift_list = []
            
            else:
                link_position_shifts = U.position
                shift_list = []
                #print(link_position_shifts)
                for shift in link_position_shifts.args:
                    if shift != n_symbol:
                        shift_list.append(str(shift))
            
            term_list.append((shift_list, str(U.direction), U.dagger))
            
        full_list.append(coeff + term_list)
    
    return full_list

def pos_to_unit_vec(position_shifts, STCoordN, mu,  rho, sigma, nu, unit_vec_list, vol):
    
    dir_dict = {'mu' : mu,
                'nu' : nu,
                'rho' : rho,
                'sigma' : sigma
               }
    
    totalShift = [0 for i in STCoordN]
    for shift in position_shifts:
        
        if shift[0] == '-':
            shift_val = dir_dict[shift[1:]]
            totalShift = [m-n for (m,n) in zip(totalShift, unit_vec_list[shift_val])]
        else:
            
            shift_val = dir_dict[shift]
            totalShift = [m+n for (m,n) in zip(totalShift, unit_vec_list[shift_val])]
    STCoordM = add_lattice_periodic(totalShift, STCoordN, vol)
    
    #print(STCoordN) 
    #print(STCoordM)
    
    m = lattice_coord_iso_reverse(STCoordM, vol)
    
    return m


def applyGaugeProduct(gFieldOperationList, gFieldArr, gFieldArrList, n, STCoordN, mu, rho, sigma, nu, unit_vec_list, vol, SUN=3):
    
    dir_dict = {'mu' : mu,
                'nu' : nu,
                'rho' : rho,
                'sigma' : sigma
               }
                
    
    final_res = np.zeros((SUN, SUN), dtype=np.complex128) 
    
    for term in gFieldOperationList:
        coeff = term[0]
        
        res = coeff * np.identity(SUN, dtype=np.complex128) 
        for U in term[1:]:
            link_position_shifts = U[0]
            link_direction = U[1]
            link_dagger = U[2]
            
            dir_no = dir_dict[link_direction]
            
            if len(link_position_shifts)==0:
                #print(link_position_shifts)
                #print(dir_no)
                m = n
            
            else:
                #print(link_position_shifts)
                #print(dir_no)
                m = pos_to_unit_vec(STCoordN=STCoordN, position_shifts=link_position_shifts, mu=mu,  rho=rho, sigma=sigma, nu=nu, unit_vec_list=unit_vec_list, vol=vol)
                
            if link_dagger:
                
                res = res @ gFieldArrList[m][dir_no].conj().T
            else:
                res = res @ gFieldArrList[m][dir_no]
            #print(gFieldArrList[m][dir_no])
            #print('')
        final_res += res
            
    return final_res


def FmuSmeared_link(gFieldArr, gFieldArrList, FmuFieldOperationList, n, STCoordN, mu, unit_vec_list, Nlinks, vol):
    
    dir_list = list(range(Nlinks))
    dir_list.remove(mu)
    symFac = factorial(len(dir_list))
    
    symmetrised_dir_list = permutations(dir_list)
    
    
    res = np.zeros((3,3), dtype=np.complex128)
    
    for dir_perm in symmetrised_dir_list:
        #print('dir_perm')
        #print(dir_perm)
        
        dir_perm_res = applyGaugeProduct(gFieldOperationList=FmuFieldOperationList, gFieldArr=gFieldArr, gFieldArrList=gFieldArrList, n=n, STCoordN=STCoordN, mu=mu, rho=dir_perm[0], sigma=dir_perm[1], nu=dir_perm[2], unit_vec_list=unit_vec_list, vol=vol)
         
        
            
        res += dir_perm_res
    
    return res / symFac

def lepage_term(gFieldArr, gFieldArrList, n, STCoordN, mu, nPlus1IndexList,  nMinus1IndexList, nPlus2IndexList,  nMinus2IndexList, unit_vec_list, Nlinks, vol):
    
    res = np.zeros((3,3), dtype = 'complex128')
    
    unitVecmu = unit_vec_list[mu]
    
    nPlus1mu = nPlus1IndexList[mu]
    nMinus1mu = nMinus1IndexList[mu]
    
    for rho in range(Nlinks):
        if rho != mu:
            
            Minus1rhoPlus1muList =  [ a-b  for a,b in zip(unitVecmu, unit_vec_list[rho])]
            nMinus1rhoPlus1muCoord = add_lattice_periodic(Minus1rhoPlus1muList, STCoordN, vol)
            nM1rhoP1mu = lattice_coord_iso_reverse(nMinus1rhoPlus1muCoord, vol)
            
            Minus2rhoPlus1muList =  [ a-2*b  for a,b in zip(unitVecmu, unit_vec_list[rho])]
            nMinus2rhoPlus1muCoord = add_lattice_periodic(Minus2rhoPlus1muList, STCoordN, vol)
            nM2rhoP1mu = lattice_coord_iso_reverse(nMinus2rhoPlus1muCoord, vol)
            
            Plus1rhoPlus1muList =  [ a+b  for a,b in zip(unitVecmu, unit_vec_list[rho])]
            nPlus1rhoPlus1muCoord = add_lattice_periodic(Plus1rhoPlus1muList, STCoordN, vol)
            nP1rhoP1mu = lattice_coord_iso_reverse(nMinus1rhoPlus1muCoord, vol)
            
            nPlus1rho = nPlus1IndexList[rho]
            nMinus1rho = nMinus1IndexList[rho]
            
            nPlus2rho = nPlus2IndexList[rho]
            nMinus2rho = nMinus2IndexList[rho]
            
        
            res += gFieldArrList[n][rho] @ gFieldArrList[nPlus1rho][rho] @ gFieldArrList[nPlus2rho][mu] @ gFieldArrList[nP1rhoP1mu][rho].conj().T @ gFieldArrList[nMinus1mu][rho].conj().T
            
            res += gFieldArrList[nMinus1rho][rho].conj().T @ gFieldArrList[nMinus2rho][rho].conj().T @ gFieldArrList[nMinus2rho][mu] @ gFieldArrList[nM2rhoP1mu][rho] @ gFieldArrList[nM1rhoP1mu][rho]
            
            res -= 2*gFieldArr
    
    return res /2

def HISQleftSmear(gFieldArr, gFieldArrList, FmuFieldOperationList, n, STCoordN, mu, nPlus1IndexList,  nMinus1IndexList, nPlus2IndexList,  nMinus2IndexList, unit_vec_list, unit_vec_list2, Nlinks, vol):
    
    lepage_correction = lepage_term(gFieldArr, gFieldArrList, n, STCoordN, mu, nPlus1IndexList,  nMinus1IndexList, nPlus2IndexList,  nMinus2IndexList, unit_vec_list, Nlinks, vol)
    
    FmuU = FmuSmeared_link(gFieldArr, gFieldArrList, FmuFieldOperationList, n, STCoordN, mu, unit_vec_list, Nlinks, vol)
    
    HISQleftSmearLink =  FmuU - lepage_correction
    
    return HISQleftSmearLink


def FHISQ(gFieldArr, gFieldArrList, FmuFieldOperationList, n, STCoordN, mu, nPlus1IndexList,  nMinus1IndexList, nPlus2IndexList,  nMinus2IndexList, unit_vec_list, unit_vec_list2, Nlinks, vol):
    
    FmuU = FmuSmeared_link(gFieldArr, gFieldArrList, FmuFieldOperationList, n, STCoordN, mu, unit_vec_list, Nlinks, vol)
    X = reuniterise(FmuU)
    W = HISQleftSmear(X, gFieldArrList, FmuFieldOperationList, n, STCoordN, mu, nPlus1IndexList,  nMinus1IndexList, nPlus2IndexList,  nMinus2IndexList, unit_vec_list, unit_vec_list2, Nlinks, vol)
    
    return W, X

def generatePeriodicNearestNeighbourGaugeSmearLists(STCoordn, unit_vec_list, unit_vec_list2, nPlus1IndexList, nMinus1IndexList,  nPlus2IndexList, nMinus2IndexList, vol):
    
    for mu, unit_vec in enumerate(unit_vec_list):

        nPlus1 = add_lattice_periodic(STCoordn, unit_vec, vol=vol)
        nPlus1Index = lattice_coord_iso_reverse(nPlus1, vol=vol)
        nPlus1IndexList.append(nPlus1Index)

        nMinus1 = subtract_lattice_periodic(STCoordn, unit_vec, vol=vol)
        nMinus1Index = lattice_coord_iso_reverse(nMinus1, vol=vol)
        nMinus1IndexList.append(nMinus1Index)
        
        nPlus2 = add_lattice_periodic(STCoordn, unit_vec_list2[mu], vol=vol)
        nPlus2Index = lattice_coord_iso_reverse(nPlus2, vol=vol)
        nPlus2IndexList.append(nPlus2Index)

        nMinus2 = subtract_lattice_periodic(STCoordn, unit_vec_list2[mu], vol=vol)
        nMinus2Index = lattice_coord_iso_reverse(nMinus2, vol=vol)
        nMinus2IndexList.append(nMinus2Index)
        
    return None

def ApplyHISQSmearing(gFieldArrList, vol):
    dim=len(vol)
    dim_range = range(dim)
    physRange = range(len(gFieldArrList))
    
    unit_vec_list = []
    unit_vec_list2 = []
    for mu in dim_range:
        unit_vec = [0 for i in dim_range]
        unit_vec[mu] = 1
        unit_vec_list.append(unit_vec[:])
        unit_vec[mu] = 2
        unit_vec_list2.append(unit_vec[:])
    
    
    WGaugeArrList = []
    XGaugeArrList = []
            
    fullFmuSmear = convert_sympy_op(applyFmuSinglePerm())
    
    for n in physRange:
        
        STCoordN = lattice_coord_iso(n, vol=vol)
        
        ### NNeighbour lists
        nPlus1IndexList = []
        nMinus1IndexList = []
        
        nPlus2IndexList = []
        nMinus2IndexList = []

        generatePeriodicNearestNeighbourGaugeSmearLists(STCoordn=STCoordN, unit_vec_list=unit_vec_list, unit_vec_list2=unit_vec_list2, nPlus1IndexList=nPlus1IndexList, nMinus1IndexList=nMinus1IndexList,  nPlus2IndexList=nPlus2IndexList, nMinus2IndexList=nMinus2IndexList, vol=vol)
        
        WGaugeArrList.append([])
        XGaugeArrList.append([])
        
        for mu in dim_range:
            W, X = FHISQ(gFieldArr=gFieldArrList[n][mu], gFieldArrList=gFieldArrList, FmuFieldOperationList=fullFmuSmear, n=n, STCoordN=STCoordN, mu=mu, nPlus1IndexList=nPlus1IndexList,  nMinus1IndexList=nMinus1IndexList, nPlus2IndexList=nPlus2IndexList, nMinus2IndexList=nMinus2IndexList, unit_vec_list=unit_vec_list, unit_vec_list2=unit_vec_list2, Nlinks=dim, vol=vol)
            
            WGaugeArrList[n].append(W)
            XGaugeArrList[n].append(X)
            
    return WGaugeArrList, XGaugeArrList



def constructHISQGFields(gField, vol):
    if gField == None:
        WGaugeArrList = XGaugeArrList = generate_unit_gauge_field_vec(vol=vol)
    elif type(gField)==type(""):
        gFieldX = gField+"_X"
        gFieldW = gField+"_W"
        try:
            XGaugeArrList, volG = load_gauge_field_vec(filename=gFieldX, vol=vol)
            WGaugeArrList, volG = load_gauge_field_vec(filename=gFieldW, vol=vol)
            print("HISQ Smeared " +gField+ " W & X Fields Loaded")
        except:
            gFieldArrList, volG = load_gauge_field_vec(filename=gField, vol=vol)
            print("Gauge Field Loaded")
            GSmearstart = time.time()
            WGaugeArrList, XGaugeArrList = ApplyHISQSmearing(gFieldArrList=gFieldArrList, vol=vol)
            GSmearstop = time.time()
            print("HISQ Smeared " +gField+ " W & X Fields Generated: " + str(GSmearstop-GSmearstart) +"s")
            #save_gauge_field(filename=gFieldX, gfieldArrList=XGaugeArrList, vol=vol)
            #save_gauge_field(filename=gFieldW, gfieldArrList=WGaugeArrList, vol=vol)
            
        if volG != vol:
                print("Gauge volum does not match")
    else:
        GSmearstart = time.time()
        WGaugeArrList, XGaugeArrList = ApplyHISQSmearing(gFieldArrList=gField, vol=vol)
        GSmearstop = time.time()
        print("HISQ Smeared W & X Fields Generated: " + str(GSmearstop-GSmearstart) +"s")
                
    return XGaugeArrList, WGaugeArrList