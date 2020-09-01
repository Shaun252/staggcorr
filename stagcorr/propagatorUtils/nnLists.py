from stagcorr.lattice.navigation import *

########################################################################################################################################
#### Functions to create nearest neighbour lists for the staggered operator in 2x2 matrix form

def generateAntiPeriodicNearestNeighbourLists(STCoordn, unit_vec_list, unit_vec_list2, unit_vec_list3, nPlus1IndexList, nMinus1IndexList, nPlus2IndexList, nMinus2IndexList, nPlus3IndexList, nMinus3IndexList, vol):
    
    
    nPlus1, bdryPhaseP1 = add_lattice_antiperiodic_txyz(STCoordn, unit_vec_list[0], vol=vol)
    nPlus1Index = lattice_coord_iso_reverse(nPlus1, vol=vol)
    nPlus1IndexList.append(nPlus1Index)
    
    nMinus1, bdryPhaseM1 = subtract_lattice_antiperiodic_txyz(STCoordn, unit_vec_list[0], vol=vol)
    nMinus1Index = lattice_coord_iso_reverse(nMinus1, vol=vol)
    nMinus1IndexList.append(nMinus1Index)
    
    nPlus2 = add_lattice_periodic(STCoordn, unit_vec_list2[0], vol=vol)
    nPlus2Index = lattice_coord_iso_reverse(nPlus2, vol=vol)
    nPlus2IndexList.append(nPlus2Index)

    nMinus2 = subtract_lattice_periodic(STCoordn, unit_vec_list2[0], vol=vol)
    nMinus2Index = lattice_coord_iso_reverse(nMinus2, vol=vol)
    nMinus2IndexList.append(nMinus2Index)
    
    nPlus3, bdryPhaseP3 = add_lattice_antiperiodic_txyz(STCoordn, unit_vec_list3[0], vol=vol)
    nPlus3Index = lattice_coord_iso_reverse(nPlus3, vol=vol)
    nPlus3IndexList.append(nPlus3Index)

    nMinus3, bdryPhaseM3 = subtract_lattice_antiperiodic_txyz(STCoordn, unit_vec_list3[0], vol=vol)
    nMinus3Index = lattice_coord_iso_reverse(nMinus3, vol=vol)
    nMinus3IndexList.append(nMinus3Index)
    
    for mu, unit_vec in enumerate(unit_vec_list[1:]):
        
        mu+=1

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

        nPlus3 = add_lattice_periodic(STCoordn, unit_vec_list3[mu], vol=vol)
        nPlus3Index = lattice_coord_iso_reverse(nPlus3, vol=vol)
        nPlus3IndexList.append(nPlus3Index)

        nMinus3 = subtract_lattice_periodic(STCoordn, unit_vec_list3[mu], vol=vol)
        nMinus3Index = lattice_coord_iso_reverse(nMinus3, vol=vol)
        nMinus3IndexList.append(nMinus3Index)
        
    return bdryPhaseP1, bdryPhaseM1, bdryPhaseP3, bdryPhaseM3

def generatePeriodicNearestNeighbourLists(STCoordn, unit_vec_list, unit_vec_list2, unit_vec_list3, nPlus1IndexList, nMinus1IndexList, nPlus2IndexList, nMinus2IndexList, nPlus3IndexList, nMinus3IndexList, vol):
    
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

        nPlus3 = add_lattice_periodic(STCoordn, unit_vec_list3[mu], vol=vol)
        nPlus3Index = lattice_coord_iso_reverse(nPlus3, vol=vol)
        nPlus3IndexList.append(nPlus3Index)

        nMinus3 = subtract_lattice_periodic(STCoordn, unit_vec_list3[mu], vol=vol)
        nMinus3Index = lattice_coord_iso_reverse(nMinus3, vol=vol)
        nMinus3IndexList.append(nMinus3Index)
        
    return 1, 1, 1, 1