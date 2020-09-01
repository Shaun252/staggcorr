import numpy as np

### Function which returns lattice coord [t,x,y,z] from a 2x2 matrix index 'position'
def lattice_coord_iso(position, vol):
    a = [0 for i in vol]
    i = 0
    while position:
        position,r = divmod(position, vol[::-1][i])
        i += 1
        a[-i] = r
    return a

### Function which returns a 2x2 matrix index from a lattice coord [t,x,y,z] 'position'
def lattice_coord_iso_reverse(position, vol):
    res = position[::-1][0]
    for i, coord in enumerate(position[::-1][1:]):
        res += coord * np.prod(vol[::-1][:i+1])
    return res


def add_lattice_periodic(coord1, coord2, vol):
    
    lattice_coord = [(a + b) % v for a, b, v in zip(coord1, coord2, vol)]
    
    return lattice_coord

def subtract_lattice_periodic(coord1, coord2, vol):
    
    lattice_coord = [(a - b) % v for a, b, v in zip(coord1, coord2, vol)]

    return lattice_coord

def add_lattice_antiperiodic_txyz(coord1, coord2, vol):

    lattice_coord = [(a + b) % v for a, b, v in zip(coord1, coord2, vol)]
    
    if (coord1[0] + coord2[0]) // vol[0] % 2 ==0:
        return lattice_coord, 1
    else:
        return lattice_coord, -1

def subtract_lattice_antiperiodic_txyz(coord1, coord2, vol):

    lattice_coord = [(a - b) % v for a, b, v in zip(coord1, coord2, vol)]
    
    if (coord1[0] - coord2[0]) // vol[0] % 2 ==0:
        return lattice_coord, 1
    else:
        return lattice_coord, -1
