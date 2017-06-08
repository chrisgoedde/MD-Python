import numpy as np

def findPotential(theta):

    k = [ 0.06450, 0.14975, 0.09458, 0.11251 ]
    n = [ 2, 3, 4, 5 ]
    delta = [ 0, np.pi, 0, 0 ]
    
    # Conversion factor from kcal/mol to J
    
    # conv = 6.9477e-21
    conv = 1
    
    # Initialization for V
    
    V = 0
    
    for i in range(0, 4):
    
        V = V + conv * k[i] * (1 + np.cos(n[i]*theta - delta[i]))
        
    return V
    
def calcDihedral(atoms1, atoms2, atoms3, atoms4, radian=False):
    """Returns the dihedral angle between atoms in degrees."""

    if not isinstance(atoms1, Atomic):
        raise TypeError('atoms1 must be an Atomic instance')
    if not isinstance(atoms2, Atomic):
        raise TypeError('atoms2 must be an Atomic instance')
    if not isinstance(atoms3, Atomic):
        raise TypeError('atoms3 must be an Atomic instance')
    if not isinstance(atoms4, Atomic):
        raise TypeError('atoms4 must be an Atomic instance')
    if not (atoms1.numAtoms() == atoms2.numAtoms() ==
            atoms3.numAtoms() == atoms4.numAtoms()):
        raise ValueError('all arguments must have same number of atoms')

    return getDihedral(atoms1._getCoords(), atoms2._getCoords(),
                       atoms3._getCoords(), atoms4._getCoords(), radian)

def getDihedral(coords1, coords2, coords3, coords4, radian=False):
    """Returns the dihedral angle in degrees."""

    a1 = coords2 - coords1
    a2 = coords3 - coords2
    a3 = coords4 - coords3

    v1 = np.cross(a1, a2)
    v1 = v1 / (v1 * v1).sum(-1)**0.5
    v2 = np.cross(a2, a3)
    v2 = v2 / (v2 * v2).sum(-1)**0.5
    porm = np.sign((v1 * a3).sum(-1))
    rad = np.arccos((v1*v2).sum(-1) / ((v1**2).sum(-1) * (v2**2).sum(-1))**0.5)
    if not porm == 0:
        rad = rad * porm
    if radian:
        return rad
    else:
        return rad * 180 / np.pi

def getDihedralChain(coords, radian=False):
    """Returns the dihedral angle in degrees for an alkane chain."""

    n = coords.shape[0]
    c1 = coords[0:(n-3)]
    c2 = coords[1:(n-2)]
    c3 = coords[2:(n-1)]
    c4 = coords[3:n]
    
    a1 = c2 - c1
    a2 = c3 - c2
    a3 = c4 - c3

    v1 = np.cross(a1, a2)
    v1 = v1 / np.linalg.norm(v1, axis = -1, keepdims = True)
    v2 = np.cross(a2, a3)
    v2 = v2 / np.linalg.norm(v2, axis = -1, keepdims = True)
    
    porm = np.sign((v1 * a3).sum(-1))
    
    temp = (v1*v2).sum(-1) / ((v1**2).sum(-1) * (v2**2).sum(-1))**0.5
    
    temp[temp < -1] = -1
    temp[temp > 1] = 1
    
    rad = np.arccos(temp)
    rad[porm == -1] = -rad[porm == -1]
    rad[rad < 0] = rad[rad < 0] + 2*np.pi

    if radian:
        return rad
    else:
        return rad * 180 / np.pi
