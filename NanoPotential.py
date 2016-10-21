from prody import *
import numpy as np

def findPotential(pdbFile, pos = ([0], [0], [0]), cutoff = 0, field = 'CHARMM'):
    
    if field == 'AMBER':
        epsilon = 9.3368e-22
        R = 3.5692
    else:
        epsilon = 7.1689e-22
        R = 3.7606
    
    nCarb, r2 = getCoords(pdbFile)
    
    numR = len(pos[0])
    numTheta = len(pos[1])
    numZ = len(pos[2])
    
    V = np.zeros((numR, numTheta, numZ))
    x = np.zeros((numR, numTheta, numZ))
    y = np.zeros((numR, numTheta, numZ))
    z = np.zeros((numR, numTheta, numZ))
    
    for rIndex in range(0, numR):
    
        for thetaIndex in range(0, numTheta):
        
            for zIndex in range(0, numZ):
        
                x[rIndex, thetaIndex, zIndex] = pos[0][rIndex] * np.cos(pos[1][thetaIndex])
                y[rIndex, thetaIndex, zIndex] = pos[0][rIndex] * np.sin(pos[1][thetaIndex])
                z[rIndex, thetaIndex, zIndex] = pos[2][zIndex]
                
                r1 = (x[rIndex, thetaIndex, zIndex], y[rIndex, thetaIndex, zIndex], z[rIndex, thetaIndex, zIndex])
                V[rIndex, thetaIndex, zIndex] = sumLJ(nCarb, r1, r2, cutoff, epsilon = epsilon, R = R)
        
    return x, y, z, V
    
def getCoords(pdbFile):

    cnt = parsePDB(pdbFile, chain = 'TUB')
    coords = cnt.getCoords()
    
    nCarb = coords.shape[0]
    r2 = (coords[:, 0], coords[:, 1], coords[:, 2])
    
    return nCarb, r2

def sumLJ(nCarb, (x1, y1, z1), (x2, y2, z2), cutoff, epsilon = 7.1689e-22, R = 3.7606):

    r = np.sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)
    V = LJ(r, epsilon, R)
    if cutoff > 0:
        V[r > cutoff] = 0
    
    return sum(V)

def LJ(r, epsilon, R):

    V = epsilon * ((R/r)**12 - 2*(R/r)**6)
    
    return V