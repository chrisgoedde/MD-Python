from prody import *
import numpy as np

def getCoordTimeSeries(paths, runFile, start, finish, atomName):

    pdbFile = paths['tfinal'] + runFile + '.pdb'
    dcdFile = paths['tfinal'] + runFile + '.dcd'

    # Now, open up the DCD file and get the ensembles for the carbons and oxygens

    structure = parsePDB(pdbFile)
    ensemble = parseDCD(dcdFile)

    atom = structure.select(atomName)

    ensemble.setAtoms(atom)
    atomC = ensemble.getCoordsets()
    
    if start == 0:
        firstTime = 0
    else:
        firstTime = start - 1
        
    lastTime = finish
    
    atomX = atomC[firstTime:lastTime, :, 0]
    atomY = atomC[firstTime:lastTime, :, 1]
    atomZ = atomC[firstTime:lastTime, :, 2]
    
    if start == 0:
    
        atomC0 = getCoords(pdbFile, atomName)
        atomX0 = atomC0[:, 0]
        atomY0 = atomC0[:, 1]
        atomZ0 = atomC0[:, 2]

        atomX = np.concatenate((np.array([atomX0]), atomX))
        atomY = np.concatenate((np.array([atomY0]), atomY))
        atomZ = np.concatenate((np.array([atomZ0]), atomZ))
    
    return atomX, atomY, atomZ

def getCoords(fileName, atomName):

    structure = parsePDB(fileName)
    
    o = structure.select(atomName)
    coords = o.getCoords()
    
    return coords
    
def getVelTimeSeries(paths, runFile, start, finish, atomName):

    velFactor = 20.45482706
    
    pdbFile = paths['tfinal'] + runFile + '.pdb'
    dcdFile = paths['tfinal'] + runFile + '.veldcd'

    # Open up the DCD file and get the ensembles for the carbons and oxygens

    structure = parsePDB(pdbFile)
    ensemble = parseDCD(dcdFile)

    atom = structure.select(atomName)

    ensemble.setAtoms(atom)
    atomC = ensemble.getCoordsets()
    
    if start == 0:
        firstTime = 0
    else:
        firstTime = start - 1
        
    lastTime = finish
    
    atomX = velFactor * atomC[firstTime:lastTime, :, 0]
    atomY = velFactor * atomC[firstTime:lastTime, :, 1]
    atomZ = velFactor * atomC[firstTime:lastTime, :, 2]
        
    return atomX, atomY, atomZ

import matplotlib.pyplot as plt

from matplotlib import animation, rc

def animateAtoms(xVar, yVar, xLabel = 'water index', yLabel = r'water offset/$\lambda$', \
    title = '', style = 'o'):

    (numFrames, numAtoms) = yVar.shape

    fig, ax = plt.subplots()

    ax.set_xlim((np.min(xVar), np.max(xVar)))
    ax.set_ylim((np.min(yVar), np.max(yVar)))
    
    plt.xlabel(xLabel, fontsize=14)
    plt.ylabel(yLabel, fontsize=14)
    # plt.axis('scaled')

    line, = ax.plot([], [], 'o')

    anim = animation.FuncAnimation(fig, animateAtomFrame, fargs=(xVar, yVar, line), \
                                   frames=numFrames, interval=20, blit=True)

    return anim
    
# animation function. This is called sequentially
def animateAtomFrame(i, xVar, yVar, line):

    if len(xVar.shape) == 1:
        x = xVar
    else:
        x = xVar[i,:]
        
    y = yVar[i,:]
    line.set_data(x, y)
    return (line,)


