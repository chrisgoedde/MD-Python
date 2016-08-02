from prody import *
import numpy as np
import myutil as my

def getCoordTimeSeries(paths, runFile, start, finish, atomName):

    pdbFile = paths['data'] + runFile + '.pdb'
    dcdFile = paths['data'] + runFile + '.dcd'

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
    
        atomC0 = getCoords(paths, runFile, atomName)
        atomX0 = atomC0[:, 0]
        atomY0 = atomC0[:, 1]
        atomZ0 = atomC0[:, 2]

        atomX = np.concatenate((np.array([atomX0]), atomX))
        atomY = np.concatenate((np.array([atomY0]), atomY))
        atomZ = np.concatenate((np.array([atomZ0]), atomZ))
    
    return np.array([ atomX, atomY, atomZ ])

def getCoords(paths, runFile, atomName):

    pdbFile = paths['data'] + runFile + '.pdb'
    
    structure = parsePDB(pdbFile)
    
    o = structure.select(atomName)
    coords = o.getCoords()
    
    return coords
    
def getVelTimeSeries(paths, runFile, start, finish, atomName):

    velFactor = 20.45482706
    
    pdbFile = paths['data'] + runFile + '.pdb'
    dcdFile = paths['data'] + runFile + '.veldcd'

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
    
    # Convert velocities to m/s
    
    atomX = 100 * velFactor * atomC[firstTime:lastTime, :, 0]
    atomY = 100 * velFactor * atomC[firstTime:lastTime, :, 1]
    atomZ = 100 * velFactor * atomC[firstTime:lastTime, :, 2]
        
    return np.array([ atomX, atomY, atomZ ])
    
def analyzeRun(pD):

    paths = my.makePaths(pD)

    # First we need the coordinates from the pdb and dcd files
    # We'll begin with the minimization trajectory, if present
    # We'll store the calculated values for the minimization trajectory
    # in the dictionary minRun
    
    minRun = {}
    
    if pD['Run Type'] == 'New':

        start = 0
        finish = pD['Min Duration']/pD['outputFreq']

        # Calculate the time of each snapshot, in ns

        minRun['time'] = np.arange(start, finish+1) * pD['dt (fs)'] * pD['outputFreq'] / 1000000.

        print('Minimization run from record ' + str(start) + ' to ' + str(finish))
        
        getValues(paths, pD, minRun, start, finish)

        print('Initial water radius = ' + str(minRun['waterMeanRadius']) + ' A')
        print('Initial carbon radius = ' + str(minRun['carbMeanRadius']) + ' A')
    
    # Now move on to the production part of the run
    # We'll save the calculated values in the dictionary proRun
    
    proRun = {}

    if pD['Run Type'] == 'New':
        start = finish + 1
        finish = finish + pD['Duration']/pD['outputFreq']
        proRun['time'] = np.arange(start, finish+1) * pD['dt (fs)'] * pD['outputFreq'] / 1000000.
    else:
        start = 1
        finish = 1000
        proRun['time'] = np.arange(start, finish+1) * pD['dt (fs)'] * 100 / 1000000.

    print('Production run from record ' + str(start) + ' to ' + str(finish))
    
    getValues(paths, pD, proRun, start, finish)

    print('Initial water radius = ' + str(proRun['waterMeanRadius']) + ' A')
    print('Initial carbon radius = ' + str(proRun['carbMeanRadius']) + ' A')
    
    if pD['Run Type'] == 'New':
        calculatePhysicalValues(pD, minRun)

    calculatePhysicalValues(pD, proRun)
    
    return minRun, proRun
    
def getValues(paths, pD, run, start, finish):

    run['waterPos'] = getCoordTimeSeries(paths, pD['File Name'], start, finish, 'oxygen')
    run['waterVel'] = getVelTimeSeries(paths, pD['File Name'], start, finish, 'oxygen')

    run['carbPos'] = getCoordTimeSeries(paths, pD['File Name'], start, finish, 'carbon')
    run['carbVel'] = getVelTimeSeries(paths, pD['File Name'], start, finish, 'carbon')

    run['waterRadius'] = np.sqrt(run['waterPos'][0]**2 + run['waterPos'][1]**2)
    run['waterMeanRadius'] = np.mean(run['waterRadius'][0,:])
    run['carbRadius'] = np.sqrt(run['carbPos'][0]**2 + run['carbPos'][1]**2)
    run['carbMeanRadius'] = np.mean(run['carbRadius'][0,:])

def calculatePhysicalValues(pD, run):

    # Length of each carbon ring, in A

    ringLength = 2.456

    # Calculate the offsets for the oxygens, and express them in dimensionless units
    # The equilibrium position for a water molecule is approximately ringLength/4

    eqZ = ringLength/2 + (np.arange(pD['N0']+pD['S'])-pD['N0']/2) * ringLength
    run['waterOffset'] = (run['waterPos'][2] - eqZ)/ringLength

    # Carbon and water masses, in kg

    carbMass = 1.994e-26
    waterMass = 2.991e-26

    # Calculate the KE and temperature of the water and carbon atoms

    run['carbKE'] = 0.5 * carbMass * run['carbVel']**2
    run['waterKE'] = 0.5 * waterMass * run['waterVel']**2
    
    kB = 1.381e-23

    run['carbTemp'] = 2 * run['carbKE'] / kB
    run['waterTemp'] = 2 * run['waterKE'] / kB
    
    # Calculate some mean values by summing over the relevant atoms or molecules

    run['waterMeanVel'] = np.mean(run['waterVel'], axis = 2)
    run['waterMeanPos'] = np.mean(run['waterPos'], axis = 2)
    
    run['carbMeanKE'] = np.mean(run['carbKE'], axis = 2)
    run['waterMeanKE'] = np.mean(run['waterKE'], axis = 2)
    
    run['carbMeanTemp'] = np.mean(run['carbTemp'], axis = 2)
    run['waterMeanTemp'] = np.mean(run['waterTemp'], axis = 2)
    
def calculateFlowRate(pD, run):

    x = run['time']
    y = run['waterMeanPos'][2]
    fit = np.polyfit(x, y, 1)
    
    return fit[0]

def calculateWork(pD, run):

    KE = np.sum(run['waterKE'], axis = 0)
    totalKE = np.sum(KE, axis = 1)
    deltaKE = np.diff(totalKE)
    deltaZ = np.diff(run['waterMeanPos'][2])
    work = (pD['N0'] + pD['S']) * pD['Force (pN)'] * deltaZ * 1e-22
    
    return deltaKE, work
    
def movingAverage(values, window, offset):

    avgV = values
    for i in range(1, window):
        leftV = np.roll(values, i)
        leftV[:i] = leftV[:i] - offset
        rightV = np.roll(values, -i)
        rightV[-i:] = rightV[-i:] + offset
        avgV = avgV + leftV + rightV
        
    return avgV/(2*window + 1)

def findPeaks(values):

    peakIndex = []
    for i in range(1, len(values)+1):
        if (values[0] < values[1]) and (values[2] < values[1]) and (values[1] > 0.004):
            peakIndex.append(np.mod(i, len(values)))
        values = np.roll(values, -1)
        
    return peakIndex
    
def findSolitons(offset):

    avgOffset = movingAverage(offset, 10, 1)
    diffOffset = np.diff(avgOffset)
    avgDiff = movingAverage(diffOffset, 10, 0)

    p = findPeaks(avgDiff)
    
    # plt.plot(np.arange(len(avgDiff)), avgDiff, 'o')

    i1 = np.argmax(avgDiff[p])
    M1 = p[i1]
    
    if len(p) > 1:
    
        del p[i1]
        i2 = np.argmax(avgDiff[p])
        M2 = p[i2]
        return np.array([ M1, M2 ])
        
    else:
    
        return np.array([ M1, M1 ])

def findPotential(pD, run, frame = 0, nPoints = 50, rings = (), cutoff = 0):

    wavelength = 2.456
    L = wavelength * pD['N0']
    kB = 1.381e-23
    
    if rings == ():
        ringRange = range(0, pD['N0'] * nPoints)
    else:
        ringRange = range(rings[0] * nPoints, rings[1] * nPoints)
    
    zVec = np.array(ringRange) * (wavelength/nPoints) - L/2

    coords = run['carbPos'][:, frame, :]
    nCarb = coords.shape[1]

    V = np.zeros(zVec.shape)
    r2 = (coords[0], coords[1], coords[2])
    
    for index in range(0, len(zVec)):
    
        r1 = (0, 0, zVec[index])
    
        V[index] = sumLJ(nCarb, r1, r2, L, cutoff)

    return zVec, V
    
def sumLJ(nCarb, (x1, y1, z1), (x2, y2, z2), L, cutoff):

    fixL = np.zeros(z2.shape)
    fixL[z2-z1 > L/2] = L
    fixL[z2-z1 < -L/2] = -L
    r = np.sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2+fixL)**2)
    V = LJ(r)
    if cutoff > 0:
        V[r > cutoff] = 0
    
    return sum(V)

def LJ(r):

    epsilon = 2.267e-21
    R = 3.7606
    
    V = epsilon * ((R/r)**12 - 2*(R/r)**6)
    
    return V