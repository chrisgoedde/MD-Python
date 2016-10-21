from prody import *
import numpy as np
import myutil as my
import cPickle as pickle
import os.path

def getCoordTimeSeries(paths, runFile, start, finish, atomName):

    pdbFile = paths['data'] + runFile + '.pdb'
    dcdFile = paths['data'] + runFile + '.dcd'

    # Now, open up the DCD file and get the ensembles for the carbons and oxygens

    structure = parsePDB(pdbFile)
    ensemble = parseDCD(dcdFile)

    atom = structure.select(atomName)
    
    if atom is None:
    
        return np.array([])
        
    else:
    
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
        
            if len(atomC0) != 0:
        
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
    if o is None:
        coords = np.array([])
    else:
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

    if atom is None:
    
        return np.array([])
        
    else:
    
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
    
def analyzeRun(pD, doPickle = False, overWrite = False):

    paths = my.makePaths(pD)
    
    pickleFileName = paths['data'] + pD['File Name'] + '.pickle'
    
    if (not overWrite) and os.path.exists(pickleFileName):
        pickleFile = open(pickleFileName, 'rb')
        minRun, proRun = pickle.load(pickleFile)
        return minRun, proRun
    
    if pD['S'] != -pD['N0']:
        selection = ('water', 'carbon')
    else:
        selection = ('carbon', )

    # First we need the coordinates from the pdb and dcd files
    # We'll begin with the minimization trajectory, if present
    # We'll store the calculated values for the minimization trajectory
    # in the dictionary minRun
    
    minRun = {}
    proRun = {}
    
    if pD['Run Type'] == 'New':

        start = 0
        finish = pD['Min Duration']/pD['outputFreq']

        # Calculate the time of each snapshot, in ns

        minRun['time'] = np.arange(start, finish+1) * pD['dt (fs)'] * pD['outputFreq'] / 1000000.

        print('Minimization run from record ' + str(start) + ' to ' + str(finish))
        
        getValues(paths, pD, selection, minRun, start, finish)
        calculatePhysicalValues(pD, selection, minRun)
        minRun['zAxis'], minRun['axisPot'] = findPotential(pD, minRun)
        minRun['initAxisPotAmp'] = (max(minRun['axisPot'])-min(minRun['axisPot']))/2
        minRun['initAxisPotMean'] = (max(minRun['axisPot'])+min(minRun['axisPot']))/2
        proRun['initAxisPotAmp'] = (max(minRun['axisPot'])-min(minRun['axisPot']))/2
        proRun['initAxisPotMean'] = (max(minRun['axisPot'])+min(minRun['axisPot']))/2
        for sel in selection:
            print('Initial ' + sel + ' radius = ' + str(minRun[sel, 'meanRadius']) + ' A')
    
    # Now move on to the production part of the run
    # We'll save the calculated values in the dictionary proRun
    
    if pD['Run Type'] == 'New':
        start = finish + 1
        finish = finish + pD['Duration']/pD['outputFreq']
        proRun['time'] = np.arange(start, finish+1) * pD['dt (fs)'] * pD['outputFreq'] / 1000000.
    else:
        start = 1
        finish = 1000
        proRun['time'] = np.arange(start, finish+1) * pD['dt (fs)'] * 100 / 1000000.

    print('Production run from record ' + str(start) + ' to ' + str(finish))
    
    getValues(paths, pD, selection, proRun, start, finish)
    calculatePhysicalValues(pD, selection, proRun)
    proRun['zAxis'], proRun['axisPot'] = findPotential(pD, proRun)
    
    for sel in selection:
        print('Initial ' + sel + ' radius = ' + str(proRun[sel, 'meanRadius']) + ' A')
    
    if doPickle:
        pickleFile = open(pickleFileName, 'wb')
        pickle.dump((minRun, proRun), pickleFile, protocol = -1)
    
    return minRun, proRun
    
def getValues(paths, pD, selection, run, start, finish):

    atomName = { 'water': 'oxygen', 'carbon': 'carbon' }

    for sel in selection:

        run[sel, 'pos'] = getCoordTimeSeries(paths, pD['File Name'], start, finish, atomName[sel])
        run[sel, 'vel'] = getVelTimeSeries(paths, pD['File Name'], start, finish, atomName[sel])

        run[sel, 'radius'] = np.sqrt(run[sel, 'pos'][0]**2 + run[sel, 'pos'][1]**2)
        run[sel, 'meanRadius'] = np.mean(run[sel, 'radius'][0,:])
        
def calculatePhysicalValues(pD, selection, run):

    # Length of each carbon ring, in A
    # Carbon and water masses, in kg

    ringLength = 2.456
    mass = { 'water': 2.991e-26, 'carbon': 1.994e-26 }
    kB = 1.381e-23

    for sel in selection:
    
        # Calculate the offsets for the oxygens, and express them in dimensionless units
        # The equilibrium position for a water molecule is approximately ringLength/4

        if sel == 'water':
            eqZ = ringLength/2 + (np.arange(pD['N0']+pD['S'])-pD['N0']/2) * ringLength
            run[sel, 'offset'] = (run[sel, 'pos'][2] - eqZ)/ringLength

        # Calculate the KE and temperature of the water and carbon atoms

        run[sel, 'KE'] = 0.5 * mass[sel] * run[sel, 'vel']**2
    
        run[sel, 'temp'] = 2 * run[sel, 'KE'] / kB
    
        # Calculate some mean values by summing over the relevant atoms or molecules

        run[sel, 'meanVel'] = np.mean(run[sel, 'vel'], axis = 2)
        run[sel, 'meanPos'] = np.mean(run[sel, 'pos'], axis = 2)    
        run[sel, 'meanKE'] = np.mean(run[sel, 'KE'], axis = 2)
        run[sel, 'meanTemp'] = np.mean(run[sel, 'temp'], axis = 2)    
    
    
def calculateFlowRate(pD, run):

    x = run['time']
    y = run['water', 'meanPos'][2]
    fit = np.polyfit(x, y, 1)
    
    return fit[0]

def calculateWork(pD, run):

    KE = np.sum(run['water', 'KE'], axis = 0)
    totalKE = np.sum(KE, axis = 1)
    deltaKE = np.diff(totalKE)
    deltaZ = np.diff(run['water', 'meanPos'][2])
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

def findPotential(pD, run, frame = 0, nPoints = 50, rings = (), cutoff = 0, field = 'CHARMM'):

    wavelength = 2.456
    L = wavelength * pD['N0']
    
    if field == 'AMBER':
        epsilon = 9.3368e-22
        R = 3.5692
    else:
        epsilon = 7.1689e-22
        R = 3.7606
    
    if rings == ():
        ringRange = range(0, pD['N0'] * nPoints)
    else:
        ringRange = range(rings[0] * nPoints, rings[1] * nPoints)
    
    zVec = np.array(ringRange) * (wavelength/nPoints) - L/2

    coords = run['carbon', 'pos'][:, frame, :]
    nCarb = coords.shape[1]

    V = np.zeros(zVec.shape)
    r2 = (coords[0], coords[1], coords[2])
    
    for index in range(0, len(zVec)):
    
        r1 = (0, 0, zVec[index])
    
        V[index] = sumLJ(nCarb, r1, r2, L, cutoff, epsilon = epsilon, R = R)
        
    amplitude = (max(V) - min(V))/2
    mean = (max(V) + min(V))/2

    return zVec, V
    
def sumLJ(nCarb, (x1, y1, z1), (x2, y2, z2), L, cutoff, epsilon = 7.1689e-22, R = 3.7606):

    fixL = np.zeros(z2.shape)
    fixL[z2-z1 > L/2] = L
    fixL[z2-z1 < -L/2] = -L
    r = np.sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2+fixL)**2)
    V = LJ(r, epsilon, R)
    if cutoff > 0:
        V[r > cutoff] = 0
    
    return sum(V)

def LJ(r, epsilon, R):

    V = epsilon * ((R/r)**12 - 2*(R/r)**6)
    
    return V