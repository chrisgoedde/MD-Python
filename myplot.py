import matplotlib
import matplotlib.pyplot as plt
from matplotlib import animation, rc
import os.path
import os
import numpy as np
from datetime import datetime

import analysis as an
import myutil as my

def animatePotential(pD, run, rings = (), cutoff = 0, save = False):

    nFrames = run['carbPos'].shape[1]
    time = run['time']
    
    V = []
    for i in range(0, nFrames):
        # startTime = datetime.now()
        z, v = an.findPotential(pD, run, frame = i, rings = rings, cutoff = cutoff)
        # print('Found Potential for Frame ' + str(i) + ' in time ' + str(datetime.now()-startTime))
        V.append(v)
    
    fig, ax = plt.subplots()

    ax.set_xlim((np.min(z), np.max(z)))
    ax.set_ylim((np.min(V), np.max(V)))

    plt.xlabel('z (A)')
    plt.ylabel('On-Axis Potential (J)')
    theTitle = makePlotTitle(pD)
    if cutoff != 0:
        theTitle = theTitle.replace('\n', ', Cutoff = ' + str(cutoff) + ' A' + '\n')
    else:
        theTitle = theTitle.replace('\n', ', No Cutoff' + '\n')
    plt.title(theTitle)    
    plt.text(0.1, 0.9, 'Time = ' + '{:.2f}'.format(time[0]) + ' ns', transform = plt.gca().transAxes)

    line, = ax.plot([], [])

    anim = animation.FuncAnimation(fig, animatePotentialFrame, fargs=(z, V, time, line), \
                                   frames=nFrames, interval=40, blit=True)

    # Note that fps must match the interval given above ...
    
    if save:
        savePath = my.makePaths(pD)['pictures']
        if not os.path.exists(savePath):
            os.makedirs(savePath)
    
        if cutoff == 0:
            cutString = '-NoCutoff'
        else:
            cutString = '-Cutoff-' + str(cutoff)
        if rings == ():
            fileName = 'Potential-All' + cutString
        else:
            fileName = 'Potential-' + str(rings[0]) + '-' + str(rings[1]) + cutString
        print('Saving ' + fileName + ' to ' + savePath)
        anim.save(savePath + fileName + '.mp4', fps=25, dpi = 300)
    
    return anim

def animatePotentialFrame(i, xVar, yVar, time, line):

    y = yVar[i]

    line.set_data(xVar, y)
    plt.gca().texts[-1].set_text('Time = ' + '{:.2f}'.format(time[i]) + ' ns')
    return (line,)

def animateWaterOffset(pD, run, style = 'o', save = False):
    
    yVar = run['waterOffset']
    time = run['time']

    (numFrames, numAtoms) = yVar.shape
    xVar = np.arange(numAtoms)

    fig, ax = plt.subplots()

    ax.set_xlim((np.min(xVar), np.max(xVar)))
    ax.set_ylim((-0.5, 1.5))
    
    plt.xlabel('Water Index', fontsize=14)
    plt.ylabel(r'Water Offset/$\lambda$', fontsize=14)
    plt.title(makePlotTitle(pD))
    plt.text(0.1, 0.9, 'Time = ' + '{:.2f}'.format(time[0]) + ' ns', transform = plt.gca().transAxes)

    line, = ax.plot([], [], 'o')

    anim = animation.FuncAnimation(fig, animateOffsetFrame, fargs=(xVar, yVar, time, line), \
                                   frames=numFrames, interval=40, blit=True)
                                   
    # Note that fps must match the interval given above ...
    
    if save:
        savePath = my.makePaths(pD)['pictures']
        if not os.path.exists(savePath):
            os.makedirs(savePath)
        fileName = 'WaterOffset.mp4'
        print('Saving ' + fileName + ' to ' + savePath)
        anim.save(savePath + fileName, fps=25, dpi = 300)
    
    return anim

# animation function. This is called sequentially

def animateOffsetFrame(i, xVar, yVar, time, line):

    if len(xVar.shape) == 1:
        x = xVar
    else:
        x = xVar[i,:]
        
    y = yVar[i,:]
    while np.mean(y) > 1:
        y = y - 1
    line.set_data(x, y)
    plt.gca().texts[-1].set_text('Time = ' + '{:.2f}'.format(time[i]) + ' ns')
    return (line,)
    
def plotPotential(pD, run, frame = 0, rings = (), cutoff = 0, save = False):

    # startTime = datetime.now()
    z, V = an.findPotential(pD, run, frame = frame, rings = rings, cutoff = cutoff)
    # print('Found Potential in time ' + str(datetime.now()-startTime))

    f = plt.figure(figsize = (8,6))
    plt.plot(z, V)
    plt.xlabel('z (A))')
    plt.ylabel('On-Axis Potential (J)')
    theTitle = makePlotTitle(pD)
    if cutoff != 0:
        theTitle = theTitle.replace('\n', ', Cutoff = ' + str(cutoff) + ' A' + '\n')
    else:
        theTitle = theTitle.replace('\n', ', No Cutoff' + '\n')
    plt.title(theTitle)
    plt.grid()

    if save:
        if cutoff == 0:
            cutString = '-NoCutoff'
        else:
            cutString = '-Cutoff-' + str(cutoff)
        if rings == ():
            fileName = 'Potential-All' + cutString
        else:
            fileName = 'Potential-' + str(rings[0]) + '-' + str(rings[1]) + cutString
        savePlot(pD, fileName)
        plt.close(f)
    
def plotSolitonPosition(pD, run, save = False):

    nT = len(run['time'])
    solitonPos = np.zeros((nT, 2))

    for frame in range(0, nT):

        s = an.findSolitons(run['waterOffset'][frame])
        solitonPos[frame, :] = s

    f = plt.figure(figsize = (8,6))
    plt.plot(run['time'], solitonPos[:,0], 'o')
    plt.plot(run['time'], solitonPos[:,1], 'o')
    plt.gca().set_ylim((0, pD['N0']))
    plt.xlabel('time (ns)')
    plt.ylabel('Soliton Position (Water Index)')
    plt.title(makePlotTitle(pD))
    plt.grid()

    if save:
        savePlot(pD, 'SolitonPosition')
        plt.close(f)
    
def plotWaterMeanTemperature(pD, run, save = False):

    f = plt.figure(figsize = (8,6))
    plt.plot(run['time'], 0.5*(run['waterMeanTemp'][0] + run['waterMeanTemp'][1]), 'o')
    plt.plot(run['time'], run['waterMeanTemp'][2], 'o')
    plt.xlabel('time (ns)')
    plt.ylabel('Water Mean Temperature (K)')
    plt.title(makePlotTitle(pD))
    plt.legend(['Transverse Temperature', 'Axial Temperature'], loc = 'best')
    plt.grid()

    if save:
        savePlot(pD, 'WaterMeanTemperature')
        plt.close(f)

def plotSystemTemperature(pD, run, save = False):

    f = plt.figure(figsize = (8,6))
    plt.plot(run['time'], run['carbTemp'], 'o')
    plt.xlabel('time (ns)')
    plt.ylabel('System Temperature (K)')
    plt.title(makePlotTitle(pD))
    plt.plot(run['time'], run['waterTemp'], 'o')
    plt.legend(['Carbon', 'Water'], loc = 'best')
    plt.grid()
    
    if save:
        savePlot(pD, 'SystemTemperature')
        plt.close(f)

def plotWaterMeanPosition(pD, run, save = False, ylim = ()):

    f = plt.figure(figsize = (8,6))
    plt.plot(run['time'], run['waterMeanPos'][2], 'o')
    plt.xlabel('time (ns)')
    plt.ylabel('Water Mean Position (A)')
    plt.title(makePlotTitle(pD))
    if ylim != ():
        plt.gca().set_ylim(ylim)
    plt.grid()
    
    if save:
        savePlot(pD, 'WaterMeanPosition')
        plt.close(f)
    
def plotWaterMeanVelocity(pD, run):

    f = plt.figure(figsize = (8,6))
    plt.plot(run['time'], 10*run['waterMeanVel'][2], 'o')
    plt.xlabel('time (ns)')
    plt.ylabel('Mean Water Vel (A/ns)')
    plt.title(makePlotTitle(pD))
    plt.grid()

def showAnimation(picPathList, fileName):

    files = ''
    for pp in picPathList:
        files = files + '"' + pp + fileName + '" '
    os.system('open ' + files)
                
def showPicture(picPathList, fileName):

    files = ''
    for pp in picPathList:
        print('Opening picture ' + '"' + pp + fileName + '"')
        files = files + '"' + pp + fileName + '" '
    os.system('open ' + files)

def makePlotTitle(pD):

    theTitle = '({0:1}, {1:1})'.format(pD['n'], pD['m']) + ' nanotube, S = ' \
        + my.num2Str(pD['S']) + ', N = ' + str(pD['N0']) \
        + '\n' + 'T = ' + str(pD['Temperature (K)']) + ' K' + ', Forcing = ' \
        + str(pD['Force (pN)']) + ' pN, Restraint = ' + str(pD['Restraint'])

    return theTitle

def savePlot(pD, fileName):

    savePath = my.makePaths(pD)['pictures']
    if not os.path.exists(savePath):
        os.makedirs(savePath)

    print('Saving ' + fileName + '.pdf' + ' to ' + savePath)
    plt.savefig(savePath + fileName + '.pdf')

