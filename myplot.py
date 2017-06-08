import matplotlib
import matplotlib.pyplot as plt
from matplotlib import animation, rc
import os.path
import os
import numpy as np
from datetime import datetime

import analysis as an
import myutil as my

def animateDihedralChain(pD, run, save = False):

def animateDihedralFrame(i, bN, phi, time, line):
    
    plt.gca().texts[-1].set_text('Time = ' + '{:.1f}'.format(float(time[i])) + ' ps')

    line.set_data(bN, phi[i,:])
    return (line,)



def animatePotential(pD, run, rings = (), cutoff = 0, save = False, field = 'CHARMM'):

    nFrames = run['carbon', 'pos'].shape[1]
    time = run['time']
    
    V = []
    for i in range(0, nFrames):
        # startTime = datetime.now()
        z, v = an.findPotential(pD, run, frame = i, rings = rings, \
            cutoff = cutoff, field = field)
        v = v / 6.9477e-21
        # print('Found Potential for Frame ' + str(i) + ' in time ' + str(datetime.now()-startTime))
        V.append(v)
    
    fig, ax = plt.subplots()

    ax.set_xlim((np.min(z), np.max(z)))
    
    if 'initAxisPotAmp' in run:
    
        yLower = (run['initPotMean']-1.5*run['initPotAmp']) / 6.9477e-21
        yUpper = (run['initPotMean']+1.5*run['initPotAmp']) / 6.9477e-21
        ax.set_ylim((yLower, yUpper))
    
    else:
    
        ax.set_ylim((np.min(V), np.max(V)))

    plt.xlabel('z (A)')
    plt.ylabel('On-Axis Potential (kcal/mol)')
    theTitle = makePotentialPlotTitle(pD, cutoff, field, run['initPotAmp'])
    plt.title(theTitle)    
    plt.text(0.1, 0.9, 'Time = ' + '{:.2f}'.format(time[0]) + ' ns', transform = plt.gca().transAxes)

    line, = ax.plot([], [])

    anim = animation.FuncAnimation(fig, animatePotentialFrame, fargs=(z, V, time, line), \
                                   frames=nFrames, interval=40, blit=True)

    # Note that fps must match the interval given above ...
    
    if save:
        savePath = my.makeWaterPaths(pD)['pictures']
        if not os.path.exists(savePath):
            os.makedirs(savePath)
    
        if cutoff == 0:
            cutString = '-NoCutoff'
        else:
            cutString = '-Cutoff-' + str(cutoff)
        if rings == ():
            fileName = 'Potential-All' + cutString + '-' + field
        else:
            fileName = 'Potential-' + str(rings[0]) + '-' + str(rings[1]) + cutString + '-' + field
        print('Saving ' + fileName + ' to ' + savePath)
        anim.save(savePath + fileName + '.mp4', fps=25, dpi = 300)
        plt.close(fig)
    
    return anim

def animatePotentialFrame(i, xVar, yVar, time, line):

    y = yVar[i]

    line.set_data(xVar, y)
    plt.gca().texts[-1].set_text('Time = ' + '{:.2f}'.format(time[i]) + ' ns')
    return (line,)

def animateWaterOffset(pD, run, style = '.', save = False):
    
    yVar = run['water', 'offset']
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

    line, = ax.plot([], [], '.')

    anim = animation.FuncAnimation(fig, animateOffsetFrame, fargs=(xVar, yVar, time, line), \
                                   frames=numFrames, interval=40, blit=True)
                                   
    # Note that fps must match the interval given above ...
    
    if save:
        savePath = my.makeWaterPaths(pD)['pictures']
        if not os.path.exists(savePath):
            os.makedirs(savePath)
        fileName = 'WaterOffset.mp4'
        print('Saving ' + fileName + ' to ' + savePath)
        anim.save(savePath + fileName, fps=25, dpi = 300)
        plt.close(fig)
    
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
    
def plotWaterOffset(pD, run, frame = (), save = False):
    
    offset = run['water', 'offset']
    time = run['time']
    
    (numFrames, numAtoms) = offset.shape
    x = np.arange(numAtoms)

    if frame == ():
        frame = (numFrames-1,)
    elif type(frame) is not tuple:
        frame = (frame, )
        
    f = plt.figure(figsize = (8,6))

    plt.gca().set_xlim((np.min(x), np.max(x)))
    plt.gca().set_ylim((-0.5, np.abs(pD['S'])+0.5))
    
    plt.xlabel('Water Index')
    plt.ylabel(r'Water Offset/$\lambda$')

    for i in range(0, len(frame)):
    
        while np.mean(offset[frame[i]]) > 1:
            offset[frame[i]] = offset[frame[i]] - 1
        while np.mean(offset[frame[i]]) < 0:
            offset[frame[i]] = offset[frame[i]] + 1

        plt.plot(x, offset[frame[i]], '.')

    plt.title(makePlotTitle(pD))
    if save:
        savePlot(pD, 'WaterOffset')
        plt.close(f)

def plotPotential(pD, run, frame = 0, rings = (), cutoff = 0, field = 'CHARMM', save = False):

    # startTime = datetime.now()
    z, V = an.findPotential(pD, run, frame = frame, rings = rings, \
        cutoff = cutoff, field = field)
    # print('Found Potential in time ' + str(datetime.now()-startTime))

    V = V / 6.9477e-21

    if 'initAxisPotAmp' in run:
    
        yLower = (run['initPotMean']-1.5*run['initPotAmp']) / 6.9477e-21
        yUpper = (run['initPotMean']+1.5*run['initPotAmp']) / 6.9477e-21
        ax.set_ylim((yLower, yUpper))
    
    else:
    
        ax.set_ylim((np.min(V), np.max(V)))
    
    f = plt.figure(figsize = (8,6))
    plt.plot(z, V)
    plt.gca().set_ylim((yLower, yUpper))
    plt.xlabel('z (A)')
    plt.ylabel('On-Axis Potential (kcal/mol)')
    theTitle = makePotentialPlotTitle(pD, cutoff, field, ampV, showAmplitude = True)
    plt.title(theTitle)
    plt.grid()

    if save:
        if cutoff == 0:
            cutString = '-NoCutoff'
        else:
            cutString = '-Cutoff-' + str(cutoff)
        if rings == ():
            fileName = 'Potential-All' + cutString + '-' + field
        else:
            fileName = 'Potential-' + str(rings[0]) + '-' + str(rings[1]) + cutString + '-' + field
        savePlot(pD, fileName)
        plt.close(f)
    
def plotSolitonPosition(pD, run, save = False):

    nT = len(run['time'])
    solitonPos = np.zeros((nT, 2))

    for frame in range(0, nT):

        s = an.findSolitons(run['water', 'offset'][frame])
        solitonPos[frame, :] = s

    f = plt.figure(figsize = (8,6))
    plt.plot(run['time'], solitonPos[:,0], '.')
    plt.plot(run['time'], solitonPos[:,1], '.')
    plt.gca().set_ylim((0, pD['N0']))
    plt.xlabel('time (ns)')
    plt.ylabel('Soliton Position (Water Index)')
    plt.title(makePlotTitle(pD))
    plt.grid()

    if save:
        savePlot(pD, 'SolitonPosition')
        plt.close(f)
    
def plotWaterMeanTemperature(pD, run, save = False, xlim = (), ylim = (), fit = False):

    f = plt.figure(figsize = (8,6))
    if fit:
        if pD['Force (pN)'] == 0:
            setBounds = False
        else:
            setBounds = True
        tXY, cXY = an.fitTS(run['time'], 0.5*(run['water', 'meanTemp'][0] + run['water', 'meanTemp'][1]), setBounds)
        tZ, cZ = an.fitTS(run['time'], run['water', 'meanTemp'][2], setBounds)
                
        plt.plot(run['time'], tXY)
        plt.plot(run['time'], tZ)
        
    else:
        plt.plot(run['time'], 0.5*(run['water', 'meanTemp'][0] + run['water', 'meanTemp'][1]), '.')
        plt.plot(run['time'], run['water', 'meanTemp'][2], '.')

    plt.xlabel('time (ns)')
    plt.ylabel('Water Mean Temperature (K)')
    plt.title(makePlotTitle(pD))
    if xlim != ():
        plt.gca().set_xlim(xlim)
    if ylim != ():
        plt.gca().set_ylim(ylim)
    if fit:
        plt.legend(['Transverse Temperature {:.1f} K'.format(cXY[2]), 'Axial Temperature {:.1f} K'.format(cZ[2])], loc = 'best')
    else:
        plt.legend(['Transverse Temperature', 'Axial Temperature'], loc = 'best')
    plt.grid()

    if save:
        if fit:
            savePlot(pD, 'WaterMeanTemperatureFitted')
        else:
            savePlot(pD, 'WaterMeanTemperature')
        plt.close(f)

def plotSystemTemperature(pD, run, save = False, ylim = (), fit = False):

    f = plt.figure(figsize = (8,6))
    if fit:
        if pD['Force (pN)'] == 0:
            setBounds = False
        else:
            setBounds = True
        tC, cC = an.fitTS(run['time'], np.mean(run['carbon', 'meanTemp'],0), setBounds)
        tW, cW = an.fitTS(run['time'], np.mean(run['water', 'meanTemp'],0), setBounds)
                
        plt.plot(run['time'], tC)
        plt.plot(run['time'], tW)
        
    else:
        plt.plot(run['time'], np.mean(run['carbon', 'meanTemp'],0), '.')
        plt.plot(run['time'], np.mean(run['water', 'meanTemp'],0), '.')
        
    plt.xlabel('time (ns)')
    plt.ylabel('System Temperature (K)')
    plt.title(makePlotTitle(pD))
    if ylim != ():
        plt.gca().set_ylim(ylim)
    plt.legend(['Carbon', 'Water'], loc = 'best')
    plt.grid()
    
    if save:
        if fit:
            savePlot(pD, 'SystemTemperatureFitted')
        else:
            savePlot(pD, 'SystemTemperature')
        
        plt.close(f)

def plotWaterMeanPosition(pD, run, save = False, ylim = ()):

    ms = dict(markersize = 3)
    f = plt.figure(figsize = (8,6))
    plt.plot(run['time'], run['water', 'meanPos'][2], '.', **ms)
    plt.xlabel('time (ns)')
    plt.ylabel('Water Mean Position (A)')
    plt.title(makePlotTitle(pD))
    if ylim != ():
        plt.gca().set_ylim(ylim)
    plt.grid()
    
    if save:
        savePlot(pD, 'WaterMeanPosition')
        plt.close(f)
    
def plotWaterMeanVelocity(pD, run, save = False):

    f = plt.figure(figsize = (8,6))
    plt.plot(run['time'], 10*run['water', 'meanVel'][2], '.')
    plt.xlabel('time (ns)')
    plt.ylabel('Water Mean Vel (A/ns)')
    plt.title(makePlotTitle(pD))
    plt.grid()
    
    if save:
        savePlot(pD, 'WaterMeanVelocity')
        plt.close(f)

def plotGroupFlowRate(pD, f, t, runs = (1, 10), frameRange = (5000, 10000), xLim = (5, 10), yLim = (-800, 800), slice = False):

    v = {}
    meanV = {}
    stdV = {}
    n = runs[1]-runs[0]+1
    
    if slice:
        sliceSize = 100
        numSlice = (frameRange[1]-frameRange[0])/sliceSize
    else:
        numSlice = 1

    for force in f:
        for temp in t:
            pD['Force (pN)'] = force
            pD['Temperature (K)'] = temp
        
            if slice:
                v[force, temp] = np.zeros((n, numSlice))
                tSlice = np.zeros(numSlice)
            else:
                v[force, temp] = np.zeros(n)
                
            f = plt.figure(figsize = (8,6))
            for num in range(runs[0]-1, runs[1]):

                pD['File Name'] = 'Run-' + str(num+1)
                mR, pR = an.analyzeRun(pD)
                
                if slice:
                    for i in range(0, numSlice):
                        [ vel, pos ] = np.polyfit(pR['time'][frameRange[0]+i*sliceSize:frameRange[0]+(i+1)*sliceSize], pR['water', 'meanPos'][2, frameRange[0]+i*sliceSize:frameRange[0]+(i+1)*sliceSize], 1)
                        v[force, temp][num, i] = vel
                        tSlice[i] = np.mean(pR['time'][frameRange[0]+i*sliceSize:frameRange[0]+(i+1)*sliceSize])
                
                    plt.plot(tSlice, v[force, temp][num,:])

                else:
                    [ vel, pos ] = np.polyfit(pR['time'][frameRange[0]:frameRange[1]], pR['water', 'meanPos'][2, frameRange[0]:frameRange[1]],1)
                    v[force, temp][num] = vel
                    plt.plot(pR['time'][frameRange[0]:frameRange[1]], pR['water', 'meanPos'][2, frameRange[0]:frameRange[1]]-np.mean(pR['water', 'meanPos'][2, frameRange[0]:frameRange[1]]))
            
            plt.xlabel('time (ns)')
            plt.ylabel('Water Mean Position (A)')
            plt.grid()
            plt.gca().set_xlim(xLim)
            plt.gca().set_ylim(yLim)

            theTitle = makePlotTitle(pD)
        
            meanV[force, temp] = np.mean(v[force, temp])
            stdV[force, temp] = np.std(v[force, temp])/np.sqrt(num*numSlice)
        
            theTitle = theTitle + ', ' + my.makeFlowRateString(meanV[force, temp], stdV[force, temp])
            plt.title(theTitle)
        
            if slice:
                fileName = 'SliceFlowRate'
            else:
                fileName = 'FlowRate'
            savePlot(pD, fileName, groupRuns = True)
            plt.close(f)
        
            print('Velocity = ' + str(meanV[force, temp]) + ' +- ' + str(stdV[force, temp]) + ' A/ns at force = ' + str(force) + ' pN, and T = ' + str(temp) + ' K')
            
def plotGroupTemperature(pD, f, t, runs = (1, 10), frameRange = (5000, 10000), xLim = (5, 10), yLim = (0, 300), atom = 'water'):

    for force in f:
        for temp in t:
            pD['Force (pN)'] = force
            pD['Temperature (K)'] = temp
        
            f = plt.figure(figsize = (8,6))
            for num in range(runs[0]-1, runs[1]):

                pD['File Name'] = 'Run-' + str(num+1)
                mR, pR = an.analyzeRun(pD)
                
                plt.plot(pR['time'][frameRange[0]:frameRange[1]], 0.5*(pR[atom, 'MeanTemp'][0][frameRange[0]:frameRange[1]] + pR[atom, 'MeanTemp'][1][frameRange[0]:frameRange[1]]))

            plt.xlabel('time (ns)')
            plt.ylabel('Transverse Water Temperature (K)')
            plt.grid()
            plt.gca().set_xlim(xLim)
            plt.gca().set_ylim(yLim)

            theTitle = makePlotTitle(pD)
            plt.title(theTitle)
            
            fileName = atom.capitalize() + 'Temperature'
            savePlot(pD, fileName, groupRuns = True)
            plt.close(f)
            
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

def makePotentialPlotTitle(pD, cutoff, field, amplitude, showAmplitude = False):

    if cutoff == 0:
        cutString = 'No Cutoff, '
    else:
        cutString = 'Cutoff = ' + str(cutoff) + ' A, '

    if showAmplitude:
        ampString = 'Amp = ' + '{:.2e}'.format(amplitude) + ' kcal/mol'
    else:
        ampString = 'S = ' + my.num2Str(pD['S'])
         
    theTitle = '({0:1}, {1:1})'.format(pD['n'], pD['m']) + ' nanotube, N = ' + str(pD['N0']) \
        + ', ' + field + ', '+ cutString \
        + ampString \
        + '\n' + 'T = ' + str(pD['Temperature (K)']) + ' K' + ', Forcing = ' \
        + str(pD['Force (pN)']) + ' pN, Restraint = ' + str(pD['Restraint'])

    return theTitle

def savePlot(pD, fileName, groupRuns = False):

    savePath = my.makeWaterPaths(pD, groupRuns)['pictures']
    if not os.path.exists(savePath):
        os.makedirs(savePath)

    print('Saving ' + fileName + '.pdf' + ' to ' + savePath)
    plt.savefig(savePath + fileName + '.pdf')

def plotFlowRateTemperature():

    kB = 1.381e-23
    mW = 2.991e-26

    T = np.linspace(0, 400, 1000)
    v = np.sqrt(3*kB*T/mW);
    
    plt.plot(T, 10*v)
    plt.xlabel('Temperature (K)')
    plt.ylabel('Corresponding Velocity (A/ns)')
    plt.grid()
    plt.savefig('VofT.pdf')