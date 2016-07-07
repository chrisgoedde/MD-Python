# _*_ coding: utf-8

import os

# Simple function to create a quoted string to protect file/path names with spaces
def quoted(theString):
    """quoted puts double quotes around any string"""

    return '"' + theString + '"'

def config(thePath):
    """config turns a path into the appropriate path for config files"""

    return thePath + 'Config Files/'

def configFile(thePath, theFile):
    """configFile takes a path and a file name and constructs the full path"""

    return config(thePath) + theFile
    
def makePaths(fileName = 'Run-1', type = 'Data', N0 = 20, S = 0, n = 5, m = 5, temp = 300, damping = 1, \
        force = 0, restraint = 0, duration = 10000, minDuration = 10000, \
        dt = 1, outputFreq = 1000, PME = 'off'):

    paths = {}
    paths['home'] = os.path.abspath('..') + '/'
    paths['type'] = paths['home'] + type + '/' + '(' + str(n) + ', ' + str(m) + ')/'
    paths['N0'] = paths['type'] + 'N0 = ' + str(N0) + '/'
    paths['pbc'] = paths['N0'] + 'PBC/'
    paths['solvate'] = paths['pbc'] + 'PME ' + PME + '/' + 'S = ' + str(S) + '/'
    paths['restraint'] = paths['solvate'] + 'R = ' + str(restraint) + '/'
    paths['forcing'] = paths['restraint'] + 'F = ' + str(force) + '/'
    paths['temperature'] = paths['forcing'] + 'Temp = ' + str(temp) \
        + ', Damping = ' + str(damping) + '/'
    paths['tfinal'] = paths['temperature'] + 'Run = ' + str(duration) \
        + ', Min = ' + str(minDuration) + ', dt = ' + str(dt) \
        + ', Out = ' + str(outputFreq) + '/'
    paths['data'] = paths['tfinal'] + fileName + '/'

    return paths
    
def num2Str(theNum):

    if theNum >= 0:
    
        return str(theNum)
        
    else:
    
        return u'−' + str(-theNum)

