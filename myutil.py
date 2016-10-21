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
    
def makePaths(paramDict, groupRuns = False):

    paths = {}
    paths['home'] = os.path.abspath('..') + '/'
    paths['type'] = paths['home'] + paramDict['Folder'] + '/' + '(' + str(paramDict['n']) + ', ' + str(paramDict['m']) + ')/'
    paths['N0'] = paths['type'] + 'N0 = ' + str(paramDict['N0']) + '/'
    paths['pbc'] = paths['N0'] + 'PBC/'
    paths['solvate'] = paths['pbc'] + 'PME ' + paramDict['PME'] + '/' + 'S = ' + str(paramDict['S']) + '/'
    paths['restraint'] = paths['solvate'] + 'R = ' + str(paramDict['Restraint']) + '/'
    paths['minimization'] = config(paths['restraint']) + 'Min = ' + str(paramDict['Min Duration']) \
         +', dt = ' + str(paramDict['dt (fs)']) \
        + ', Out = ' + str(paramDict['outputFreq']) +  '/'
    paths['forcing'] = paths['restraint'] + 'F = ' + str(paramDict['Force (pN)']) + ' pN' + '/'
    if paramDict['Thermostat'] == 'On':
        paths['temperature'] = paths['forcing'] + 'Temp = ' + str(paramDict['Temperature (K)']) \
            + ' K, Damping = ' + str(paramDict['Damping']) + '/'
    else:
        paths['temperature'] = paths['forcing'] + 'Temp = ' + str(paramDict['Temperature (K)']) \
            + ' K, Thermostat = ' + paramDict['Thermostat'] + '/'    
    paths['tfinal'] = paths['temperature'] + 'Run = ' + str(paramDict['Duration']) \
        + ', dt = ' + str(paramDict['dt (fs)']) \
        + ', Out = ' + str(paramDict['outputFreq']) + '/'
    if paramDict['Run Type'] == 'New':
        paths['data'] = paths['tfinal'] + paramDict['File Name'] + '/'
    else:
        paths['parent data'] = paths['tfinal'] + paramDict['File Name'] + '/'
        paths['data'] = paths['tfinal'] + paramDict['File Name'] + '/' + 'Extend-1' + '/'

    if not groupRuns:
    
        paths['pictures'] = paths['data'] + 'Pictures' + '/'
        
    else:
    
        paths['pictures'] = paths['tfinal'] + 'Pictures' + '/'

    return paths
    
def num2Str(theNum):

    if theNum >= 0:
    
        return str(theNum)
        
    else:
    
        return u'−' + str(-theNum)

def makeFlowRateString(mean, std):

    return u'Flow rate = {0:.0f} ± {1:.0f}'.format(mean, std) + ' A/ns'
        
def setParamDefaults():

    paramDict = {}
    paramDict['File Name'] = 'Run-1'
    paramDict['Folder'] = 'Data'
    paramDict['Run Type'] = 'New'
    paramDict['N0'] = 200
    paramDict['S'] = -1
    paramDict['n'] = 4
    paramDict['m'] = 4
    paramDict['Temperature (K)'] = 5
    paramDict['Damping'] = 1
    paramDict['Thermostat'] = 'On'
    paramDict['Force (pN)'] = 0.2
    paramDict['PME'] = 'on'
    paramDict['Restraint'] = 600
    paramDict['Duration'] = 1000000
    paramDict['Min Duration'] = 10000
    paramDict['dt (fs)'] = 1
    paramDict['outputFreq'] = 1000
    paramDict['Config File Name'] = 'Config'

    return paramDict
    
import csv
import pprint 

def readCSV(cvsName, inputParams):

    inputItems = inputParams.items()
    
    pp = pprint.PrettyPrinter(indent=2)
    runList = [];
    with open(cvsName, 'r') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
                
            convert(row, 'Damping')
            convert(row, 'Duration')
            convert(row, 'Min Duration')
            convert(row, 'Force (pN)')
            convert(row, 'N0')
            convert(row, 'S')
            convert(row, 'm')
            convert(row, 'n')
            convert(row, 'outputFreq')
            convert(row, 'dt (fs)')
            convert(row, 'Temperature (K)')
            convert(row, 'Restraint')
            convert(row, 'Run Time (h)')

            allItems = row.items()
            if all(m in allItems for m in inputItems):
                runList.append(row)

    return runList

def convert(theDict, theKey):

    if '.' in theDict[theKey]:
        theDict[theKey] = float(theDict[theKey])
    else:
        theDict[theKey] = int(theDict[theKey])
