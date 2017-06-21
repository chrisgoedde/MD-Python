# _*_ coding: utf-8

import os
import subprocess

def quoted(theString):
    """quoted puts double quotes around any string"""

    return '"' + theString + '"'

def config(thePath):
    """config turns a path into the appropriate path for config files"""

    return thePath + 'Config Files/'

def configFile(thePath, theFile):
    """configFile takes a path and a file name and constructs the full path"""

    return config(thePath) + theFile
    
def makeWaterPaths(paramDict, groupRuns = False):

    paths = {}
    paths['home'] = os.path.abspath('..') + '/' + 'Water' + '/'
    paths['type'] = paths['home'] + paramDict['Top Folder'] + '/' + '(' + str(paramDict['n']) + ', ' + str(paramDict['m']) + ')/'
    paths['N0'] = paths['type'] + 'N0 = ' + str(paramDict['N0']) + '/'
    paths['pbc'] = paths['N0'] + 'PBC/'
    paths['solvate'] = paths['pbc'] + 'PME ' + paramDict['PME'] + '/' + 'S = ' + str(paramDict['S']) + '/'
    paths['restraint'] = paths['solvate'] + 'R = ' + str(paramDict['Restraint']) + '/'
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
        paths['data'] = paths['tfinal'] + paramDict['Data Folder'] + '/'
        paths['minimization'] = config(paths['restraint']) + 'Min = ' + str(paramDict['Min Duration']) \
            +', dt = ' + str(paramDict['dt (fs)']) \
            + ', Out = ' + str(paramDict['outputFreq']) +  '/'
    else:
        paths['data'] = config(paths['restraint']) + 'Min = ' + str(paramDict['Min Duration']) \
            +', dt = ' + str(paramDict['dt (fs)']) \
            + ', Out = ' + str(paramDict['outputFreq']) +  '/'

    if not groupRuns:
        paths['pictures'] = paths['data'] + 'Pictures' + '/'
    else:
        paths['pictures'] = paths['tfinal'] + 'Pictures' + '/'
        
    paths['param'] = os.path.abspath('..') + '/Source/Parameter Files' + '/'

    return paths
    
def makePolymerPaths(paramDict):

    paths = {}
    paths['home'] = os.path.abspath('..') + '/' + 'Polymers' + '/'
    paths['type'] = paths['home'] + paramDict['Top Folder'] + '/' + '(' + str(paramDict['n']) + ', ' + str(paramDict['m']) + ')/'
    paths['L'] = paths['type'] + 'L = ' + str(paramDict['L']) + '/'
    paths['polymer'] = paths['L'] + paramDict['Polymer PDB'] + '/'

    paths['restraint'] = paths['polymer'] + 'R = ' + str(paramDict['Restraint']) + '/'
    paths['driving amplitude'] = paths['restraint'] + 'A = ' + str(paramDict['Driving Amplitude (A)']) + ' A/'
    paths['driving period'] = paths['driving amplitude'] + 'P = ' + str(paramDict['Driving Period (ps)']) + ' ps/'
    paths['tfinal'] = paths['driving period'] + 'Run = ' + str(paramDict['Duration']) \
        + ', dt = ' + str(paramDict['dt (fs)']) \
        + ', Out = ' + str(paramDict['outputFreq']) + '/'

    if paramDict['Run Type'] == 'New':
        paths['data'] = paths['tfinal'] + paramDict['Data Folder'] + '/'
        paths['minimization'] = config(paths['restraint']) + 'Min = ' + str(paramDict['Min Duration']) \
            +', dt = ' + str(paramDict['dt (fs)']) \
            + ', Out = ' + str(paramDict['outputFreq']) +  '/'
    else:
        paths['data'] = config(paths['restraint']) + 'Min = ' + str(paramDict['Min Duration']) \
            +', dt = ' + str(paramDict['dt (fs)']) \
            + ', Out = ' + str(paramDict['outputFreq']) +  '/'
    
    paths['param'] = os.path.abspath('..') + '/Source/Parameter Files' + '/'
    paths['pdb'] = os.path.abspath('..') + '/Source/Polymer PDB Files' + '/'

    return paths
    
def num2Str(theNum):

    if theNum >= 0:
    
        return str(theNum)
        
    else:
    
        return u'−' + str(-theNum)

def makeFlowRateString(mean, std):

    return u'Flow rate = {0:.0f} ± {1:.0f}'.format(mean, std) + ' A/ns'
        
def setWaterParamDefaults():

    paramDict = {}
    paramDict['System'] = 'Water'
    paramDict['Data Folder'] = 'Run-1'
    paramDict['Top Folder'] = 'Data'
    paramDict['Run Type'] = 'Minimize' # Alternate: 'New'
    paramDict['N0'] = 200
    paramDict['S'] = -1
    paramDict['n'] = 4
    paramDict['m'] = 4
    paramDict['Temperature (K)'] = 5
    paramDict['Damping'] = 1
    paramDict['Thermostat'] = 'On'
    paramDict['Force (pN)'] = 0.02
    paramDict['PME'] = 'on'
    paramDict['Restraint'] = 600
    paramDict['Duration'] = 1000000
    paramDict['Min Duration'] = 10000
    paramDict['Preheat Duration'] = 10000
    paramDict['dt (fs)'] = 1
    paramDict['outputFreq'] = 1000
    paramDict['Config File Name'] = 'Config'
    paramDict['Data File Name'] = 'Data'

    return paramDict
    
def setPolymerParamDefaults():

    paramDict = {}
    paramDict['System'] = 'Polymers'
    paramDict['Data Folder'] = 'Run-1'
    paramDict['Top Folder'] = 'Data'
    paramDict['Run Type'] = 'Minimize'
    paramDict['L'] = 200
    paramDict['n'] = 8
    paramDict['m'] = 8
    paramDict['Polymer PSF'] = 'PE50'
    paramDict['Polymer PDB'] = 'PE50m'
    paramDict['Temperature (K)'] = 0
    paramDict['Damping'] = 1
    paramDict['Thermostat'] = 'Off'
    paramDict['Restraint'] = 600
    paramDict['Duration'] = 20000
    paramDict['Min Duration'] = 1000
    paramDict['dt (fs)'] = 1
    paramDict['outputFreq'] = 10
    paramDict['Config File Name'] = 'Config'
    paramDict['Data File Name'] = 'Data'
    paramDict['tcl Forces'] = True
    paramDict['tcl Script'] = 'forces.tcl'
    paramDict['Driving Amplitude (A)'] = 2
    paramDict['Driving Period (ps)'] = 2

    return paramDict

import csv
import pprint 

def readWaterCSV(cvsName, inputParams):

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

def pdbWrite(folder, fileName, package, commands):
    """ pdbWrite creates pdb and/or psf files by calling VMD.
    
    This function is called from nanomer.run() or nanogen.run() and is not meant to be
    called directly.
    
    It writes the files in the appropriate place in the data hierarchy."""

    if not needConfigFiles(folder, fileName):
        return
        
    # Open up the log file for the VMD run.

    log = open(configFile(folder, fileName + ".log"), "w")

    # Feed everything to driveVMD.
    
    driveVMD(package, commands, logFile = log)

    # If we get here we've "succeeded", even if VMD threw an error.
    
    print("################################################################\n" \
        + "Finished writing configuration file(s) in\n" + config(folder) + ".\n" \
        + "################################################################\n")

def driveVMD(package, commands, logFile = subprocess.DEVNULL):
    """ driveVMD feeds tcl commands to VMD to manipulate pdb and/or psf files."""
    
    # Open a pipe to VMD in the shell.
    
    VMDin=subprocess.Popen(["vmd", "-dispdev", "none"], stdin=subprocess.PIPE, \
                           stdout=logFile)

    # Feed the required package to VMD.
    # This assumes that "Research/CNT/MD Simulations/Source" has been added to VMD's
    # auto_path variable in the .vmdrc file through
    # lappend auto_path [file join $env(HOME) "Dropbox/Work/Research/CNT/MD Simulations/Source" ]

    VMDin.stdin.write(package.encode())
    
    # Feed the commands, one by one, to VMD.
    
    for com in commands:
    
        VMDin.stdin.write(com.encode())

    # Clean up the pipe.
    
    VMDin.stdin.flush()
    VMDin.stdin.close
    VMDin.communicate()

def needConfigFiles(folder, fileName):

    # If the required pdb file already exists, we don't need to recreate it,
    # so we just return. Otherwise, we will run VMD to generate it.
        
    if os.path.isfile(configFile(folder, fileName + ".pdb")):
        
        print("################################################################\n" \
              + "Configuration file(s) already exist in:\n" + config(folder) + ".\n" \
              + "################################################################\n")
              
        return False

    print("################################################################\n" \
          + "Writing configuration file(s) in:\n" + config(folder) + ".\n" \
          + "################################################################\n")
        
    # Create the folder in the data hierarchy for the pdb and psf files, if needed.
        
    if not os.path.exists(config(folder)):
        os.makedirs(config(folder))
        
    return True
